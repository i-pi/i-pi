"""Classes that deal with constant pressure simulations.

Contains the algorithms which propagate the position and momenta steps in the
constant pressure ensemble. Holds the properties directly related to
these ensembles, such as the internal and external pressure and stress.

Note that this file also contains a 'BaroMHT' class, that follows more closely the
Martyna, Hughes, Tuckerman implementation of a PIMD barostat. However it is so
close to the BZP implementation that we disabled it for the sake of simplicity.
The original reference is:
G. Martyna, A. Hughes and M. Tuckerman, J. Chem. Phys., 110, 3275.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np

from ipi.utils.depend import *
from ipi.utils.units import *
from ipi.utils.mathtools import (
    matrix_exp,
    sinch,
    mat_taylor,
)
from ipi.engine.thermostats import Thermostat
from ipi.engine.cell import Cell


__all__ = ["Barostat", "BaroBZP", "BaroRGB", "BaroSCBZP", "BaroMTK"]


def mask_from_fix(fix):
    """Builds a mask for the cell velocities based on a list of strings
    that indicates the entries that should be held fixed"""

    m2f = dict(
        xx=(0, 0),
        xy=(0, 1),
        xz=(0, 2),
        yx=(0, 1),
        yy=(1, 1),
        # NB: we re-map to the upper-triangular section, which is the only meaningful one
        yz=(1, 2),
        zx=(0, 2),
        zy=(1, 2),
        zz=(2, 2),
        offdiagonal=[(0, 1), (0, 2), (1, 2)],
    )

    hmask = np.ones(shape=(3, 3))
    # hmask[(1,0)] = 0; hmask[(2,0)] = 0; hmask[(2,1)] = 0
    for f in fix:
        if type(m2f[f]) is list:
            for i in m2f[f]:
                hmask[i] = 0
        else:
            hmask[m2f[f]] = 0
    return hmask


class Barostat(dobject):
    """Base barostat class.

    Gives the standard methods and attributes needed in all the barostat classes.

    Attributes:
       beads: A beads object giving the atoms positions
       cell: A cell object giving the system box.
       forces: A forces object giving the virial and the forces acting on
          each bead.
       nm: An object to do the normal mode transformation.
       thermostat: A thermostat coupled to the barostat degrees of freedom.
       mdof: The number of atomic degrees of freedom

    Depend objects:
       dt: The time step used in the algorithms. Depends on the simulation dt.
       temp: The (classical) simulation temperature. Higher than the physical
          temperature by a factor of the number of beads.
       tau: The timescale associated with the piston
       pext: The external pressure
       stressext: The external stress
       ebaro: The conserved quantity associated with the barostat.
       pot: The potential energy associated with the barostat.
       kstress: The system kinetic stress tensor.
       stress: The system stress tensor.
       press: The system pressure.
    """

    def __init__(
        self, dt=None, temp=None, tau=None, ebaro=None, thermostat=None, nmts=None
    ):
        """Initialises base barostat class.

        Note that the external stress and the external pressure are synchronized.
        This makes most sense going from the stress to the pressure, but if you
        must go in the other direction the stress is assumed to be isotropic.

        Args:
           dt: Optional float giving the time step for the algorithms. Defaults
              to the simulation dt.
           temp: Optional float giving the temperature for the thermostat.
              Defaults to the simulation temp.
           tau: Optional float giving the time scale associated with the barostat.
           ebaro: Optional float giving the conserved quantity already stored
              in the barostat initially. Used on restart.
           thermostat: The thermostat connected to the barostat degree of freedom.
        """

        dself = dd(self)
        dself.dt = depend_value(name="dt")
        if dt is not None:
            self.dt = dt
        else:
            self.dt = 1.0

        dself.temp = depend_value(name="temp")
        if temp is not None:
            self.temp = temp
        else:
            self.temp = 1.0

        dself.tau = depend_value(name="tau")
        if tau is not None:
            self.tau = tau
        else:
            self.tau = 1.0

        dself.ebaro = depend_value(name="ebaro")
        if ebaro is not None:
            self.ebaro = ebaro
        else:
            self.ebaro = 0.0

        if thermostat is None:
            thermostat = Thermostat()
        self.thermostat = thermostat

        # temperature to the thermostat
        dpipe(dself.temp, dd(self.thermostat).temp)
        dself.pext = depend_value(name="pext", value=-1.0)
        dself.stressext = depend_array(name="stressext", value=-np.ones((3, 3), float))

    def bind(self, beads, nm, cell, forces, bias=None, prng=None, fixdof=None, nmts=1):
        """Binds beads, cell and forces to the barostat.

        This takes a beads object, a cell object and a forcefield object and
        makes them members of the barostat. It also then creates the objects that
        will hold the data needed in the barostat algorithms and the dependency
        network.

        Args:
           beads: The beads object from which the bead positions are taken.
           nm: The normal modes propagator object
           cell: The cell object from which the system box is taken.
           forces: The forcefield object from which the force and virial are
              taken.
           prng: The parent PRNG to bind the thermostat to
           fixdof: The number of blocked degrees of freedom.
        """

        self.beads = beads
        self.cell = cell
        self.forces = forces
        self.bias = bias
        self.nm = nm

        dself = dd(self)

        dself.kstress = depend_value(
            name="kstress",
            func=self.get_kstress,
            dependencies=[dd(beads).q, dd(beads).qc, dd(beads).pc, dd(forces).f],
        )
        dself.stress = depend_value(
            name="stress",
            func=self.get_stress,
            dependencies=[dself.kstress, dd(cell).V, dd(forces).vir],
        )

        dself.pot = depend_value(name="pot", value=0.0)

        dself.kin = depend_value(name="kin", value=0.0)

        dself.cell_jacobian = depend_value(name="kin", value=0.0)

        # Stress depend objects for Suzuki-Chin PIMD
        dself.kstress_sc = depend_value(
            name="kstress_sc",
            func=self.get_kstress_sc,
            dependencies=[
                dd(beads).q,
                dd(beads).qc,
                dd(forces).fsc_part_2,
                dd(forces).f,
            ],
        )

        dself.stress_sc = depend_value(
            name="stress_sc",
            func=self.get_stress_sc,
            dependencies=[
                dself.kstress_sc,
                dd(self.cell).V,
                dd(forces).vir,
                dd(forces).virssc_part_2,
            ],
        )

        if fixdof is None:
            self.mdof = float(self.beads.natoms) * 3.0
        else:
            self.mdof = float(self.beads.natoms) * 3.0 - float(fixdof)

        # creates and connects timesteps for the different parts of the propagator
        self.nmtslevels = nmts
        dself.qdt = depend_value(name="qdt", value=self.dt)
        dself.pdt = depend_array(name="pdt", value=np.zeros(nmts, float))
        dself.tdt = depend_value(name="tdt", value=self.dt)
        dpipe(dself.tdt, dd(self.thermostat).dt)

    # THESE SHOULD NOT BE USED ANYMORE
    def get_kstress(self):
        """Calculates the quantum centroid virial kinetic stress tensor
        estimator.
        """

        kst = np.zeros((3, 3), float)
        q = dstrip(self.beads.q)
        qc = dstrip(self.beads.qc)
        pc = dstrip(self.beads.pc)
        m = dstrip(self.beads.m)
        na3 = 3 * self.beads.natoms
        fall = dstrip(self.forces.f)
        if self.bias is None:
            ball = fall * 0.00
        else:
            ball = dstrip(self.bias.f)

        for b in range(self.beads.nbeads):
            for i in range(3):
                for j in range(i, 3):
                    kst[i, j] -= np.dot(
                        q[b, i:na3:3] - qc[i:na3:3], fall[b, j:na3:3] + ball[b, j:na3:3]
                    )

        # NOTE: In order to have a well-defined conserved quantity, the Nf kT term in the
        # diagonal stress estimator must be taken from the centroid kinetic energy.
        for i in range(3):
            kst[i, i] += np.dot(pc[i:na3:3], pc[i:na3:3] / m) * self.beads.nbeads

        return kst

    def kstress_mts_sc(self, level):
        """Calculates the Suzuki-Chin quantum centroid virial kinetic stress tensor
        associated with the forces at a MTS level.
        """

        kst = np.zeros((3, 3), float)
        q = dstrip(self.beads.q)
        qc = dstrip(self.beads.qc)
        pc = dstrip(self.beads.pc)
        m = dstrip(self.beads.m)
        na3 = 3 * self.beads.natoms
        fall = dstrip(self.forces.forces_mts(level)) * (1 + self.forces.coeffsc_part_1)
        if self.bias is None or level != 0:
            ball = fall * 0.00
        else:
            ball = dstrip(self.bias.f)

        for b in range(self.beads.nbeads):
            for i in range(3):
                for j in range(i, 3):
                    kst[i, j] -= np.dot(
                        q[b, i:na3:3] - qc[i:na3:3], fall[b, j:na3:3] + ball[b, j:na3:3]
                    )

        if level == self.nmtslevels - 1:
            for i in range(3):
                kst[i, i] += np.dot(pc[i:na3:3], pc[i:na3:3] / m) * self.beads.nbeads

        return kst

    def kstress_mts(self, level):
        """Calculates the quantum centroid virial kinetic stress tensor
        associated with the forces at a MTS level.
        """

        kst = np.zeros((3, 3), float)
        q = dstrip(self.beads.q)
        qc = dstrip(self.beads.qc)
        pc = dstrip(self.beads.pc)
        m = dstrip(self.beads.m)
        na3 = 3 * self.beads.natoms
        fall = dstrip(self.forces.forces_mts(level))
        if self.bias is None or level != 0:
            ball = fall * 0.00
        else:
            ball = dstrip(self.bias.f)

        for b in range(self.beads.nbeads):
            for i in range(3):
                for j in range(i, 3):
                    kst[i, j] -= np.dot(
                        q[b, i:na3:3] - qc[i:na3:3], fall[b, j:na3:3] + ball[b, j:na3:3]
                    )

        if level == self.nmtslevels - 1:
            for i in range(3):
                kst[i, i] += np.dot(pc[i:na3:3], pc[i:na3:3] / m) * self.beads.nbeads

        return kst

    def get_kstress_sc(self):
        """Calculates the high order part of the Suzuki-Chin
        quantum centroid virial kinetic stress tensor
        associated with the forces at a MTS level.
        """

        kst = np.zeros((3, 3), float)
        q = dstrip(self.beads.q)
        qc = dstrip(self.beads.qc)
        # pc = dstrip(self.beads.pc)
        # m = dstrip(self.beads.m)
        na3 = 3 * self.beads.natoms
        fall = dstrip(self.forces.fsc_part_2)

        for b in range(self.beads.nbeads):
            for i in range(3):
                for j in range(i, 3):
                    kst[i, j] -= np.dot(q[b, i:na3:3] - qc[i:na3:3], fall[b, j:na3:3])
        return kst

    def get_stress(self):
        """Calculates the internal stress tensor."""

        bvir = np.zeros((3, 3), float)
        if self.bias is not None:
            bvir[:] = self.bias.vir

        return (self.kstress + self.forces.virs + bvir) / self.cell.V

    def get_stress_sc(self):
        """Calculates the high order part of the Suzuki-Chin internal stress tensor."""
        return (
            self.kstress_sc + np.sum(dstrip(self.forces.virssc_part_2), axis=0)
        ) / self.cell.V

    def stress_mts_sc(self, level):
        """Calculates the internal Suzuki-Chin stress tensor
        associated with the forces at a MTS level.
        """

        bvir = np.zeros((3, 3), float)
        if self.bias is not None and level == 0:
            bvir[:] = self.bias.vir

        return (
            self.kstress_mts_sc(level)
            + np.sum(
                self.forces.virs_mts(level)
                * (1 + self.forces.coeffsc_part_1).reshape((self.beads.nbeads, 1, 1)),
                axis=0,
            )
            + bvir
        ) / self.cell.V

    def stress_mts(self, level):
        """Calculates the internal stress tensor
        associated with the forces at a MTS level.
        """

        bvir = np.zeros((3, 3), float)
        if self.bias is not None and level == 0:
            bvir[:] = self.bias.vir

        return (
            self.kstress_mts(level) + self.forces.vir_mts(level) + bvir
        ) / self.cell.V

    def pstep(self, level=0):
        """Dummy momenta propagator step."""

        pass

    def qcstep(self, level=0):
        """Dummy centroid position propagator step."""

        pass


class BaroBZP(Barostat):
    """Bussi-Zykova-Parrinello barostat class.

    Just extends the standard class adding finite-dt propagators for the
    barostat velocities, positions, piston.

    Generates dynamics with a stochastic barostat. Implementation details:
    Ceriotti, More, Manolopoulos, Comp. Phys. Comm. 185, 1019, (2013)

    Depend objects:
       p: The momentum associated with the volume degree of freedom.
       m: The mass associated with the volume degree of freedom.
    """

    def __init__(
        self,
        dt=None,
        temp=None,
        tau=None,
        ebaro=None,
        thermostat=None,
        pext=None,
        p=None,
    ):
        """Initializes BZP barostat.

        Args:
           dt: Optional float giving the time step for the algorithms. Defaults
              to the simulation dt.
           temp: Optional float giving the temperature for the thermostat.
              Defaults to the simulation temp.
           pext: Optional float giving the external pressure.
           tau: Optional float giving the time scale associated with the barostat.
           ebaro: Optional float giving the conserved quantity already stored
              in the barostat initially. Used on restart.
           thermostat: The thermostat connected to the barostat degree of freedom.
           p: Optional initial volume conjugate momentum. Defaults to 0.
        """

        super(BaroBZP, self).__init__(dt, temp, tau, ebaro, thermostat)

        dself = dd(self)  # direct access
        dself.p = depend_array(name="p", value=np.atleast_1d(0.0))

        if p is not None:
            self.p = np.asarray([p])
        else:
            self.p = 0.0

        if pext is not None:
            self.pext = pext
        else:
            self.pext = -1.0

    def bind(self, beads, nm, cell, forces, bias=None, prng=None, fixdof=None, nmts=1):
        """Binds beads, cell and forces to the barostat.

        This takes a beads object, a cell object and a forcefield object and
        makes them members of the barostat. It also then creates the objects that
        will hold the data needed in the barostat algorithms and the dependency
        network.

        Args:
           beads: The beads object from which the bead positions are taken.
           nm: The normal modes propagator object
           cell: The cell object from which the system box is taken.
           forces: The forcefield object from which the force and virial are
              taken.
           prng: The parent PRNG to bind the thermostat to
           fixdof: The number of blocked degrees of freedom.
        """

        super(BaroBZP, self).bind(beads, nm, cell, forces, bias, prng, fixdof, nmts)
        dself = dd(self)

        # obtain the thermostat mass from the given time constant
        # note that the barostat temperature is nbeads times the physical T
        dself.m = depend_array(
            name="m",
            value=np.atleast_1d(0.0),
            func=(
                lambda: np.asarray(
                    [self.tau ** 2 * 3 * self.beads.natoms * Constants.kb * self.temp]
                )
            ),
            dependencies=[dself.tau, dself.temp],
        )

        # binds the thermostat to the piston degrees of freedom
        self.thermostat.bind(pm=[self.p, self.m], prng=prng)

        # barostat elastic energy
        dself.pot = depend_value(
            name="pot", func=self.get_pot, dependencies=[dd(cell).V, dself.pext]
        )

        dself.kin = depend_value(
            name="kin",
            func=(lambda: 0.5 * self.p[0] ** 2 / self.m[0]),
            dependencies=[dself.p, dself.m],
        )

        # defines the term that accounts for the explicit dependence of the volume on the ensemble
        dself.cell_jacobian = depend_value(
            name="cell_jacobian",
            func=self.get_cell_jacobian,
            dependencies=[dd(self.cell).V, dself.temp],
        )

        # the barostat energy must be computed from bits & pieces (overwrite the default)
        dself.ebaro = depend_value(
            name="ebaro",
            func=self.get_ebaro,
            dependencies=[
                dself.kin,
                dself.pot,
                dd(self.cell).V,
                dself.temp,
                dd(self.thermostat).ethermo,
            ],
        )

    def get_pot(self):
        """Calculates the elastic strain energy of the cell."""

        # NOTE: since there are nbeads replicas of the unit cell, the enthalpy contains a nbeads factor
        return self.cell.V * self.pext * self.beads.nbeads

    def get_cell_jacobian(self):
        """Calculates the energy term that accounts for the size of the box"""

        return -1.0 * np.log(self.cell.V) * Constants.kb * self.temp

    def get_ebaro(self):
        """Calculates the barostat conserved quantity."""

        return self.thermostat.ethermo + self.kin + self.pot + self.cell_jacobian

    def pstep(self, level=0):
        """Propagates the momentum of the barostat."""

        # we are assuming then that p the coupling between p^2 and dp/dt only involves the fast force
        dt = self.pdt[
            level
        ]  # this is already set to be half a time step at the specified MTS depth
        dt2 = dt ** 2
        dt3 = dt ** 3 / 3.0

        # computes the pressure associated with the forces at each MTS level.
        press = np.trace(self.stress_mts(level)) / 3.0
        self.p += dt * 3.0 * (self.cell.V * press)

        # integerates the kinetic part of the pressure with the force at the inner-most level.
        if level == self.nmtslevels - 1:
            press = 0
            self.p += (
                dt
                * 3.0
                * (
                    self.cell.V * (press - self.beads.nbeads * self.pext)
                    + Constants.kb * self.temp
                )
            )

            pc = dstrip(self.beads.pc)
            fc = (
                np.sum(dstrip(self.forces.forces_mts(level)), axis=0)
                / self.beads.nbeads
            )
            m = dstrip(self.beads.m3)[0]

            self.p += (
                dt2 * np.dot(pc, fc / m) + dt3 * np.dot(fc, fc / m)
            ) * self.beads.nbeads

    def qcstep(self):
        """Propagates the centroid position and momentum and the volume."""

        v = self.p[0] / self.m[0]
        halfdt = (
            self.qdt
        )  # this is set to half the inner loop in all integrators that use a barostat
        expq, expp = (np.exp(v * halfdt), np.exp(-v * halfdt))

        m = dstrip(self.beads.m3)[0]

        self.nm.qnm[0, :] *= expq
        self.nm.qnm[0, :] += ((expq - expp) / (2.0 * v)) * (
            dstrip(self.nm.pnm)[0, :] / m
        )
        self.nm.pnm[0, :] *= expp

        self.cell.h *= expq


class BaroSCBZP(Barostat):
    """The Suzuki Chin Bussi-Zykova-Parrinello barostat class.

    Just extends the standard class adding finite-dt propagators for the
    barostat velocities, positions, piston.

    Generates dynamics with a stochastic barostat. Implementation details:
    Ceriotti, More, Manolopoulos, Comp. Phys. Comm. 185, 1019, (2013)

    Depend objects:
       p: The momentum associated with the volume degree of freedom.
       m: The mass associated with the volume degree of freedom.
    """

    def __init__(
        self,
        dt=None,
        temp=None,
        tau=None,
        ebaro=None,
        thermostat=None,
        pext=None,
        p=None,
    ):
        """Initializes SC barostat.

        Args:
           dt: Optional float giving the time step for the algorithms. Defaults
              to the simulation dt.
           temp: Optional float giving the temperature for the thermostat.
              Defaults to the simulation temp.
           pext: Optional float giving the external pressure.
           tau: Optional float giving the time scale associated with the barostat.
           ebaro: Optional float giving the conserved quantity already stored
              in the barostat initially. Used on restart.
           thermostat: The thermostat connected to the barostat degree of freedom.
           p: Optional initial volume conjugate momentum. Defaults to 0.
        """

        super(BaroSCBZP, self).__init__(dt, temp, tau, ebaro, thermostat)
        dself = dd(self)

        dself.p = depend_array(name="p", value=np.atleast_1d(0.0))

        if p is not None:
            self.p = np.asarray([p])
        else:
            self.p = 0.0

        if pext is not None:
            self.pext = pext
        else:
            self.pext = -1.0

    def bind(self, beads, nm, cell, forces, bias=None, prng=None, fixdof=None, nmts=1):
        """Binds beads, cell and forces to the barostat.

        This takes a beads object, a cell object and a forcefield object and
        makes them members of the barostat. It also then creates the objects that
        will hold the data needed in the barostat algorithms and the dependency
        network.

        Args:
           beads: The beads object from which the bead positions are taken.
           nm: The normal modes propagator object
           cell: The cell object from which the system box is taken.
           forces: The forcefield object from which the force and virial are
              taken.
           prng: The parent PRNG to bind the thermostat to
           fixdof: The number of blocked degrees of freedom.
        """

        super(BaroSCBZP, self).bind(beads, nm, cell, forces, bias, prng, fixdof, nmts)
        dself = dd(self)

        # obtain the thermostat mass from the given time constant
        dself.m = depend_array(
            name="m",
            value=np.atleast_1d(0.0),
            func=(
                lambda: np.asarray(
                    [self.tau ** 2 * 3 * self.beads.natoms * Constants.kb * self.temp]
                )
            ),
            dependencies=[dself.tau, dself.temp],
        )

        # binds the thermostat to the piston degrees of freedom
        self.thermostat.bind(pm=[self.p, self.m], prng=prng)

        # barostat elastic energy
        dself.pot = depend_value(
            name="pot", func=self.get_pot, dependencies=[dd(cell).V, dself.pext]
        )

        dself.kin = depend_value(
            name="kin",
            func=(lambda: 0.5 * self.p[0] ** 2 / self.m[0]),
            dependencies=[dself.p, dself.m],
        )

        # defines the term that accounts for the explicit dependence of the ensemble probability on the volume
        dself.cell_jacobian = depend_value(
            name="cell_jacobian",
            func=self.get_cell_jacobian,
            dependencies=[dd(self.cell).V, dself.temp],
        )

        # the barostat energy must be computed from bits & pieces (overwrite the default)
        dself.ebaro = depend_value(
            name="ebaro",
            func=self.get_ebaro,
            dependencies=[
                dself.kin,
                dself.pot,
                dself.cell_jacobian,
                dd(self.thermostat).ethermo,
            ],
        )

    def get_pot(self):
        """Calculates the elastic strain energy of the cell."""

        # NOTE: since there are nbeads replicas of the unit cell, the enthalpy contains a nbeads factor
        return self.cell.V * self.pext * self.beads.nbeads

    def get_cell_jacobian(self):
        """Calculates the energy term that accounts for the size of the box"""

        return -1.0 * np.log(self.cell.V) * Constants.kb * self.temp

    def get_ebaro(self):
        """Calculates the barostat conserved quantity."""

        return self.thermostat.ethermo + self.kin + self.pot + self.cell_jacobian

    def pscstep(self):
        """Propagates the momentum of the barostat with respect to the
        high order part of the Suzuki-Chin stress"""

        # integrates with respect to the "high order" part of the stress with a timestep of dt /2
        press = np.trace(self.stress_sc) / 3.0
        self.p += self.dt / 2 * 3.0 * (self.cell.V * press)

    def pstep(self, level=0):
        """Propagates the momentum of the barostat."""

        # self.pscstep()

        # we are assuming then that p the coupling between p^2 and dp/dt only involves the fast force
        dt = self.pdt[
            level
        ]  # this is already set to be half a time step at the specified MTS depth
        dt2 = dt ** 2
        dt3 = dt ** 3 / 3.0

        # computes the pressure associated with the forces at each MTS level and adds the +- 1/3 SC correction.
        press = np.trace(self.stress_mts_sc(level)) / 3.0
        self.p += dt * 3.0 * (self.cell.V * press)

        # integerates the kinetic part of the pressure with the force at the inner-most level.
        if level == self.nmtslevels - 1:
            press = 0
            self.p += (
                dt
                * 3.0
                * (
                    self.cell.V * (press - self.beads.nbeads * self.pext)
                    + Constants.kb * self.temp
                )
            )
            pc = dstrip(self.beads.pc)
            fc = (
                np.sum(
                    dstrip(self.forces.forces_mts(level))
                    * (1 + self.forces.coeffsc_part_1),
                    axis=0,
                )
                / self.beads.nbeads
            )
            m = dstrip(self.beads.m3)[0]

            self.p += (
                dt2 * np.dot(pc, fc / m) + dt3 * np.dot(fc, fc / m)
            ) * self.beads.nbeads

    def qcstep(self):
        """Propagates the centroid position and momentum and the volume."""

        v = self.p[0] / self.m[0]
        halfdt = (
            self.qdt
        )  # this is set to half the inner loop in all integrators that use a barostat
        expq, expp = (np.exp(v * halfdt), np.exp(-v * halfdt))

        m = dstrip(self.beads.m3)[0]

        self.nm.qnm[0, :] *= expq
        self.nm.qnm[0, :] += ((expq - expp) / (2.0 * v)) * (
            dstrip(self.nm.pnm)[0, :] / m
        )
        self.nm.pnm[0, :] *= expp

        self.cell.h *= expq


class BaroRGB(Barostat):
    """Raiteri-Gale-Bussi constant stress barostat class (JPCM 23, 334213, 2011).

    Just extends the standard class adding finite-dt propagators for the barostat
    velocities, positions, piston.

    Depend objects:
    p: The momentum matrix associated with the cell degrees of freedom.
    m: The mass associated with the cell degree of freedom.
    """

    def __init__(
        self,
        dt=None,
        temp=None,
        tau=None,
        ebaro=None,
        thermostat=None,
        stressext=None,
        h0=None,
        p=None,
        hfix=None,
    ):
        """Initializes RGB barostat.

        Args:
        dt: Optional float giving the time step for the algorithms. Defaults
        to the simulation dt.
        temp: Optional float giving the temperature for the thermostat.
        Defaults to the simulation temp.
        stressext: Optional float giving the external pressure.
        tau: Optional float giving the time scale associated with the barostat.
        ebaro: Optional float giving the conserved quantity already stored
        in the barostat initially. Used on restart.
        thermostat: The thermostat connected to the barostat degree of freedom.
        p: Optional initial volume conjugate momentum. Defaults to 0.
        """

        super(BaroRGB, self).__init__(dt, temp, tau, ebaro, thermostat)

        # non-zero elements of the cell momentum are only
        # pxx pyy pzz pxy pxz pyz, but we want to access it either as a
        # 6-vector or as a 3x3 upper triangular tensor.
        # we use a synchronizer to achieve that

        dself = dd(self)

        sync_baro = synchronizer()
        dself.p6 = depend_array(
            name="p6",
            value=np.zeros(6, float),
            synchro=sync_baro,
            func={"p": self.get_3x3to6},
        )
        dself.p = depend_array(
            name="p",
            value=np.zeros((3, 3), float),
            synchro=sync_baro,
            func={"p6": self.get_6to3x3},
        )

        if p is not None:
            self.p = p
        else:
            self.p = 0.0

        if h0 is not None:
            self.h0 = h0
        else:
            self.h0 = Cell()

        if hfix is None:
            hfix = []
        self.hfix = hfix

        hmask = mask_from_fix(self.hfix)

        # mask to zero out components of the cell velocity, to implement cell-boundary constraints
        dself.hmask = depend_array(name="hmask", value=hmask)
        # number of ones in the UT part of the mask
        self.L = np.diag([hmask[0].sum(), hmask[1, 1:].sum(), hmask[2, 2:].sum()])

        if stressext is not None:
            self.stressext = stressext
        else:
            self.stressext[:] = -1.0

    def bind(self, beads, nm, cell, forces, bias=None, prng=None, fixdof=None, nmts=1):
        """Binds beads, cell and forces to the barostat.

        This takes a beads object, a cell object and a forcefield object and
        makes them members of the barostat. It also then creates the objects that
        will hold the data needed in the barostat algorithms and the dependency
        network.

        Args:
        beads: The beads object from which the bead positions are taken.
        nm: The normal modes propagator object
        cell: The cell object from which the system box is taken.
        forces: The forcefield object from which the force and virial are
        taken.
        prng: The parent PRNG to bind the thermostat to
        fixdof: The number of blocked degrees of freedom.
        """

        super(BaroRGB, self).bind(beads, nm, cell, forces, bias, prng, fixdof, nmts)

        # obtain the thermostat mass from the given time constant (1/3 of what used for the corresponding NPT case)
        # note that the barostat temperature is nbeads times the physical T
        dself = dd(self)

        dself.m = depend_array(
            name="m",
            value=np.atleast_1d(0.0),
            func=(
                lambda: np.asarray(
                    [self.tau ** 2 * self.beads.natoms * Constants.kb * self.temp]
                )
            ),
            dependencies=[dself.tau, dself.temp],
        )

        dself.m6 = depend_array(
            name="m6",
            value=np.zeros(6, float),
            func=(lambda: np.asarray([1, 1, 1, 1, 1, 1]) * self.m[0]),
            dependencies=[dself.m],
        )

        # overrides definition of pot to depend on the many things it depends on for anisotropic cell
        dself.pot = depend_value(
            name="pot",
            func=self.get_pot,
            dependencies=[
                dd(self.cell).h,
                dd(self.h0).h,
                dd(self.h0).V,
                dd(self.h0).ih,
                dself.stressext,
            ],
        )

        # binds the thermostat to the piston degrees of freedom
        self.thermostat.bind(pm=[self.p6, self.m6], prng=prng)

        dself.kin = depend_value(
            name="kin",
            func=(lambda: 0.5 * np.trace(np.dot(self.p.T, self.p)) / self.m[0]),
            dependencies=[dself.p, dself.m],
        )

        # defines the term that accounts for the explicit dependence of the ensemble on the cell
        dself.cell_jacobian = depend_value(
            name="cell_jacobian",
            func=self.get_cell_jacobian,
            dependencies=[dd(self.cell).h, dself.temp],
        )

        # the barostat energy must be computed from bits & pieces (overwrite the default)
        dself.ebaro = depend_value(
            name="ebaro",
            func=self.get_ebaro,
            dependencies=[
                dself.kin,
                dself.pot,
                dself.cell_jacobian,
                dd(self.thermostat).ethermo,
            ],
        )

    def get_3x3to6(self):
        rp = np.zeros(6, float)
        rp[0] = self.p[0, 0]
        rp[1] = self.p[1, 1]
        rp[2] = self.p[2, 2]
        rp[3] = self.p[0, 1]
        rp[4] = self.p[0, 2]
        rp[5] = self.p[1, 2]
        return rp

    def get_6to3x3(self):
        rp = np.zeros((3, 3), float)
        rp[0, 0] = self.p6[0]
        rp[1, 1] = self.p6[1]
        rp[2, 2] = self.p6[2]
        rp[0, 1] = self.p6[3]
        rp[0, 2] = self.p6[4]
        rp[1, 2] = self.p6[5]
        return rp

    def get_pot(self):
        """Calculates the elastic strain energy of the cell."""

        # NOTE: since there are nbeads replicas of the unit cell, the enthalpy contains a nbeads factor
        eps = np.dot(self.cell.h, self.h0.ih)
        eps = np.dot(eps.T, eps)
        eps = 0.5 * (eps - np.identity(3))

        return self.h0.V * np.trace(np.dot(self.stressext, eps)) * self.beads.nbeads
        # p = np.trace(self.stressext)/3
        # return (p*(self.cell.V-self.h0.V) +
        #    self.h0.V * np.trace(np.dot(self.stressext-np.eye(3)*p, eps)) ) * self.beads.nbeads

    def get_cell_jacobian(self):
        """Calculates the energy term that accounts for the size of the box."""

        cj = np.sum([(3 - i) * np.log(self.cell.h[i][i]) for i in range(3)])
        return -1.0 * Constants.kb * self.temp * cj

    def get_ebaro(self):
        """Calculates the barostat conserved quantity."""

        return self.thermostat.ethermo + self.kin + self.pot + self.cell_jacobian

    def pstep(self, level=0):
        """Propagates the momenta for half a time step."""

        dt = self.pdt[level]
        dt2 = dt ** 2
        dt3 = dt ** 3 / 3.0

        hh0 = np.dot(self.cell.h, self.h0.ih)
        pi_ext = np.dot(hh0, np.dot(self.stressext, hh0.T)) * self.h0.V / self.cell.V

        stress = dstrip(self.stress_mts(level))
        self.p += dt * (self.cell.V * np.triu(stress))

        # integerates the kinetic part of the stress with the force at the inner-most level.
        if level == self.nmtslevels - 1:

            self.p += dt * (
                self.cell.V * np.triu(-self.beads.nbeads * pi_ext)
                + Constants.kb * self.temp * self.L
            )

            pc = dstrip(self.beads.pc).reshape(self.beads.natoms, 3)
            fc = (
                np.sum(dstrip(self.forces.forces_mts(level)), axis=0).reshape(
                    self.beads.natoms, 3
                )
                / self.beads.nbeads
            )
            fcTonm = (fc / dstrip(self.beads.m3)[0].reshape(self.beads.natoms, 3)).T

            self.p += (
                np.triu(dt2 * np.dot(fcTonm, pc) + dt3 * np.dot(fcTonm, fc))
                * self.beads.nbeads
            )

        # now apply the mask (and accumulate the associated change in conserved quantity)
        # we use the thermostat conserved quantity accumulator, so we don't need to create a new one
        self.thermostat.ethermo += self.kin
        self.p *= self.hmask
        self.thermostat.ethermo -= self.kin

    def qcstep(self):
        """Propagates the centroid position and momentum and the volume."""

        v = self.p / self.m[0]
        halfdt = self.qdt
        # compute in one go dt sinh v/v to handle zero-velocity cases
        # check eigen vector exists
        if np.count_nonzero(v.diagonal()) == 0:
            # compute sinhv/v by directly taylor, if eigen vector not exists
            sinh = mat_taylor(v * halfdt, function="sinhx/x")
        else:
            eigvals, eigvecs = np.linalg.eig(v)
            ieigvecs = np.linalg.inv(eigvecs)
            sinh = halfdt * np.dot(
                eigvecs, np.dot(np.diag(sinch(halfdt * eigvals)), ieigvecs)
            )

        expq, expp = (matrix_exp(v * halfdt), matrix_exp(-v * halfdt))
        # oldsinh = np.dot(invert_ut3x3(v), (expq - expp) / (2.0))

        q = dstrip(self.nm.qnm)[0].copy().reshape((self.beads.natoms, 3))
        p = dstrip(self.nm.pnm)[0].copy()

        q = np.dot(q, expq.T)
        q += np.dot((p / self.beads.m3[0]).reshape((self.beads.natoms, 3)), sinh.T)
        p = np.dot(p.reshape((self.beads.natoms, 3)), expp.T)

        # now apply the mask (and accumulate the associated change in conserved quantity)
        # we use the thermostat conserved quantity accumulator, so we don't need to create a new one
        self.thermostat.ethermo += self.kin
        self.p *= self.hmask
        self.thermostat.ethermo -= self.kin

        self.nm.qnm[0] = q.reshape((self.beads.natoms * 3))
        self.nm.pnm[0] = p.reshape((self.beads.natoms * 3))

        self.cell.h = np.dot(expq, self.cell.h)


class BaroMTK(Barostat):
    """Martyna-Tobias-Klein flexible cell, constant pressure barostat class (JCP 101, 1994).

    Just extends the standard class adding finite-dt propagators for the barostat
    velocities, positions, piston.

    Depend objects:
    p: The momentum matrix associated with the cell degrees of freedom.
    m: The mass associated with the cell degree of freedom.
    """

    def __init__(
        self,
        dt=None,
        temp=None,
        tau=None,
        ebaro=None,
        thermostat=None,
        pext=None,
        p=None,
        hfix=None,
    ):
        """Initializes RGB barostat.

        Args:
        dt: Optional float giving the time step for the algorithms. Defaults
        to the simulation dt.
        temp: Optional float giving the temperature for the thermostat.
        Defaults to the simulation temp.
        stressext: Optional float giving the external pressure.
        tau: Optional float giving the time scale associated with the barostat.
        ebaro: Optional float giving the conserved quantity already stored
        in the barostat initially. Used on restart.
        thermostat: The thermostat connected to the barostat degree of freedom.
        p: Optional initial volume conjugate momentum. Defaults to 0.
        """

        super(BaroMTK, self).__init__(dt, temp, tau, ebaro, thermostat)

        # non-zero elements of the cell momentum are only
        # pxx pyy pzz pxy pxz pyz, but we want to access it either as a
        # 6-vector or as a 3x3 upper triangular tensor.
        # we use a synchronizer to achieve that

        dself = dd(self)

        sync_baro = synchronizer()
        dself.p6 = depend_array(
            name="p6",
            value=np.zeros(6, float),
            synchro=sync_baro,
            func={"p": self.get_3x3to6},
        )
        dself.p = depend_array(
            name="p",
            value=np.zeros((3, 3), float),
            synchro=sync_baro,
            func={"p6": self.get_6to3x3},
        )

        if p is not None:
            self.p = p
        else:
            self.p = 0.0

        if hfix is None:
            hfix = []
        self.hfix = hfix
        hmask = mask_from_fix(hfix)

        # mask to zero out components of the cell velocity, to implement cell-boundary constraints
        dself.hmask = depend_array(name="hmask", value=hmask)
        # number of ones in the UT part of the mask
        self.L = np.diag([hmask[0].sum(), hmask[1, 1:].sum(), hmask[2, 2:].sum()])

        if pext is not None:
            self.pext = pext
        else:
            self.pext = -1.0

    def bind(self, beads, nm, cell, forces, bias=None, prng=None, fixdof=None, nmts=1):
        """Binds beads, cell and forces to the barostat.

        This takes a beads object, a cell object and a forcefield object and
        makes them members of the barostat. It also then creates the objects that
        will hold the data needed in the barostat algorithms and the dependency
        network.

        Args:
        beads: The beads object from which the bead positions are taken.
        nm: The normal modes propagator object
        cell: The cell object from which the system box is taken.
        forces: The forcefield object from which the force and virial are
        taken.
        prng: The parent PRNG to bind the thermostat to
        fixdof: The number of blocked degrees of freedom.
        """

        super(BaroMTK, self).bind(beads, nm, cell, forces, bias, prng, fixdof, nmts)

        # obtain the thermostat mass from the given time constant (1/3 of what used for the corresponding NPT case)
        # note that the barostat temperature is nbeads times the physical T
        dself = dd(self)

        dself.m = depend_array(
            name="m",
            value=np.atleast_1d(0.0),
            func=(
                lambda: np.asarray(
                    [self.tau ** 2 * self.beads.natoms * Constants.kb * self.temp]
                )
            ),
            dependencies=[dself.tau, dself.temp],
        )

        dself.m6 = depend_array(
            name="m6",
            value=np.zeros(6, float),
            func=(lambda: np.asarray([1, 1, 1, 1, 1, 1]) * self.m[0]),
            dependencies=[dself.m],
        )

        # overrides definition of pot to depend on the many things it depends on for anisotropic cell
        dself.pot = depend_value(
            name="pot", func=self.get_pot, dependencies=[dd(self.cell).h, dself.pext]
        )

        # binds the thermostat to the piston degrees of freedom
        self.thermostat.bind(pm=[self.p6, self.m6], prng=prng)

        dself.kin = depend_value(
            name="kin",
            func=(lambda: 0.5 * np.trace(np.dot(self.p.T, self.p)) / self.m[0]),
            dependencies=[dself.p, dself.m],
        )

        # defines the term that accounts for the explicit dependence of the ensemble on the cell
        dself.cell_jacobian = depend_value(
            name="cell_jacobian",
            func=self.get_cell_jacobian,
            dependencies=[dd(self.cell).h, dself.temp],
        )

        # the barostat energy must be computed from bits & pieces (overwrite the default)
        dself.ebaro = depend_value(
            name="ebaro",
            func=self.get_ebaro,
            dependencies=[
                dself.kin,
                dself.pot,
                dself.cell_jacobian,
                dd(self.thermostat).ethermo,
            ],
        )

    def get_3x3to6(self):
        rp = np.zeros(6, float)
        rp[0] = self.p[0, 0]
        rp[1] = self.p[1, 1]
        rp[2] = self.p[2, 2]
        rp[3] = self.p[0, 1]
        rp[4] = self.p[0, 2]
        rp[5] = self.p[1, 2]
        return rp

    def get_6to3x3(self):
        rp = np.zeros((3, 3), float)
        rp[0, 0] = self.p6[0]
        rp[1, 1] = self.p6[1]
        rp[2, 2] = self.p6[2]
        rp[0, 1] = self.p6[3]
        rp[0, 2] = self.p6[4]
        rp[1, 2] = self.p6[5]
        return rp

    def get_pot(self):
        """Calculates the elastic strain energy of the cell."""

        # return self.h0.V * np.trace(np.dot(self.stressext, eps)) * self.beads.nbeads
        return self.cell.V * self.pext * self.beads.nbeads
        # p = np.trace(self.stressext)/3
        # return (p*(self.cell.V-self.h0.V) +
        #    self.h0.V * np.trace(np.dot(self.stressext-np.eye(3)*p, eps)) ) * self.beads.nbeads

    def get_cell_jacobian(self):
        """Calculates the energy term that accounts for the size of the box."""

        cj = np.sum([(3 - i) * np.log(self.cell.h[i][i]) for i in range(3)])
        return -1.0 * Constants.kb * self.temp * cj

    def get_ebaro(self):
        """Calculates the barostat conserved quantity."""

        return self.thermostat.ethermo + self.kin + self.pot + self.cell_jacobian

    def pstep(self, level=0):
        """Propagates the momenta for half a time step."""

        dt = self.pdt[level]
        dt2 = dt ** 2
        dt3 = dt ** 3 / 3.0

        pi_ext = (
            np.eye(3) * self.pext
        )  # np.dot(hh0, np.dot(self.stressext, hh0.T)) * self.h0.V / self.cell.V
        #        L = np.diag([3, 2, 1])

        stress = dstrip(self.stress_mts(level))
        self.p += dt * (self.cell.V * np.triu(stress))

        # integerates the kinetic part of the stress with the force at the inner-most level.
        if level == self.nmtslevels - 1:

            self.p += dt * (
                self.cell.V * np.triu(-self.beads.nbeads * pi_ext)
                + Constants.kb * self.temp * self.L
            )

            pc = dstrip(self.beads.pc).reshape(self.beads.natoms, 3)
            fc = (
                np.sum(dstrip(self.forces.forces_mts(level)), axis=0).reshape(
                    self.beads.natoms, 3
                )
                / self.beads.nbeads
            )
            fcTonm = (fc / dstrip(self.beads.m3)[0].reshape(self.beads.natoms, 3)).T

            self.p += (
                np.triu(dt2 * np.dot(fcTonm, pc) + dt3 * np.dot(fcTonm, fc))
                * self.beads.nbeads
            )

        # now apply the mask (and accumulate the associated change in conserved quantity)
        # we use the thermostat conserved quantity accumulator, so we don't need to create a new one
        self.thermostat.ethermo += self.kin
        self.p *= self.hmask
        self.thermostat.ethermo -= self.kin

    def qcstep(self):
        """Propagates the centroid position and momentum and the volume."""

        v = self.p / self.m[0]
        halfdt = self.qdt
        # compute in one go dt sinh v/v to handle zero-velocity cases
        # check eigen vector exists
        if np.count_nonzero(v.diagonal()) == 0:
            # compute sinhv/v by directly taylor, if eigen vector not exists
            sinh = mat_taylor(v * halfdt, function="sinhx/x")
        else:
            eigvals, eigvecs = np.linalg.eig(v)
            ieigvecs = np.linalg.inv(eigvecs)
            sinh = halfdt * np.dot(
                eigvecs, np.dot(np.diag(sinch(halfdt * eigvals)), ieigvecs)
            )
        expq, expp = (matrix_exp(v * halfdt), matrix_exp(-v * halfdt))

        #        oldsinh = np.dot(invert_ut3x3(v), (expq - expp) / (2.0))

        q = dstrip(self.nm.qnm)[0].copy().reshape((self.beads.natoms, 3))
        p = dstrip(self.nm.pnm)[0].copy()

        q = np.dot(q, expq.T)
        q += np.dot((p / self.beads.m3[0]).reshape((self.beads.natoms, 3)), sinh.T)
        p = np.dot(p.reshape((self.beads.natoms, 3)), expp.T)

        # now apply the mask (and accumulate the associated change in conserved quantity)
        # we use the thermostat conserved quantity accumulator, so we don't need to create a new one
        self.thermostat.ethermo += self.kin
        self.p *= self.hmask
        self.thermostat.ethermo -= self.kin

        self.nm.qnm[0] = q.reshape((self.beads.natoms * 3))
        self.nm.pnm[0] = p.reshape((self.beads.natoms * 3))

        self.cell.h = np.dot(expq, self.cell.h)

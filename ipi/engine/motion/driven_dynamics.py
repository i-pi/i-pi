"""Contains the classes that deal with the different dynamics required in
different types of ensembles.

Holds the algorithms required for normal mode propagators, and the objects to
do the constant temperature and pressure algorithms. Also calculates the
appropriate conserved energy quantity for the ensemble of choice.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

import numpy as np

from ipi.utils.depend import *
from ipi.engine.motion.eda import EDA
from ipi.utils.units import Constants
from ipi.engine.motion.dynamics import NVEIntegrator, Dynamics


class DrivenDynamics(Dynamics):
    """self (path integral) molecular dynamics class.

    Gives the standard methods and attributes needed in all the
    dynamics classes.

    Attributes:
        beads: A beads object giving the atoms positions.
        cell: A cell object giving the system box.
        forces: A forces object giving the virial and the forces acting on
            each bead.
        prng: A random number generator object.
        nm: An object which does the normal modes transformation.

    Depend objects:
        econs: The conserved energy quantity appropriate to the given
            ensemble. Depends on the various energy terms which make it up,
            which are different depending on the ensemble.he
        temp: The system temperature.
        dt: The timestep for the algorithms.
        ntemp: The simulation temperature. Will be nbeads times higher than
            the system temperature as PIMD calculations are done at this
            effective classical temperature.
    """

    def __init__(
        self,
        efield=None,
        bec=None,
        *argc,
        **argv,
    ):
        """Initialises a "dynamics" motion object.

        Args:
            dt: The timestep of the simulation algorithms.
            fixcom: An optional boolean which decides whether the centre of mass
                motion will be constrained or not. Defaults to False.
        """

        super().__init__(*argc, **argv)

        if self.enstype == "eda-nve":
            # NVE integrator with an external time-dependent driving (electric field)
            self.integrator = EDANVEIntegrator()
        else:
            self.integrator = DummyIntegrator()

        # if the dynamics is driven, allocate necessary objects
        self.efield = efield
        self.bec = bec
        self.eda = EDA(self.efield, self.bec)

    def bind(self, ens, beads, nm, cell, bforce, prng, omaker):

        super().bind(ens, beads, nm, cell, bforce, prng, omaker)

        self.eda.bind(self.ensemble, self.enstype)

        # now that the timesteps are decided, we proceed to bind the integrator.
        self.integrator.bind(self)

        # applies constraints immediately after initialization.
        self.integrator.pconstraints()

    def step(self, step=None):
        """Advances the dynamics by one time step"""

        super().step(step)

        # Check that these variable are the same.
        # If they are not the same, then there is a bug in the code
        dt = abs(self.ensemble.time - self.integrator.mts_time)
        if dt > 1e-12:
            raise ValueError(
                "The time at which the Electric Field is evaluated is not properly updated!"
            )


dproperties(DrivenDynamics, ["dt", "nmts", "splitting", "ntemp"])


class DummyIntegrator:
    """No-op integrator for (PI)MD"""

    def __init__(self):
        pass

    def get_qdt(self):
        return self.dt * 0.5 / self.inmts

    def get_pdt(self):
        dtl = 1.0 / self.nmts
        for i in range(1, len(dtl)):
            dtl[i] *= dtl[i - 1]
        dtl *= self.dt * 0.5
        return dtl

    def get_tdt(self):
        if self.splitting == "obabo":
            return self.dt * 0.5
        elif self.splitting == "baoab":
            return self.dt
        else:
            raise ValueError(
                "Invalid splitting requested. Only OBABO and BAOAB are supported."
            )

    def bind(self, motion):
        """Reference all the variables for simpler access."""

        self.beads = motion.beads
        self.bias = motion.ensemble.bias
        self.ensemble = motion.ensemble
        self.forces = motion.forces
        self.prng = motion.prng
        self.nm = motion.nm
        self.thermostat = motion.thermostat
        self.barostat = motion.barostat
        self.fixcom = motion.fixcom
        self.fixatoms = motion.fixatoms
        self.enstype = motion.enstype
        if motion.enstype in EDA.integrators:
            self.eda = motion.eda

        # no need to dpipe these are really just references
        self._splitting = motion._splitting
        self._dt = motion._dt
        self._nmts = motion._nmts

        # total number of iteration in the inner-most MTS loop
        self._inmts = depend_value(name="inmts", func=lambda: np.prod(self.nmts))
        self._nmtslevels = depend_value(name="nmtslevels", func=lambda: len(self.nmts))
        # these are the time steps to be used for the different parts of the integrator
        self._qdt = depend_value(
            name="qdt",
            func=self.get_qdt,
            dependencies=[self._splitting, self._dt, self._inmts],
        )  # positions
        self._qdt_on_m = depend_array(
            name="qdt_on_m",
            value=np.zeros(3 * self.beads.natoms),
            func=lambda: self.qdt / dstrip(self.beads.m3)[0],
        )
        self._pdt = depend_array(
            name="pdt",
            func=self.get_pdt,
            value=np.zeros(len(self.nmts)),
            dependencies=[self._splitting, self._dt, self._nmts],
        )  # momenta
        self._tdt = depend_value(
            name="tdt",
            func=self.get_tdt,
            dependencies=[self._splitting, self._dt, self._nmts],
        )  # thermostat

        dpipe(self._qdt, self.nm._dt)
        dpipe(self._dt, self.barostat._dt)
        dpipe(self._qdt, self.barostat._qdt)
        dpipe(self._pdt, self.barostat._pdt)
        dpipe(self._tdt, self.barostat._tdt)
        dpipe(self._tdt, self.thermostat._dt)

        if motion.enstype == "sc" or motion.enstype == "scnpt":
            # coefficients to get the (baseline) trotter to sc conversion
            self.coeffsc = np.ones((self.beads.nbeads, 3 * self.beads.natoms), float)
            self.coeffsc[::2] /= -3.0
            self.coeffsc[1::2] /= 3.0

        # check stress tensor
        self._stresscheck = True

    def pstep(self):
        """Dummy momenta propagator which does nothing."""
        pass

    def qcstep(self):
        """Dummy centroid position propagator which does nothing."""
        pass

    def step(self, step=None):
        """Dummy simulation time step which does nothing."""
        pass

    def pconstraints(self):
        """This removes the centre of mass contribution to the kinetic energy.

        Calculates the centre of mass momenta, then removes the mass weighted
        contribution from each atom. If the ensemble defines a thermostat, then
        the contribution to the conserved quantity due to this subtraction is
        added to the thermostat heat energy, as it is assumed that the centre of
        mass motion is due to the thermostat.

        If there is a choice of thermostats, the thermostat
        connected to the centroid is chosen.
        """

        if self.fixcom:
            na3 = self.beads.natoms * 3
            nb = self.beads.nbeads
            p = dstrip(self.beads.p)
            m = dstrip(self.beads.m3)[:, 0:na3:3]
            M = self.beads[0].M
            Mnb = M * nb

            dens = 0
            for i in range(3):
                pcom = p[:, i:na3:3].sum()
                dens += pcom**2
                pcom /= Mnb
                self.beads.p[:, i:na3:3] -= m * pcom

            self.ensemble.eens += dens * 0.5 / Mnb

        if len(self.fixatoms) > 0:
            for bp in self.beads.p:
                m = dstrip(self.beads.m)
                self.ensemble.eens += 0.5 * np.dot(
                    bp[self.fixatoms * 3], bp[self.fixatoms * 3] / m[self.fixatoms]
                )
                self.ensemble.eens += 0.5 * np.dot(
                    bp[self.fixatoms * 3 + 1],
                    bp[self.fixatoms * 3 + 1] / m[self.fixatoms],
                )
                self.ensemble.eens += 0.5 * np.dot(
                    bp[self.fixatoms * 3 + 2],
                    bp[self.fixatoms * 3 + 2] / m[self.fixatoms],
                )
                bp[self.fixatoms * 3] = 0.0
                bp[self.fixatoms * 3 + 1] = 0.0
                bp[self.fixatoms * 3 + 2] = 0.0


dproperties(
    DummyIntegrator,
    ["splitting", "nmts", "dt", "inmts", "nmtslevels", "qdt", "pdt", "tdt", "qdt_on_m"],
)


class EDAIntegrator(DummyIntegrator):
    """Integrator object for simulations using the Electric Dipole Approximation (EDA)
    when an external electric field is applied.
    """

    def __init__(self):
        super().__init__()

    def bind(self, motion):
        """bind variables"""
        super().bind(motion)

        self._time = self.eda._time
        self._mts_time = self.eda._mts_time

        dep = [
            self._time,
            self._mts_time,
            self.eda.Born_Charges._bec,
            self.eda.Electric_Field._Efield,
        ]
        self._EDAforces = depend_array(
            name="EDAforces",
            func=self._eda_forces,
            value=np.zeros((self.beads.nbeads, self.beads.natoms * 3)),
            dependencies=dep,
        )
        pass

    def pstep(self, level=0):
        """Velocity Verlet momentum propagator."""
        self.beads.p += self.EDAforces * self.pdt[level]
        if dstrip(self.mts_time) == dstrip(self.time):
            # it's the first time that 'pstep' is called
            # then we need to update 'mts_time'
            self.mts_time += dstrip(self.dt)
            # the next time this condition will be 'False'
            # so we will avoid to re-compute the EDAforces
        pass

    def _eda_forces(self):
        """Compute the EDA contribution to the forces, i.e. `q_e Z^* @ E(t)`"""
        Z = dstrip(self.eda.Born_Charges.bec)  # tensor of shape (nbeads,3xNatoms,3)
        E = dstrip(self.eda.Electric_Field.Efield)  # vector of shape (3)
        forces = Constants.e * Z @ E  # array of shape (nbeads,3xNatoms)
        return forces

    def step(self, step=None):
        if len(self.nmts) > 1:
            raise ValueError(
                "EDAIntegrator is not implemented with the Multiple Time Step algorithm (yet)."
            )
        super().step(step)


dproperties(EDAIntegrator, ["EDAforces", "mts_time", "time"])


class EDANVEIntegrator(EDAIntegrator, NVEIntegrator):
    """Integrator object for simulations with constant Number of particles, Volume, and Energy (NVE)
    using the Electric Dipole Approximation (EDA) when an external electric field is applied.
    """

    def pstep(self, level):
        # NVEIntegrator does not use 'super()' within 'pstep'
        # then we can not use 'super()' here.
        # We need to call the 'pstep' methods explicitly.
        NVEIntegrator.pstep(
            self, level
        )  # the driver is called here: add nuclear and electronic forces (DFT)
        EDAIntegrator.pstep(self, level)  # add the driving forces, i.e. q_e Z @ E(t)
        pass

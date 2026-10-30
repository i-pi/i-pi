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
from ipi.utils.array_backend import xp

from ipi.utils.depend import *
from ipi.utils.units import *
from ipi.utils.mathtools import matrix_exp
from ipi.utils.jit import compiled
from ipi.utils.messages import warning, verbosity
from ipi.engine.thermostats import Thermostat
from ipi.engine.cell import Cell

__all__ = ["Barostat", "BaroSCR", "BaroBZP", "BaroRGB", "BaroSCBZP", "BaroMTK"]


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


@compiled
def cell_propagators(v, dt):
    """Returns exp(v dt), exp(-v dt) and dt sinh(v dt) / (v dt) for a 3x3 cell
    velocity v, which propagate the positions and momenta under a constant
    cell velocity."""

    # exp([[x, 1], [0, 0]]) = [[exp(x), (exp(x) - 1) / x], [0, 1]], and the
    # average of (exp(x) - 1) / x for x = v dt and x = -v dt is sinh(v dt) / (v dt)
    vdt = v * dt
    eye, zeros = xp.eye(3), xp.zeros((3, 6))
    blocks = xp.stack(
        [xp.concat([xp.concat([x, eye], axis=1), zeros]) for x in (vdt, -vdt)]
    )
    # v dt is small (~1e-4), so a short series is accurate to rounding errors
    exps = matrix_exp(blocks, ntaylor=8, nsquare=4)
    sinh = 0.5 * dt * (exps[0, :3, 3:] + exps[1, :3, 3:])
    return exps[0, :3, :3], exps[1, :3, :3], sinh


class Barostat:
    """Base barostat class.

    Gives the standard methods and attributes needed in all the barostat classes.

    Used as barostat for mode "dummy".

    Attributes:
       beads: A beads object giving the atoms positions
       cell: A cell object giving the system box.
       forces: A forces object giving the virial and the forces acting on
          each bead.
       nm: An object to do the normal mode transformation.
       thermostat: A thermostat coupled to the barostat degrees of freedom.
       mdof: The number of atomic degrees of freedom
       vol_constraint: A boolean flag indicating whether the volume is constrained

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
        self,
        dt=None,
        temp=None,
        tau=None,
        ebaro=None,
        thermostat=None,
        vol_constraint=None,
        nmts=None,
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

        self.dt_ = depend_value(name="dt")
        if dt is not None:
            self.dt = dt
        else:
            self.dt = 1.0

        self.temp_ = depend_value(name="temp")
        if temp is not None:
            self.temp = temp
        else:
            self.temp = 1.0

        self.tau_ = depend_value(name="tau")
        if tau is not None:
            self.tau = tau
        else:
            self.tau = 1.0

        self.ebaro_ = depend_value(name="ebaro")
        if ebaro is not None:
            self.ebaro = ebaro
        else:
            self.ebaro = 0.0

        if thermostat is None:
            thermostat = Thermostat()
        self.thermostat = thermostat

        if vol_constraint is not None:
            self.vol_constraint = vol_constraint
        else:
            self.vol_constraint = False

        # temperature to the thermostat
        dpipe(self.temp_, self.thermostat.temp_)
        self.pext_ = depend_value(name="pext", value=-1.0)
        self.stressext_ = depend_array(name="stressext", value=-np.ones((3, 3), float))

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

        self.kstress_ = depend_value(
            name="kstress",
            func=self.get_kstress,
            dependencies=[beads.q_, beads.qc_, beads.pc_, forces.f_],
        )
        self.stress_ = depend_value(
            name="stress",
            func=self.get_stress,
            dependencies=[self.kstress_, cell.V_, forces.vir_],
        )

        self.pot_ = depend_value(name="pot", value=0.0)

        self.kin_ = depend_value(name="kin", value=0.0)

        self.cell_jacobian_ = depend_value(name="kin", value=0.0)

        # Stress depend objects for Suzuki-Chin PIMD
        self.kstress_sc_ = depend_value(
            name="kstress_sc",
            func=self.get_kstress_sc,
            dependencies=[
                beads.q_,
                beads.qc_,
                forces.fsc_part_2_,
                forces.f_,
            ],
        )

        self.stress_sc_ = depend_value(
            name="stress_sc",
            func=self.get_stress_sc,
            dependencies=[
                self.kstress_sc_,
                self.cell.V_,
                forces.vir_,
                forces.virssc_part_2_,
            ],
        )

        if fixdof is None:
            self.mdof = float(self.beads.natoms) * 3.0
        else:
            self.mdof = float(self.beads.natoms) * 3.0 - float(fixdof)

        # creates and connects timesteps for the different parts of the propagator
        self.nmtslevels = nmts
        self.qdt_ = depend_value(name="qdt", value=self.dt)
        self.pdt_ = depend_array(name="pdt", value=np.zeros(nmts, float), host=True)
        self.tdt_ = depend_value(name="tdt", value=self.dt)
        dpipe(self.tdt_, self.thermostat.dt_)

    # THESE SHOULD NOT BE USED ANYMORE
    def get_kstress(self):
        """Calculates the quantum centroid virial kinetic stress tensor
        estimator.
        """

        kst = xp.zeros((3, 3))
        q = self.beads.q
        qc = self.beads.qc
        pc = self.beads.pc
        m = self.beads.m
        na3 = 3 * self.beads.natoms
        fall = self.forces.f
        if self.bias is None:
            ball = fall * 0.00
        else:
            ball = self.bias.f

        for b in range(self.beads.nbeads):
            for i in range(3):
                for j in range(i, 3):
                    kst[i, j] -= (q[b, i:na3:3] - qc[i:na3:3]) @ (
                        fall[b, j:na3:3] + ball[b, j:na3:3]
                    )

        # NOTE: In order to have a well-defined conserved quantity, the Nf kT term in the
        # diagonal stress estimator must be taken from the centroid kinetic energy.
        for i in range(3):
            kst[i, i] += (pc[i:na3:3] @ (pc[i:na3:3] / m)) * self.beads.nbeads

        return kst

    def kstress_mts_sc(self, level):
        """Calculates the Suzuki-Chin quantum centroid virial kinetic stress tensor
        associated with the forces at a MTS level.
        """

        kst = xp.zeros((3, 3))
        q = self.beads.q
        qc = self.beads.qc
        pc = self.beads.pc
        m = self.beads.m
        na3 = 3 * self.beads.natoms
        fall = self.forces.mts_forces[level].f * (1 + self.forces.coeffsc_part_1)
        if self.bias is None or level != 0:
            ball = fall * 0.00
        else:
            ball = self.bias.f

        for b in range(self.beads.nbeads):
            for i in range(3):
                for j in range(i, 3):
                    kst[i, j] -= (q[b, i:na3:3] - qc[i:na3:3]) @ (
                        fall[b, j:na3:3] + ball[b, j:na3:3]
                    )

        if level == self.nmtslevels - 1:
            for i in range(3):
                kst[i, i] += (pc[i:na3:3] @ (pc[i:na3:3] / m)) * self.beads.nbeads

        return kst

    def kstress_mts(self, level):
        """Calculates the quantum centroid virial kinetic stress tensor
        associated with the forces at a MTS level.
        """

        fall = self.forces.mts_forces[level].f
        if level == 0 and self.bias is not None:
            # adds the bias. NB: don't += we don't want to overwrite the force
            fall = fall + self.bias.f

        beads = self.beads
        qqc = (beads.q - beads.qc).reshape(-1, 3)
        fall = fall.reshape(-1, 3)

        kst = -(qqc.T @ fall)
        # kst = np.zeros((3,3))
        # kst[0] = -np.dot(qqc[:,0], fall)
        # kst[1,1:] = -np.dot(qqc[:,1], fall[:,1:])
        # kst[2,2] = -np.dot(qqc[:,2], fall[:,2])

        if level == self.nmtslevels - 1:
            pc = beads.pc.reshape(-1, 3)
            m = beads.m
            nbeads = beads.nbeads

            # adds the centroid kinetic energy to the diagonal
            kst = kst + xp.eye(3) * ((nbeads / m) @ (pc * pc))

        return kst

    def get_kstress_sc(self):
        """Calculates the high order part of the Suzuki-Chin
        quantum centroid virial kinetic stress tensor
        associated with the forces at a MTS level.
        """

        kst = xp.zeros((3, 3))
        q = self.beads.q
        qc = self.beads.qc
        na3 = 3 * self.beads.natoms
        fall = self.forces.fsc_part_2

        for b in range(self.beads.nbeads):
            for i in range(3):
                for j in range(i, 3):
                    kst[i, j] -= (q[b, i:na3:3] - qc[i:na3:3]) @ fall[b, j:na3:3]
        return kst

    def get_stress(self):
        """Calculates the internal stress tensor."""

        bvir = xp.zeros((3, 3))
        if self.bias is not None:
            bvir[:] = self.bias.vir

        return (self.kstress + self.forces.virs + bvir) / self.cell.V

    def get_stress_sc(self):
        """Calculates the high order part of the Suzuki-Chin internal stress tensor."""
        return (self.kstress_sc + self.forces.virssc_part_2.sum(axis=0)) / self.cell.V

    def stress_mts_sc(self, level):
        """Calculates the internal Suzuki-Chin stress tensor
        associated with the forces at a MTS level.
        """

        bvir = xp.zeros((3, 3))
        if self.bias is not None and level == 0:
            bvir[:] = self.bias.vir

        return (
            self.kstress_mts_sc(level)
            + (
                self.forces.mts_forces[level].virs
                * (1 + self.forces.coeffsc_part_1).reshape((self.beads.nbeads, 1, 1))
            ).sum(axis=0)
            + bvir
        ) / self.cell.V

    def stress_mts(self, level):
        """Calculates the internal stress tensor
        associated with the forces at a MTS level.
        """

        stress = self.kstress_mts(level) + self.forces.mts_forces[level].vir
        if self.bias is not None and level == 0:
            stress = stress + self.bias.vir

        return stress / self.cell.V

    def pstep(self, level=0):
        """Dummy momenta propagator step."""

        pass

    def qcstep(self, level=0):
        """Dummy centroid position propagator step."""

        pass


dproperties(
    Barostat,
    [
        "dt",
        "temp",
        "tau",
        "ebaro",
        "pext",
        "stressext",
        "kstress",
        "stress",
        "pot",
        "kin",
        "cell_jacobian",
        "kstress_sc",
        "stress_sc",
        "qdt",
        "pdt",
        "tdt",
    ],
)


class BaroSCR(Barostat):
    """Bernetti--Bussi stochastic cell-rescaling barostat.

    Mode ``stochastic-rescaling`` implements the reversible, isotropic
    Trotter scheme with instantaneous kinetic pressure described in
    Bernetti and Bussi, J. Chem. Phys. 153, 114107 (2020)
    [doi:10.1063/5.0020514]. Positions and cell vectors are scaled by the
    same factor, while momenta are scaled by its inverse.

    This implementation is restricted to classical dynamics. The particle
    thermostat remains responsible for temperature control; stochastic cell
    rescaling does not have a separate piston or cell thermostat.
    """

    def __init__(
        self,
        dt=None,
        temp=None,
        tau=None,
        ebaro=None,
        thermostat=None,
        pext=None,
        compressibility=None,
    ):
        """Initializes the stochastic cell-rescaling barostat.

        Args:
            dt: Optional simulation time step.
            temp: Optional simulation temperature.
            tau: Pressure relaxation time.
            ebaro: Accumulated effective-energy correction, used on restart.
            thermostat: Must be the dummy thermostat. Kept for input
                compatibility with other barostats.
            pext: Optional external pressure.
            compressibility: Estimate of the isothermal compressibility.
        """

        super(BaroSCR, self).__init__(dt, temp, tau, ebaro, thermostat)

        self.compressibility_ = depend_value(
            name="compressibility",
            value=-1.0 if compressibility is None else compressibility,
        )
        self.prng = None
        self._scr_state = None

        if pext is not None:
            self.pext = pext

    def bind(self, beads, nm, cell, forces, bias=None, prng=None, fixdof=None, nmts=1):
        """Binds the physical system and random-number generator."""

        super(BaroSCR, self).bind(beads, nm, cell, forces, bias, prng, fixdof, nmts)

        if self.tau <= 0.0:
            raise ValueError(
                "The stochastic cell-rescaling relaxation time must be positive."
            )
        if self.compressibility <= 0.0:
            raise ValueError(
                "The stochastic cell-rescaling compressibility must be positive."
            )
        if self.beads.nbeads != 1:
            raise ValueError(
                "Stochastic cell rescaling currently supports classical "
                "dynamics only (nbeads=1)."
            )
        if type(self.thermostat) is not Thermostat:
            raise ValueError(
                "Stochastic cell rescaling does not use a separate cell thermostat."
            )
        if prng is None:
            raise ValueError(
                "Stochastic cell rescaling requires the simulation random-number "
                "generator."
            )

        self.prng = prng
        self.pot_ = depend_value(
            name="pot", func=self.get_pot, dependencies=[cell.V_, self.pext_]
        )
        self.cell_jacobian_ = depend_value(
            name="cell_jacobian",
            func=self.get_cell_jacobian,
            dependencies=[cell.V_, self.temp_],
        )

    def get_pot(self):
        """Returns the external-pressure contribution to the enthalpy."""

        return self.pext * self.cell.V

    def get_cell_jacobian(self):
        """Returns the Jacobian contribution for lambda = sqrt(volume)."""

        return -0.5 * Constants.kb * self.temp * xp.log(self.cell.V)

    def get_lambda_force(self):
        """Returns the drift force acting on lambda = sqrt(volume)."""

        volume = float(self.cell.V)
        if not np.isfinite(volume) or volume <= 0.0:
            raise ValueError(
                "Stochastic cell rescaling encountered a non-positive or "
                "non-finite volume."
            )
        pressure = (
            sum(
                float(xp.linalg.trace(self.stress_mts(level)))
                for level in range(self.nmtslevels)
            )
            / 3.0
        )
        if not np.isfinite(pressure):
            raise ValueError(
                "Stochastic cell rescaling encountered a non-finite internal "
                "pressure."
            )

        lam = np.sqrt(volume)
        return (
            -2.0
            * lam
            * (self.pext - pressure - Constants.kb * self.temp / (2.0 * volume))
        )

    def prepare(self):
        """Draws a cell move and returns its isotropic length scale."""

        lambda_force = self.get_lambda_force()
        volume = float(self.cell.V)
        lam = np.sqrt(volume)
        coupling_dt = 2.0 * self.qdt
        kbt = Constants.kb * self.temp
        diffusion = kbt * self.compressibility / (4.0 * self.tau)
        delta_lambda = (
            diffusion * lambda_force * coupling_dt / kbt
            + np.sqrt(2.0 * diffusion * coupling_dt) * self.prng.g
        )
        new_lambda = lam + delta_lambda

        if not np.isfinite(new_lambda) or new_lambda <= 0.0:
            raise ValueError(
                "Stochastic cell rescaling proposed a non-positive or non-finite "
                "sqrt(volume). Increase tau, reduce stride or the MD time step, "
                "and check the compressibility."
            )

        delta_log_lambda = np.log1p(delta_lambda / lam)
        scale = np.exp((2.0 / 3.0) * delta_log_lambda)
        if scale < 0.99 or scale > 1.01:
            warning(
                "Stochastic cell rescaling changed a cell length by more than "
                "1%. Consider increasing tau or reducing stride.",
                verbosity.low,
            )

        self._scr_state = {
            "volume": volume,
            "lambda": lam,
            "lambda_force": lambda_force,
            "delta_lambda": delta_lambda,
            "delta_log_lambda": delta_log_lambda,
            "diffusion": diffusion,
            "coupling_dt": coupling_dt,
        }
        return scale

    def qcstep(self):
        """Propagates one half of the reversible cell-rescaling map.

        The standard NPT integrator calls this method twice around the free
        ring-polymer propagation. The first call draws the stochastic cell
        move and applies drift--rescale; the second applies rescale--drift and
        closes the effective-energy bookkeeping.
        """

        masses = self.nm.dynm3[0]

        if self._scr_state is None:
            scale = self.prepare()
            half_scale = np.sqrt(scale)
            self._scr_state["half_scale"] = half_scale

            self.nm.qnm_[0, :] += self.nm.pnm[0] * self.qdt / masses
            self.nm.qnm_[0, :] *= half_scale
            self.nm.pnm_[0, :] *= 1.0 / half_scale
            self.cell.h_ *= half_scale
        else:
            half_scale = self._scr_state["half_scale"]

            self.nm.qnm_[0, :] *= half_scale
            self.nm.pnm_[0, :] *= 1.0 / half_scale
            self.cell.h_ *= half_scale
            self.nm.qnm_[0, :] += self.nm.pnm[0] * self.qdt / masses
            self.finalize()

    def finalize(self):
        """Accumulates the reversible effective-energy correction."""

        if self._scr_state is None:
            raise RuntimeError(
                "Stochastic cell-rescaling finalize called without a pending move."
            )

        state = self._scr_state
        new_volume = float(self.cell.V)
        new_lambda_force = self.get_lambda_force()
        kbt = Constants.kb * self.temp

        correction = self.pext * (new_volume - state["volume"])
        correction -= kbt * state["delta_log_lambda"]
        correction += (
            0.5 * state["delta_lambda"] * (state["lambda_force"] + new_lambda_force)
        )
        correction += (
            state["diffusion"]
            * state["coupling_dt"]
            / (4.0 * kbt)
            * (new_lambda_force**2 - state["lambda_force"] ** 2)
        )
        self.ebaro += correction
        self._scr_state = None


dproperties(BaroSCR, ["compressibility"])


class BaroBZP(Barostat):
    """Bussi-Zykova-Parrinello barostat class.

    Mode "isotropic" implements the Bussi-Zykova-Parrinello barostat.

    Isotropically scales the unit cell volume while sampling the
    isothermal isobaric ensemble.

    Bussi, Zykova-Timan, Parinello, J. Chem. Phys. 130, 074101 (2009)
        [doi:10.1063/1.3073889]

    Just extends the standard class adding finite-dt propagators for the
    barostat velocities, positions, piston.

    Generates dynamics with a stochastic barostat. Implementation details:
    Ceriotti, More, Manolopoulos, Comp. Phys. Comm. 185, 1019 (2013)

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

        self.p_ = depend_array(name="p", value=xp.zeros(1))

        if p is not None:
            self.p = xp.asarray([p])
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

        # obtain the thermostat mass from the given time constant
        # note that the barostat temperature is nbeads times the physical T
        self.m_ = depend_array(
            name="m",
            value=xp.zeros(1),
            func=(
                lambda: xp.asarray(
                    [self.tau**2 * 3 * self.beads.natoms * Constants.kb * self.temp]
                )
            ),
            dependencies=[self.tau_, self.temp_],
        )

        # binds the thermostat to the piston degrees of freedom
        self.thermostat.bind(pm=[self.p_, self.m_], prng=prng)

        # barostat elastic energy
        self.pot_ = depend_value(
            name="pot", func=self.get_pot, dependencies=[cell.V_, self.pext_]
        )

        self.kin_ = depend_value(
            name="kin",
            func=(lambda: 0.5 * self.p[0] ** 2 / self.m[0]),
            dependencies=[self.p_, self.m_],
        )

        # defines the term that accounts for the explicit dependence of the volume on the ensemble
        self.cell_jacobian_ = depend_value(
            name="cell_jacobian",
            func=self.get_cell_jacobian,
            dependencies=[self.cell.V_, self.temp_],
        )

        # the barostat energy must be computed from bits & pieces (overwrite the default)
        self.ebaro_ = depend_value(
            name="ebaro",
            func=self.get_ebaro,
            dependencies=[
                self.kin_,
                self.pot_,
                self.cell.V_,
                self.temp_,
                self.thermostat.ethermo_,
            ],
        )

    def get_pot(self):
        """Calculates the elastic strain energy of the cell."""

        # NOTE: since there are nbeads replicas of the unit cell, the enthalpy contains a nbeads factor
        return self.cell.V * self.pext * self.beads.nbeads

    def get_cell_jacobian(self):
        """Calculates the energy term that accounts for the size of the box"""

        return -1.0 * xp.log(self.cell.V) * Constants.kb * self.temp

    def get_ebaro(self):
        """Calculates the barostat conserved quantity."""

        return self.thermostat.ethermo + self.kin + self.pot + self.cell_jacobian

    def pstep(self, level=0):
        """Propagates the momentum of the barostat."""

        # we are assuming then that p the coupling between p^2 and dp/dt only involves the fast force
        dt = self.pdt[
            level
        ]  # this is already set to be half a time step at the specified MTS depth

        # computes the pressure associated with the forces at each MTS level.
        press = xp.linalg.trace(self.stress_mts(level)) / 3.0
        self.p_ += dt * 3.0 * (self.cell.V * press)

        # integerates the kinetic part of the pressure with the force at the inner-most level.
        if level == self.nmtslevels - 1:
            nbeads = self.beads.nbeads
            press = 0
            self.p_ += (
                3.0
                * dt
                * (
                    self.cell.V * (press - nbeads * self.pext)
                    + Constants.kb * self.temp
                )
            )

            pc = self.beads.pc
            fc = self.forces.mts_forces[level].f.sum(axis=0) / nbeads
            fc_m = fc / self.beads.m3[0]

            dt2 = dt * dt * self.beads.nbeads
            dt3 = dt / 3.0
            self.p_ += ((pc @ fc_m) + dt3 * (fc @ fc_m)) * dt2

    def qcstep(self):
        """Propagates the centroid position and momentum and the volume."""

        v = self.p[0] / self.m[0]
        halfdt = (
            self.qdt
        )  # this is set to half the inner loop in all integrators that use a barostat
        expq, expp = (xp.exp(v * halfdt), xp.exp(-v * halfdt))

        m = self.beads.m3[0]

        self.nm.qnm_[0, :] = (self.nm.qnm[0, :] * expq) + (
            (expq - expp) / (2.0 * v)
        ) * (self.nm.pnm[0, :] / m)
        self.nm.pnm_[0, :] = expp * self.nm.pnm[0, :]

        self.cell.h_ *= expq


dproperties(BaroBZP, ["p", "m"])


class BaroSCBZP(Barostat):
    """The Suzuki Chin Bussi-Zykova-Parrinello barostat class.

    Mode "sc-isotropic" implements the Bussi-Zykova-Parrinello
    barostat for Suzuki-Chin path integral molecular dynamics.

    Isotropically scales the unit cell volume while sampling the
    isothermal isobaric ensemble.

    This barostat is suitable for simulating liquids.

    Kapil et al., J. Chem. Theory Comput. 15, 5 (2019)
        [10.1021/acs.jctc.8b01297]

    Just extends the standard class adding finite-dt propagators for the
    barostat velocities, positions, piston.

    Generates dynamics with a stochastic barostat. Implementation details:
    Ceriotti, More, Manolopoulos, Comp. Phys. Comm. 185, 1019 (2013)

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

        self.p_ = depend_array(name="p", value=xp.zeros(1))

        if p is not None:
            self.p = xp.asarray([p])
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

        # obtain the thermostat mass from the given time constant
        self.m_ = depend_array(
            name="m",
            value=xp.zeros(1),
            func=(
                lambda: xp.asarray(
                    [self.tau**2 * 3 * self.beads.natoms * Constants.kb * self.temp]
                )
            ),
            dependencies=[self.tau_, self.temp_],
        )

        # binds the thermostat to the piston degrees of freedom
        self.thermostat.bind(pm=[self.p_, self.m_], prng=prng)

        # barostat elastic energy
        self.pot_ = depend_value(
            name="pot", func=self.get_pot, dependencies=[cell.V_, self.pext_]
        )

        self.kin_ = depend_value(
            name="kin",
            func=(lambda: 0.5 * self.p[0] ** 2 / self.m[0]),
            dependencies=[self.p_, self.m_],
        )

        # defines the term that accounts for the explicit dependence of the ensemble probability on the volume
        self.cell_jacobian_ = depend_value(
            name="cell_jacobian",
            func=self.get_cell_jacobian,
            dependencies=[self.cell.V_, self.temp_],
        )

        # the barostat energy must be computed from bits & pieces (overwrite the default)
        self.ebaro_ = depend_value(
            name="ebaro",
            func=self.get_ebaro,
            dependencies=[
                self.kin_,
                self.pot_,
                self.cell_jacobian_,
                self.thermostat.ethermo_,
            ],
        )

    def get_pot(self):
        """Calculates the elastic strain energy of the cell."""

        # NOTE: since there are nbeads replicas of the unit cell, the enthalpy contains a nbeads factor
        return self.cell.V * self.pext * self.beads.nbeads

    def get_cell_jacobian(self):
        """Calculates the energy term that accounts for the size of the box"""

        return -1.0 * xp.log(self.cell.V) * Constants.kb * self.temp

    def get_ebaro(self):
        """Calculates the barostat conserved quantity."""

        return self.thermostat.ethermo + self.kin + self.pot + self.cell_jacobian

    def pscstep(self):
        """Propagates the momentum of the barostat with respect to the
        high order part of the Suzuki-Chin stress"""

        # integrates with respect to the "high order" part of the stress with a timestep of dt /2
        press = xp.linalg.trace(self.stress_sc) / 3.0
        self.p_ += self.dt / 2 * 3.0 * (self.cell.V * press)

    def pstep(self, level=0):
        """Propagates the momentum of the barostat."""

        # we are assuming then that p the coupling between p^2 and dp/dt only involves the fast force
        dt = self.pdt[
            level
        ]  # this is already set to be half a time step at the specified MTS depth
        dt2 = dt**2
        dt3 = dt**3 / 3.0

        # computes the pressure associated with the forces at each MTS level and adds the +- 1/3 SC correction.
        press = xp.linalg.trace(self.stress_mts_sc(level)) / 3.0
        self.p_ += dt * 3.0 * (self.cell.V * press)

        # integerates the kinetic part of the pressure with the force at the inner-most level.
        if level == self.nmtslevels - 1:
            press = 0
            self.p_ += (
                dt
                * 3.0
                * (
                    self.cell.V * (press - self.beads.nbeads * self.pext)
                    + Constants.kb * self.temp
                )
            )
            pc = self.beads.pc
            fc = (
                self.forces.mts_forces[level].f * (1 + self.forces.coeffsc_part_1)
            ).sum(axis=0) / self.beads.nbeads
            m = self.beads.m3[0]

            self.p_ += (
                dt2 * (pc @ (fc / m)) + dt3 * (fc @ (fc / m))
            ) * self.beads.nbeads

    def qcstep(self):
        """Propagates the centroid position and momentum and the volume."""

        v = self.p[0] / self.m[0]
        halfdt = (
            self.qdt
        )  # this is set to half the inner loop in all integrators that use a barostat
        expq, expp = (xp.exp(v * halfdt), xp.exp(-v * halfdt))

        m = self.beads.m3[0]

        self.nm.qnm_[0, :] *= expq
        self.nm.qnm_[0, :] += ((expq - expp) / (2.0 * v)) * (self.nm.pnm[0, :] / m)
        self.nm.pnm_[0, :] *= expp

        self.cell.h_ *= expq


dproperties(BaroSCBZP, ["p", "m"])


class BaroRGB(Barostat):
    """Raiteri-Gale-Bussi constant stress barostat class.

    Mode "anisotropic" implements the Raiteri-Gale-Bussi barostat which enables
    cell fluctuations (lengths and angles) at constant external stress.

    It is suitable for simulating solids at given external (non-diagonal) stresses
    and requires specifying a reference cell for estimating strain.

     For diagonal stresses (or external pressures) the "flexible" and the
     "anisotropic" modes should give very similar results.

    Note that this ensemble is valid only within the elastic limit of small strains.

    Raiteri, Gale, Bussi, J. Phys. Condens. Matter, 23, 334213 (2011)
        [doi:10.1088/0953-8984/23/33/334213]

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
        # 6-vector or as a 3x3 upper triangular tensor, which are synchronized.
        # The rows of `select6` pick these elements from the flattened tensor
        self.select6 = xp.take(xp.eye(9), xp.asarray([0, 4, 8, 1, 2, 5]), axis=0)
        self.p_ = depend_array(name="p", value=xp.zeros((3, 3)))
        self.p6_ = depend_array(name="p6", value=xp.zeros(6))
        dsynchro(self.p_, self.p6_, self.get_6to3x3, self.get_3x3to6)

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
        self.hmask_ = depend_array(name="hmask", value=hmask)
        # number of ones in the UT part of the mask
        self.L = xp.asarray(
            np.diag([hmask[0].sum(), hmask[1, 1:].sum(), hmask[2, 2:].sum()])
        )

        if stressext is not None:
            self.stressext = stressext
        else:
            self.stressext_[:] = -1.0

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

        self.m_ = depend_array(
            name="m",
            value=xp.zeros(1),
            func=(
                lambda: xp.asarray(
                    [self.tau**2 * self.beads.natoms * Constants.kb * self.temp]
                )
            ),
            dependencies=[self.tau_, self.temp_],
        )

        self.m6_ = depend_array(
            name="m6",
            value=xp.zeros(6),
            func=(lambda: xp.ones(6) * self.m[0]),
            dependencies=[self.m_],
        )

        # overrides definition of pot to depend on the many things it depends on for anisotropic cell
        self.pot_ = depend_value(
            name="pot",
            func=self.get_pot,
            dependencies=[
                self.cell.h_,
                self.h0.h_,
                self.h0.V_,
                self.h0.ih_,
                self.stressext_,
            ],
        )

        # binds the thermostat to the piston degrees of freedom
        self.thermostat.bind(pm=[self.p6_, self.m6_], prng=prng)

        self.kin_ = depend_value(
            name="kin",
            func=(lambda: 0.5 * xp.linalg.trace(self.p.T @ self.p) / self.m[0]),
            dependencies=[self.p_, self.m_],
        )

        # defines the term that accounts for the explicit dependence of the ensemble on the cell
        self.cell_jacobian_ = depend_value(
            name="cell_jacobian",
            func=self.get_cell_jacobian,
            dependencies=[self.cell.h_, self.temp_],
        )

        # the barostat energy must be computed from bits & pieces (overwrite the default)
        self.ebaro_ = depend_value(
            name="ebaro",
            func=self.get_ebaro,
            dependencies=[
                self.kin_,
                self.pot_,
                self.cell_jacobian_,
                self.thermostat.ethermo_,
            ],
        )

    def get_3x3to6(self):
        return self.select6 @ self.p.reshape(9)

    def get_6to3x3(self):
        return (self.p6 @ self.select6).reshape(3, 3)

    def get_pot(self):
        """Calculates the elastic strain energy of the cell."""

        # NOTE: since there are nbeads replicas of the unit cell, the enthalpy contains a nbeads factor
        eps = self.cell.h @ self.h0.ih
        eps = eps.T @ eps
        eps = 0.5 * (eps - xp.eye(3))

        return self.h0.V * xp.linalg.trace(self.stressext @ eps) * self.beads.nbeads

    def get_cell_jacobian(self):
        """Calculates the energy term that accounts for the size of the box."""

        cj = sum((3 - i) * xp.log(self.cell.h[i, i]) for i in range(3))
        return -1.0 * Constants.kb * self.temp * cj

    def get_ebaro(self):
        """Calculates the barostat conserved quantity."""

        return self.thermostat.ethermo + self.kin + self.pot + self.cell_jacobian

    def pstep(self, level=0):
        """Propagates the momenta for half a time step."""

        dt = self.pdt[level]
        dt2 = dt**2
        dt3 = dt**3 / 3.0

        hh0 = self.cell.h @ self.h0.ih
        pi_ext = hh0 @ self.stressext @ hh0.T * self.h0.V / self.cell.V

        stress = self.stress_mts(level)
        self.p_ += dt * (self.cell.V * xp.triu(stress))

        # integerates the kinetic part of the stress with the force at the inner-most level.
        if level == self.nmtslevels - 1:
            self.p_ += dt * (
                self.cell.V * xp.triu(-self.beads.nbeads * pi_ext)
                + Constants.kb * self.temp * self.L
            )

            pc = self.beads.pc.reshape(self.beads.natoms, 3)
            fc = (
                self.forces.mts_forces[level]
                .f.sum(axis=0)
                .reshape(self.beads.natoms, 3)
                / self.beads.nbeads
            )
            fcTonm = (fc / self.beads.m3[0].reshape(self.beads.natoms, 3)).T

            self.p_ += (
                xp.triu(dt2 * (fcTonm @ pc) + dt3 * (fcTonm @ fc)) * self.beads.nbeads
            )

        # now apply the mask (and accumulate the associated change in conserved quantity)
        # we use the thermostat conserved quantity accumulator, so we don't need to create a new one
        self.thermostat.ethermo += self.kin
        self.p_ *= self.hmask
        self.thermostat.ethermo -= self.kin

    def qcstep(self):
        """Propagates the centroid position and momentum and the volume."""

        expq, expp, sinh = cell_propagators(self.p / self.m[0], self.qdt)

        q = self.nm.qnm[0].reshape(self.beads.natoms, 3)
        p = self.nm.pnm[0]

        q = q @ expq.T
        q += (p / self.beads.m3[0]).reshape(self.beads.natoms, 3) @ sinh.T
        p = p.reshape(self.beads.natoms, 3) @ expp.T

        self.nm.qnm_[0] = q.reshape(self.beads.natoms * 3)
        self.nm.pnm_[0] = p.reshape(self.beads.natoms * 3)

        self.cell.h = expq @ self.cell.h


dproperties(BaroRGB, ["p", "m", "p6", "hmask", "m6"])


class BaroMTK(Barostat):
    """Martyna-Tobias-Klein flexible cell, constant pressure barostat class.

    Mode "flexible" implements the path integral version of Martyna-Tuckerman-Tobias-Klein
    barostat, which incorporates full cell fluctuations (lengths and angles) while
    sampling the isothermal isobaric ensemble.

    This is suitable for anisotropic systems such as molecular solids.

    For diagonal stresses (or external pressures) the "flexible" and the
    "anisotropic" modes should give very similar results.

    Martyna and Hughes, J. Chem. Phys 110, 3275 (1999)
        [doi:10.1063/1.478193]

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
        vol_constraint=None,
        pext=None,
        p=None,
        hfix=None,
    ):
        """Initializes MTK barostat.

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
        vol_constraint: Optional boolean specifying whether or not the volume is constrained.
        p: Optional initial volume conjugate momentum. Defaults to 0.
        """

        super(BaroMTK, self).__init__(dt, temp, tau, ebaro, thermostat, vol_constraint)

        # non-zero elements of the cell momentum are only
        # pxx pyy pzz pxy pxz pyz, but we want to access it either as a
        # 6-vector or as a 3x3 upper triangular tensor, which are synchronized.
        # The rows of `select6` pick these elements from the flattened tensor
        self.select6 = xp.take(xp.eye(9), xp.asarray([0, 4, 8, 1, 2, 5]), axis=0)
        self.p_ = depend_array(name="p", value=xp.zeros((3, 3)))
        self.p6_ = depend_array(name="p6", value=xp.zeros(6))
        dsynchro(self.p_, self.p6_, self.get_6to3x3, self.get_3x3to6)

        if p is not None:
            self.p = p
        else:
            self.p = 0.0

        if hfix is None:
            hfix = []
        self.hfix = hfix
        hmask = mask_from_fix(hfix)

        # mask to zero out components of the cell velocity, to implement cell-boundary constraints
        self.hmask_ = depend_array(name="hmask", value=hmask)
        # number of ones in the UT part of the mask
        self.L = xp.asarray(
            np.diag([hmask[0].sum(), hmask[1, 1:].sum(), hmask[2, 2:].sum()])
        )

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

        self.m_ = depend_array(
            name="m",
            value=xp.zeros(1),
            func=(
                lambda: xp.asarray(
                    [self.tau**2 * self.beads.natoms * Constants.kb * self.temp]
                )
            ),
            dependencies=[self.tau_, self.temp_],
        )

        self.m6_ = depend_array(
            name="m6",
            value=xp.zeros(6),
            func=(lambda: xp.ones(6) * self.m[0]),
            dependencies=[self.m_],
        )

        # overrides definition of pot to depend on the many things it depends on for anisotropic cell
        self.pot_ = depend_value(
            name="pot", func=self.get_pot, dependencies=[self.cell.h_, self.pext_]
        )

        # binds the thermostat to the piston degrees of freedom
        self.thermostat.bind(pm=[self.p6_, self.m6_], prng=prng)

        self.kin_ = depend_value(
            name="kin",
            func=(lambda: 0.5 * xp.linalg.trace(self.p.T @ self.p) / self.m[0]),
            dependencies=[self.p_, self.m_],
        )

        # defines the term that accounts for the explicit dependence of the ensemble on the cell
        self.cell_jacobian_ = depend_value(
            name="cell_jacobian",
            func=self.get_cell_jacobian,
            dependencies=[self.cell.h_, self.temp_],
        )

        # the barostat energy must be computed from bits & pieces (overwrite the default)
        self.ebaro_ = depend_value(
            name="ebaro",
            func=self.get_ebaro,
            dependencies=[
                self.kin_,
                self.pot_,
                self.cell_jacobian_,
                self.thermostat.ethermo_,
            ],
        )

    def get_3x3to6(self):
        return self.select6 @ self.p.reshape(9)

    def get_6to3x3(self):
        return (self.p6 @ self.select6).reshape(3, 3)

    def get_pot(self):
        """Calculates the elastic strain energy of the cell."""

        # return self.h0.V * np.trace(np.dot(self.stressext, eps)) * self.beads.nbeads
        return self.cell.V * self.pext * self.beads.nbeads
        # p = np.trace(self.stressext)/3
        # return (p*(self.cell.V-self.h0.V) +
        #    self.h0.V * np.trace(np.dot(self.stressext-np.eye(3)*p, eps)) ) * self.beads.nbeads

    def get_cell_jacobian(self):
        """Calculates the energy term that accounts for the size of the box."""

        cj = sum((3 - i) * xp.log(self.cell.h[i, i]) for i in range(3))
        return -1.0 * Constants.kb * self.temp * cj

    def get_ebaro(self):
        """Calculates the barostat conserved quantity."""

        return self.thermostat.ethermo + self.kin + self.pot + self.cell_jacobian

    def pstep(self, level=0):
        """Propagates the momenta for half a time step."""

        dt = self.pdt[level]
        dt2 = dt**2
        dt3 = dt**3 / 3.0

        pi_ext = (
            xp.eye(3) * self.pext
        )  # np.dot(hh0, np.dot(self.stressext, hh0.T)) * self.h0.V / self.cell.V
        #        L = np.diag([3, 2, 1])

        stress = self.stress_mts(level)
        self.p_ += dt * (self.cell.V * xp.triu(stress))

        # integerates the kinetic part of the stress with the force at the inner-most level.
        if level == self.nmtslevels - 1:
            self.p_ += dt * (
                self.cell.V * xp.triu(-self.beads.nbeads * pi_ext)
                + Constants.kb * self.temp * self.L
            )

            pc = self.beads.pc.reshape(self.beads.natoms, 3)
            fc = (
                self.forces.mts_forces[level]
                .f.sum(axis=0)
                .reshape(self.beads.natoms, 3)
                / self.beads.nbeads
            )
            fcTonm = (fc / self.beads.m3[0].reshape(self.beads.natoms, 3)).T

            self.p_ += (
                xp.triu(dt2 * (fcTonm @ pc) + dt3 * (fcTonm @ fc)) * self.beads.nbeads
            )

        # now apply the mask (and accumulate the associated change in conserved quantity)
        # we use the thermostat conserved quantity accumulator, so we don't need to create a new one
        self.thermostat.ethermo += self.kin
        self.p_ *= self.hmask
        if self.vol_constraint is True:
            # a traceless momentum tensor will not change the volume, see
            # Rogge, S. M. J. et al. Theory Comput. 11, 5583–5597 (2015) DOI: 10.1021/acs.jctc.5b00748
            self.p_ -= xp.eye(3) * xp.linalg.trace(self.p) / 3.0
        self.thermostat.ethermo -= self.kin

    def qcstep(self):
        """Propagates the centroid position and momentum and the volume."""

        expq, expp, sinh = cell_propagators(self.p / self.m[0], self.qdt)

        q = self.nm.qnm[0].reshape(self.beads.natoms, 3)
        p = self.nm.pnm[0]

        q = q @ expq.T
        q += (p / self.beads.m3[0]).reshape(self.beads.natoms, 3) @ sinh.T
        p = p.reshape(self.beads.natoms, 3) @ expp.T

        self.nm.qnm_[0] = q.reshape(self.beads.natoms * 3)
        self.nm.pnm_[0] = p.reshape(self.beads.natoms * 3)

        self.cell.h = expq @ self.cell.h


dproperties(BaroMTK, ["p", "m", "p6", "hmask", "m6"])

"""Classes that deal with constant pressure and constant stress simulations.

Contains the algorithms which propagate the positions, the momenta and the cell
in the constant pressure and constant stress ensembles, and the properties
directly related to them, such as the internal and external pressure and stress.
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

    # the lower-triangular elements are mapped to the upper-triangular ones,
    # the only ones that change
    m2f = dict(
        xx=(0, 0),
        xy=(0, 1),
        xz=(0, 2),
        yx=(0, 1),
        yy=(1, 1),
        yz=(1, 2),
        zx=(0, 2),
        zy=(1, 2),
        zz=(2, 2),
        offdiagonal=[(0, 1), (0, 2), (1, 2)],
    )

    hmask = np.ones(shape=(3, 3))
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
       bias: A forces object giving the bias.
       nm: An object to do the normal mode transformation.
       thermostat: A thermostat coupled to the barostat degrees of freedom.
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
       kin: The kinetic energy of the barostat.
       cell_jacobian: The energy term that accounts for the measure of the cell.
    """

    def __init__(
        self,
        dt=None,
        temp=None,
        tau=None,
        ebaro=None,
        thermostat=None,
        vol_constraint=None,
    ):
        """Initialises base barostat class.

        Args:
           dt: Optional float giving the time step for the algorithms. Defaults
              to the simulation dt.
           temp: Optional float giving the temperature for the thermostat.
              Defaults to the simulation temp.
           tau: Optional float giving the time scale associated with the barostat.
           ebaro: Optional float giving the conserved quantity already stored
              in the barostat initially. Used on restart.
           thermostat: The thermostat connected to the barostat degree of freedom.
           vol_constraint: Optional boolean, True if the volume is constrained.
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

    def bind(self, beads, nm, cell, forces, bias=None, prng=None, nmts=1):
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
           bias: The forces object of the bias.
           prng: The parent PRNG to bind the thermostat to
           nmts: The number of levels of multiple time stepping.
        """

        self.beads = beads
        self.cell = cell
        self.forces = forces
        self.bias = bias
        self.nm = nm

        self.pot_ = depend_value(name="pot", value=0.0)

        self.kin_ = depend_value(name="kin", value=0.0)

        self.cell_jacobian_ = depend_value(name="cell_jacobian", value=0.0)

        # creates and connects timesteps for the different parts of the propagator
        self.nmtslevels = nmts
        self.qdt_ = depend_value(name="qdt", value=self.dt)
        self.pdt_ = depend_array(name="pdt", value=np.zeros(nmts, float), host=True)
        self.tdt_ = depend_value(name="tdt", value=self.dt)
        dpipe(self.tdt_, self.thermostat.dt_)

    def mts_force(self, level):
        """Returns the force at the MTS level `level`."""

        return self.forces.mts_forces[level].f

    def kick_force(self, level):
        """Returns the force that kicks the beads at the MTS level `level`: that
        of the level and, at level 0, the bias."""

        force = self.mts_force(level)
        if level == 0 and self.bias is not None:
            force = force + self.bias.f
        return force

    def mts_virial(self, level):
        """Returns the virial at the MTS level `level`."""

        return self.forces.mts_forces[level].vir

    def kstress_mts(self, level):
        """Calculates the quantum centroid virial kinetic stress tensor
        associated with the forces at a MTS level.
        """

        beads = self.beads
        qqc = (beads.q - beads.qc).reshape(-1, 3)
        fall = self.kick_force(level).reshape(-1, 3)

        kst = -(qqc.T @ fall)

        if level == self.nmtslevels - 1:
            pc = beads.pc.reshape(-1, 3)
            m = beads.m
            nbeads = beads.nbeads

            # adds the centroid kinetic energy to the diagonal
            kst = kst + xp.eye(3) * ((nbeads / m) @ (pc * pc))

        return kst

    def stress_mts(self, level):
        """Calculates the internal stress tensor
        associated with the forces at a MTS level.
        """

        stress = self.kstress_mts(level) + self.mts_virial(level)
        if self.bias is not None and level == 0:
            stress = stress + self.bias.vir

        return stress / self.cell.V

    def centroid_arrays(self):
        """Returns the depend arrays whose first rows hold the centroid
        positions and momenta: those of the beads if there is only one, which
        avoids the transformations to normal modes and back, and those of the
        normal modes otherwise."""

        if self.beads.nbeads == 1:
            return self.beads.q_, self.beads.p_
        return self.nm.qnm_, self.nm.pnm_

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
        "pot",
        "kin",
        "cell_jacobian",
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
        bulk_modulus=None,
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
            bulk_modulus: Estimate of the isothermal bulk modulus.
        """

        super(BaroSCR, self).__init__(dt, temp, tau, ebaro, thermostat)

        self.bulk_modulus_ = depend_value(
            name="bulk_modulus",
            value=-1.0 if bulk_modulus is None else bulk_modulus,
        )
        self.prng = None
        self._scr_state = None

        if pext is not None:
            self.pext = pext

    def bind(self, beads, nm, cell, forces, bias=None, prng=None, nmts=1):
        """Binds the physical system and random-number generator."""

        super(BaroSCR, self).bind(beads, nm, cell, forces, bias, prng, nmts)

        if self.tau <= 0.0:
            raise ValueError(
                "The stochastic cell-rescaling relaxation time must be positive."
            )
        if self.bulk_modulus <= 0.0:
            raise ValueError(
                "The stochastic cell-rescaling bulk modulus must be positive."
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
        diffusion = kbt / (4.0 * self.tau * self.bulk_modulus)
        delta_lambda = (
            diffusion * lambda_force * coupling_dt / kbt
            + np.sqrt(2.0 * diffusion * coupling_dt) * self.prng.g
        )
        new_lambda = lam + delta_lambda

        if not np.isfinite(new_lambda) or new_lambda <= 0.0:
            raise ValueError(
                "Stochastic cell rescaling proposed a non-positive or non-finite "
                "sqrt(volume). Increase tau, reduce stride or the MD time step, "
                "and check the bulk modulus."
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
        q_, p_ = self.centroid_arrays()

        if self._scr_state is None:
            scale = self.prepare()
            half_scale = np.sqrt(scale)
            self._scr_state["half_scale"] = half_scale

            q_[0, :] += p_.get()[0] * self.qdt / masses
            q_[0, :] *= half_scale
            p_[0, :] *= 1.0 / half_scale
            self.cell.h_ *= half_scale
        else:
            half_scale = self._scr_state["half_scale"]

            q_[0, :] *= half_scale
            p_[0, :] *= 1.0 / half_scale
            self.cell.h_ *= half_scale
            q_[0, :] += p_.get()[0] * self.qdt / masses
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


dproperties(BaroSCR, ["bulk_modulus"])


class BaroBZP(Barostat):
    """Bussi-Zykova-Parrinello barostat class.

    Mode "isotropic" implements the Bussi-Zykova-Parrinello barostat, which
    scales the cell isotropically, changing its volume but not its shape, while
    sampling the isothermal-isobaric ensemble.

    Bussi, Zykova-Timan, Parrinello, J. Chem. Phys. 130, 074101 (2009)
        [doi:10.1063/1.3073889]
    Implementation details: Ceriotti, More, Manolopoulos, Comp. Phys. Comm.
        185, 1019 (2014) [doi:10.1016/j.cpc.2013.10.027]

    With an estimate of the bulk modulus B, the mass of the piston is
    9 V B (tau / 2 pi)^2, with V the reference volume v0, which defaults to the
    volume at the first bind, so that the volume oscillates with a period of
    about tau for any material and temperature. Otherwise it is 3 N kB T tau^2,
    as in Bussi et al., and the period is 2 pi tau / sqrt(3 B / (rho kB T)) for
    a material of bulk modulus B and atomic density rho.

    Depend objects:
       p: The momentum associated with the volume degree of freedom.
       m: The mass associated with the volume degree of freedom.
       bulk_modulus: The estimate of the bulk modulus, -1 for an ideal gas.
       v0: The reference volume that sets the mass with a given bulk modulus.
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
        bulk_modulus=None,
        v0=None,
    ):
        """Initializes BZP barostat.

        Args:
           dt, temp, tau, ebaro, thermostat: See `Barostat`.
           pext: Optional float giving the external pressure.
           p: Optional initial volume conjugate momentum. Defaults to 0.
           bulk_modulus: Optional estimate of the bulk modulus, which sets the
              mass of the piston. Defaults to that of an ideal gas.
           v0: Optional reference volume for the mass of the piston. Defaults
              to the volume of the cell at bind.
        """

        super(BaroBZP, self).__init__(dt, temp, tau, ebaro, thermostat)

        self.bulk_modulus_ = depend_value(
            name="bulk_modulus",
            value=-1.0 if bulk_modulus is None else bulk_modulus,
        )
        self.v0_ = depend_value(name="v0", value=-1.0 if v0 is None else v0)
        self.p_ = depend_array(name="p", value=xp.zeros(1))

        if p is not None:
            self.p = xp.asarray([p])
        else:
            self.p = 0.0

        if pext is not None:
            self.pext = pext

    def bind(self, beads, nm, cell, forces, bias=None, prng=None, nmts=1):
        """Binds beads, cell and forces to the barostat (see `Barostat.bind`)."""

        super(BaroBZP, self).bind(beads, nm, cell, forces, bias, prng, nmts)

        if self.v0 <= 0.0:
            self.v0 = float(self.cell.V)
        self.m_ = depend_array(
            name="m",
            value=xp.zeros(1),
            func=self.get_m,
            dependencies=[self.tau_, self.temp_, self.bulk_modulus_, self.v0_],
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
                self.cell_jacobian_,
                self.thermostat.ethermo_,
            ],
        )

    def get_m(self):
        """Calculates the mass of the piston (see the class)."""

        # the force on the piston is a sum over the beads, so V B carries a
        # factor nbeads, as does the temperature of the barostat
        if self.bulk_modulus > 0.0:
            vb = self.beads.nbeads * self.v0 * self.bulk_modulus
            return xp.asarray([9.0 * vb * (self.tau / (2.0 * np.pi)) ** 2])
        nkt = self.beads.natoms * Constants.kb * self.temp
        return xp.asarray([3.0 * nkt * self.tau**2])

    def get_pot(self):
        """Calculates the P V term of the enthalpy, for the nbeads replicas of
        the cell."""

        return self.cell.V * self.pext * self.beads.nbeads

    def get_cell_jacobian(self):
        """Calculates the energy term that accounts for the size of the box"""

        return -1.0 * xp.log(self.cell.V) * Constants.kb * self.temp

    def get_ebaro(self):
        """Calculates the barostat conserved quantity."""

        return self.thermostat.ethermo + self.kin + self.pot + self.cell_jacobian

    def pstep(self, level=0):
        """Propagates the momentum of the barostat."""

        # half a time step at the MTS level
        dt = self.pdt[level]

        # computes the pressure associated with the forces at each MTS level.
        press = xp.linalg.trace(self.stress_mts(level)) / 3.0
        self.p_ += dt * 3.0 * (self.cell.V * press)

        # integrates the kinetic part of the pressure with the force at the inner-most level.
        if level == self.nmtslevels - 1:
            nbeads = self.beads.nbeads
            self.p_ += (
                3.0
                * dt
                * (Constants.kb * self.temp - self.cell.V * (nbeads * self.pext))
            )

            pc = self.beads.pc
            fc = self.kick_force(level).sum(axis=0) / nbeads
            fc_m = fc / self.beads.m3[0]

            dt2 = dt * dt * nbeads
            dt3 = dt / 3.0
            self.p_ += ((pc @ fc_m) + dt3 * (fc @ fc_m)) * dt2

    def qcstep(self):
        """Propagates the centroid position and momentum and the volume."""

        v = self.p[0] / self.m[0]
        # half the time step of the inner MTS loop
        halfdt = self.qdt
        expq, expp = (xp.exp(v * halfdt), xp.exp(-v * halfdt))

        m = self.beads.m3[0]

        q_, p_ = self.centroid_arrays()
        q, p = q_.get()[0], p_.get()[0]
        q_[0, :] = (q * expq) + ((expq - expp) / (2.0 * v)) * (p / m)
        p_[0, :] = expp * p

        self.cell.h_ *= expq


dproperties(BaroBZP, ["p", "m", "bulk_modulus", "v0"])


class BaroSCBZP(BaroBZP):
    """The Suzuki Chin Bussi-Zykova-Parrinello barostat class.

    Mode "sc-isotropic" implements the Bussi-Zykova-Parrinello barostat (see
    `BaroBZP`) for Suzuki-Chin path integral molecular dynamics, with the
    'scnpt' dynamics.

    Kapil et al., J. Chem. Theory Comput. 15, 5 (2019)
        [10.1021/acs.jctc.8b01297]

    The forces and virials of each MTS level carry the +-1/3 Suzuki-Chin
    corrections, and the momentum is also propagated with the high-order
    part of the stress (`pscstep`).

    Depend objects (besides those of `BaroBZP`):
       kstress_sc: The high-order part of the kinetic stress.
       stress_sc: The high-order part of the stress.
    """

    def bind(self, beads, nm, cell, forces, bias=None, prng=None, nmts=1):
        """Binds beads, cell and forces to the barostat (see `BaroBZP.bind`)."""

        super().bind(beads, nm, cell, forces, bias, prng, nmts)

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

    def mts_force(self, level):
        """Returns the force at the MTS level `level`, with the +-1/3
        Suzuki-Chin correction."""

        return self.forces.mts_forces[level].f * (1 + self.forces.coeffsc_part_1)

    def mts_virial(self, level):
        """Returns the virial at the MTS level `level`, with the +-1/3
        Suzuki-Chin correction."""

        return (
            self.forces.mts_forces[level].virs
            * (1 + self.forces.coeffsc_part_1).reshape((self.beads.nbeads, 1, 1))
        ).sum(axis=0)

    def kstress_mts(self, level):
        """Calculates the Suzuki-Chin quantum centroid virial kinetic stress
        tensor associated with the forces at a MTS level."""

        return xp.triu(super().kstress_mts(level))

    def get_kstress_sc(self):
        """Calculates the high order part of the Suzuki-Chin
        quantum centroid virial kinetic stress tensor
        associated with the forces at a MTS level.
        """

        qqc = (self.beads.q - self.beads.qc).reshape(-1, 3)
        return xp.triu(-(qqc.T @ self.forces.fsc_part_2.reshape(-1, 3)))

    def get_stress_sc(self):
        """Calculates the high order part of the Suzuki-Chin internal stress tensor."""
        return (self.kstress_sc + self.forces.virssc_part_2.sum(axis=0)) / self.cell.V

    def pscstep(self):
        """Propagates the momentum of the barostat with respect to the
        high order part of the Suzuki-Chin stress"""

        # integrates with respect to the "high order" part of the stress with a timestep of dt /2
        press = xp.linalg.trace(self.stress_sc) / 3.0
        self.p_ += self.dt / 2 * 3.0 * (self.cell.V * press)


dproperties(BaroSCBZP, ["kstress_sc", "stress_sc"])


class BaroFlexibleCell(Barostat):
    """Base class of the barostats that change the whole cell (lengths and
    angles), RGB and MTK.

    The momentum of the cell is an upper-triangular 3x3 tensor, also accessed
    as a 6-vector of its independent elements, and the elements of the cell
    fixed by `hfix` do not change. The subclasses define the elastic energy of
    the cell (`get_pot`, which depends on `pot_dependencies`) and the external
    stress (`external_stress`).

    The mass of each element of the cell momentum is a third of the mass of the
    isotropic piston, 3 V B (tau / 2 pi)^2 or N kB T tau^2 (as in Martyna,
    Tobias, Klein, J. Chem. Phys. 101, 4177 (1994) [doi:10.1063/1.467468]), so
    that a uniform dilation of the cell oscillates as the volume of `BaroBZP`,
    with a period of about tau with an estimate of B. The other modes oscillate
    with frequencies set by the corresponding elastic constants instead of 3 B:
    C11 - C12 for a tetragonal and C44 for a shear deformation of a cubic
    crystal.

    The cell moves as dh/dt = v h, with v = p / m upper triangular and zero in
    the elements fixed by `hfix`, so the dynamics conserves the measure of the
    cells that is invariant under these transformations applied from the left.
    The force kB T L on the diagonal of p turns it into V times the measure
    invariant under the transformations applied from the right, i.e. under a
    change of the vectors that describe the same lattice, which for the full
    cell is the ensemble of Martyna, Tobias, Klein (1994). L_kk is one plus the
    number of free elements in row k minus that in column k, e.g. (3, 1, -1)
    for the full cell and (1, 1, 1) for a diagonal one. The free elements must
    form a group, so xz cannot be fixed while xy and yz change.

    Depend objects:
       p: The momentum matrix associated with the cell degrees of freedom.
       p6: The independent elements of p: xx, yy, zz, xy, xz, yz.
       m: The mass associated with the cell degree of freedom.
       m6: The masses of the elements of p6.
       hmask: The mask of the elements of the cell that can change.
       bulk_modulus: The estimate of the bulk modulus, -1 for an ideal gas.
       v0: The reference volume that sets the mass with a given bulk modulus.
    """

    def __init__(
        self,
        dt=None,
        temp=None,
        tau=None,
        ebaro=None,
        thermostat=None,
        vol_constraint=None,
        p=None,
        hfix=None,
        bulk_modulus=None,
        v0=None,
    ):
        """Initializes the barostat (see `Barostat`).

        Args:
           p: Optional initial cell momentum. Defaults to 0.
           hfix: Optional list of the elements of the cell that are kept fixed
              (see `mask_from_fix`).
           bulk_modulus: Optional estimate of the bulk modulus, which sets the
              mass of the cell. Defaults to that of an ideal gas.
           v0: Optional reference volume for the mass of the cell. Defaults to
              `reference_volume()` at bind.
        """

        super().__init__(dt, temp, tau, ebaro, thermostat, vol_constraint)

        self.bulk_modulus_ = depend_value(
            name="bulk_modulus",
            value=-1.0 if bulk_modulus is None else bulk_modulus,
        )
        self.v0_ = depend_value(name="v0", value=-1.0 if v0 is None else v0)

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
        free = np.triu(hmask)
        if free[0, 1] and free[1, 2] and not free[0, 2]:
            raise ValueError(
                "Fixing the cell element xz requires fixing xy or yz too, as xz "
                "changes when both xy and yz do."
            )

        # mask to zero out components of the cell velocity, to implement cell-boundary constraints
        self.hmask_ = depend_array(name="hmask", value=hmask)
        # the measure term of the free diagonal elements (see the class)
        self.L = xp.asarray(
            np.diag((free.sum(axis=1) - free.sum(axis=0) + 1) * np.diag(free))
        )

    def bind(self, beads, nm, cell, forces, bias=None, prng=None, nmts=1):
        """Binds beads, cell and forces to the barostat (see `Barostat.bind`)."""

        super().bind(beads, nm, cell, forces, bias, prng, nmts)

        if self.v0 <= 0.0:
            self.v0 = self.reference_volume()
        self.m_ = depend_array(
            name="m",
            value=xp.zeros(1),
            func=self.get_m,
            dependencies=[self.tau_, self.temp_, self.bulk_modulus_, self.v0_],
        )

        self.m6_ = depend_array(
            name="m6",
            value=xp.zeros(6),
            func=(lambda: xp.ones(6) * self.m[0]),
            dependencies=[self.m_],
        )

        # the elastic energy of the cell
        self.pot_ = depend_value(
            name="pot", func=self.get_pot, dependencies=self.pot_dependencies()
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

    def reference_volume(self):
        """Returns the default reference volume for the mass, that of the cell."""

        return float(self.cell.V)

    def get_m(self):
        """Calculates the mass of the elements of the cell momentum, a third of
        that of `BaroBZP.get_m` (see the class)."""

        if self.bulk_modulus > 0.0:
            vb = self.beads.nbeads * self.v0 * self.bulk_modulus
            return xp.asarray([3.0 * vb * (self.tau / (2.0 * np.pi)) ** 2])
        nkt = self.beads.natoms * Constants.kb * self.temp
        return xp.asarray([nkt * self.tau**2])

    def get_3x3to6(self):
        return self.select6 @ self.p.reshape(9)

    def get_6to3x3(self):
        return (self.p6 @ self.select6).reshape(3, 3)

    def get_cell_jacobian(self):
        """Calculates the energy term that accounts for the size of the box."""

        cj = sum(self.L[i, i] * xp.log(self.cell.h[i, i]) for i in range(3))
        return -1.0 * Constants.kb * self.temp * cj

    def get_ebaro(self):
        """Calculates the barostat conserved quantity."""

        return self.thermostat.ethermo + self.kin + self.pot + self.cell_jacobian

    def kstress_mts(self, level):
        """Calculates the kinetic stress tensor associated with the forces at a
        MTS level, with all the elements that the cell velocity v exchanges with
        the beads as q_c -> exp(v dt) q_c and p_c -> exp(-v dt) p_c: the sum of
        -f (q - q_c)^T over the beads and atoms and, at the inner-most level,
        nbeads p_c p_c^T / m over the atoms.
        """

        beads = self.beads
        qqc = (beads.q - beads.qc).reshape(-1, 3)
        kst = -(self.kick_force(level).reshape(-1, 3).T @ qqc)

        if level == self.nmtslevels - 1:
            pc = beads.pc.reshape(-1, 3)
            kst = kst + (beads.nbeads * pc.T / beads.m) @ pc

        return kst

    def pstep(self, level=0):
        """Propagates the momenta for half a time step."""

        dt = self.pdt[level]
        dt2 = dt**2
        dt3 = dt**3 / 3.0

        stress = self.stress_mts(level)
        dp = dt * (self.cell.V * xp.triu(stress))

        # integerates the kinetic part of the stress with the force at the inner-most level.
        if level == self.nmtslevels - 1:
            dp = dp + dt * (
                self.cell.V * xp.triu(-self.beads.nbeads * self.external_stress())
                + Constants.kb * self.temp * self.L
            )

            pc = self.beads.pc.reshape(self.beads.natoms, 3)
            fc = (
                self.kick_force(level).sum(axis=0).reshape(self.beads.natoms, 3)
                / self.beads.nbeads
            )
            fcTonm = (fc / self.beads.m3[0].reshape(self.beads.natoms, 3)).T

            # the kinetic stress integrated over the kick p -> p + f t
            fp = fcTonm @ pc
            dp = dp + (
                xp.triu(dt2 * 0.5 * (fp + fp.T) + dt3 * (fcTonm @ fc))
                * self.beads.nbeads
            )

        # the fixed cell elements, and the volume if it is constrained, receive
        # no momentum from the stress
        dp = dp * self.hmask
        if self.vol_constraint is True:
            dp = dp - xp.eye(3) * xp.linalg.trace(dp) / 3.0
        self.p_ += dp

        # applies the same constraints to the momentum from the thermostat, and
        # accumulates the change of kinetic energy in the energy of the
        # thermostat, so that the conserved quantity is unchanged
        if len(self.hfix) > 0 or self.vol_constraint is True:
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

        q_, p_ = self.centroid_arrays()
        q = q_.get()[0].reshape(self.beads.natoms, 3)
        p = p_.get()[0]

        q = q @ expq.T
        q += (p / self.beads.m3[0]).reshape(self.beads.natoms, 3) @ sinh.T
        p = p.reshape(self.beads.natoms, 3) @ expp.T

        q_[0] = q.reshape(self.beads.natoms * 3)
        p_[0] = p.reshape(self.beads.natoms * 3)

        self.cell.h = expq @ self.cell.h


dproperties(BaroFlexibleCell, ["p", "m", "p6", "hmask", "m6", "bulk_modulus", "v0"])


class BaroRGB(BaroFlexibleCell):
    """Raiteri-Gale-Bussi constant stress barostat class.

    Mode "anisotropic" implements the Raiteri-Gale-Bussi barostat, which
    changes the lengths and angles of the cell at a constant external stress,
    for solids under (non-diagonal) stresses. The elastic energy is that of the
    strain of the cell relative to a reference cell h0, so the ensemble is valid
    within the elastic limit of small strains. For diagonal stresses (or
    external pressures) and small strains, the "flexible" and "anisotropic"
    modes give very similar results.

    Raiteri, Gale, Bussi, J. Phys. Condens. Matter, 23, 334213 (2011)
        [doi:10.1088/0953-8984/23/33/334213]
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
        bulk_modulus=None,
    ):
        """Initializes RGB barostat.

        Args:
           dt, temp, tau, ebaro, thermostat: See `Barostat`.
           stressext: Optional 3x3 external stress.
           h0: The reference cell for the strain, which must be set before bind.
           p, hfix: See `BaroFlexibleCell`.
           bulk_modulus: Optional estimate of the bulk modulus, which sets the
              mass of the cell with the volume of the reference cell.
        """

        super().__init__(
            dt, temp, tau, ebaro, thermostat, p=p, hfix=hfix, bulk_modulus=bulk_modulus
        )

        if h0 is not None:
            self.h0 = h0
        else:
            self.h0 = Cell()

        if stressext is not None:
            self.stressext = stressext

    def reference_volume(self):
        """Returns the reference volume for the mass, that of the reference cell."""

        return float(self.h0.V)

    def pot_dependencies(self):
        """Returns the depend objects that the elastic energy depends on."""

        return [self.cell.h_, self.h0.h_, self.h0.V_, self.h0.ih_, self.stressext_]

    def get_pot(self):
        """Calculates the elastic energy V0 Tr(stress eps) of the Lagrangian
        strain eps relative to the reference cell, for the nbeads replicas of
        the cell."""

        eps = self.cell.h @ self.h0.ih
        eps = eps.T @ eps
        eps = 0.5 * (eps - xp.eye(3))

        return self.h0.V * xp.linalg.trace(self.stressext @ eps) * self.beads.nbeads

    def external_stress(self):
        """Returns the external stress, transformed from the reference cell to
        the current one."""

        hh0 = self.cell.h @ self.h0.ih
        return hh0 @ self.stressext @ hh0.T * self.h0.V / self.cell.V


class BaroMTK(BaroFlexibleCell):
    """Martyna-Tuckerman-Tobias-Klein flexible cell, constant pressure barostat
    class.

    Mode "flexible" implements the path integral version of the
    Martyna-Tuckerman-Tobias-Klein barostat, which changes the lengths and
    angles of the cell while sampling the isothermal-isobaric ensemble, for
    anisotropic systems such as molecular solids. With `vol_constraint` the
    volume is constant, sampling the (N, V, sigma_a = 0, T) ensemble of Rogge
    et al., J. Chem. Theory Comput. 11, 5583 (2015) [doi:10.1021/acs.jctc.5b00748].

    Martyna, Hughes, Tuckerman, J. Chem. Phys. 110, 3275 (1999)
        [doi:10.1063/1.478193]
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
        bulk_modulus=None,
        v0=None,
    ):
        """Initializes MTK barostat.

        Args:
           dt, temp, tau, ebaro, thermostat, vol_constraint: See `Barostat`.
           pext: Optional float giving the external pressure.
           p, hfix, bulk_modulus, v0: See `BaroFlexibleCell`.
        """

        super().__init__(
            dt,
            temp,
            tau,
            ebaro,
            thermostat,
            vol_constraint,
            p,
            hfix,
            bulk_modulus,
            v0,
        )

        if pext is not None:
            self.pext = pext

    def pot_dependencies(self):
        """Returns the depend objects that the elastic energy depends on."""

        return [self.cell.h_, self.pext_]

    def get_pot(self):
        """Calculates the P V term of the enthalpy, for the nbeads replicas of
        the cell."""

        return self.cell.V * self.pext * self.beads.nbeads

    def external_stress(self):
        """Returns the external stress, an isotropic pressure."""

        return xp.eye(3) * self.pext

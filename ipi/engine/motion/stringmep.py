"""Contains the algorithms to perform String minimal energy path calculations.

(C) Karen Fidanyan 2022

This file is part of i-PI.
i-PI Copyright (C) 2014-2022 i-PI developers
See the "licenses" directory for full license information.
"""

import numpy as np
from numpy.linalg import norm as npnorm
import time

from ipi.engine.motion import Motion
from ipi.utils.depend import dstrip
from ipi.utils.softexit import softexit
from ipi.utils.mintools import Damped_BFGS, BFGSTRM, FIRE
from ipi.utils.messages import verbosity as vrb, info  # , warning
from ipi.engine.beads import Beads

try:
    import scipy
    from scipy.interpolate import make_interp_spline, splev
except Exception as e:
    scipy = None
    scipy_exception = e

np.set_printoptions(threshold=10000, linewidth=1000)  # Remove in cleanup

decor = 60 * "=" + "\n"

__all__ = [
    "StringGradientMapper",
    "GradientMapper",
    "StringClimbGrMapper",
    "StringMover",
]


class StringGradientMapper(object):
    """Creation of the multi-dimensional function that will be minimized.
    Functional analog of a GradientMapper in geop.py,
    but since the endpoints don't move, we construct another
    object 'rbeads' with (N-2) beads to avoid constant recalculation
    of the forces on the end-point beads.

    This version returns the orthogonal component of the force, as required in the
    original 2002 version of the string method.

    Fixed atoms are excluded via boolean mask. 1 = moving, 0 = fixed.
    """

    def __init__(self):
        self.allpots = None

    def bind(self, ens):
        if scipy is None:
            info(" @NM: scipy import failed", vrb.low)
            raise scipy_exception

        # In principle, there is no need in dforces within the Mapper,
        # BUT dbeads are needed to calculate tangents for the endpoints,
        # and dforces are needed outside the Mapper to construct the "main" forces.
        self.dbeads = ens.beads.copy()
        self.dcell = ens.cell.copy()
        self.dforces = ens.forces.copy(self.dbeads, self.dcell)
        self.fixatoms = ens.fixatoms.copy()

        # Mask to exclude fixed atoms from 3N-arrays
        self.fixmask = np.ones(3 * ens.beads.natoms, dtype=bool)
        if len(ens.fixatoms) > 0:
            self.fixmask[3 * ens.fixatoms] = 0
            self.fixmask[3 * ens.fixatoms + 1] = 0
            self.fixmask[3 * ens.fixatoms + 2] = 0

        # Create reduced bead and force object
        self.rbeads = Beads(ens.beads.natoms, ens.beads.nbeads - 2)
        self.rbeads.q[:] = ens.beads.q[1:-1]
        self.rforces = ens.forces.copy(self.rbeads, self.dcell)

    def __call__(self, x):
        """Returns the potential for all beads and the gradient."""

        # Bead positions
        # Touch positions only if they have changed (to avoid triggering forces)
        # I need both dbeads and rbeads because of the endpoint tangents.
        if (self.rbeads.q[:, self.fixmask] != x).any():
            self.rbeads.q[:, self.fixmask] = x
        rbq = self.rbeads.q[:, self.fixmask]
        # We need all beads, but only reduced number is being updated.
        self.dbeads.q[1:-1, self.fixmask] = rbq

        # Forces and energies
        rbf = dstrip(self.rforces.f).copy()[:, self.fixmask]
        be = dstrip(self.rforces.pots).copy()

        # Calculate the force component perpendicular to the spline.
        tangent = spline_derv(self.dbeads.q[:, self.fixmask], self.dbeads.nbeads)[1:-1]
        rbf_perp = rbf - tangent * np.sum(rbf * tangent, axis=1)[:, np.newaxis]

        # Return potential energy of the whole string and minus f_perpendicular
        e = be.sum()
        g = -rbf_perp
        return e, g


class GradientMapper(object):
    """Creation of the multi-dimensional function that will be minimized.
    Functional analog of a GradientMapper in geop.py,
    but since the endpoints don't move, we construct another
    object 'rbeads' with (N-2) beads to avoid constant recalculation
    of the forces on the end-point beads.

    This version does not make any projections and returns just physical forces
    for the intermediate beads (except fixed atoms). Its purpose here is to be used
    in the "Simplified and improved string method" (JCP 126, 164103 (2007),
    http://dx.doi.org/10.1063/1.2720838), which is not yet implemented.

    Fixed atoms are excluded via boolean mask. 1 = moving, 0 = fixed.
    """

    def __init__(self):
        self.allpots = None

    def bind(self, ens):
        if scipy is None:
            info(" @NM: scipy import failed", vrb.low)
            raise scipy_exception

        # In principle, there is no need in dforces within the Mapper,
        # BUT dbeads are needed to calculate tangents for the endpoints,
        # and dforces are needed outside the Mapper to construct the "main" forces.
        self.dbeads = ens.beads.copy()
        self.dcell = ens.cell.copy()
        self.dforces = ens.forces.copy(self.dbeads, self.dcell)
        self.fixatoms = ens.fixatoms.copy()

        # Mask to exclude fixed atoms from 3N-arrays
        self.fixmask = np.ones(3 * ens.beads.natoms, dtype=bool)
        if len(ens.fixatoms) > 0:
            self.fixmask[3 * ens.fixatoms] = 0
            self.fixmask[3 * ens.fixatoms + 1] = 0
            self.fixmask[3 * ens.fixatoms + 2] = 0

        # Create reduced bead and force object
        self.rbeads = Beads(ens.beads.natoms, ens.beads.nbeads - 2)
        self.rbeads.q[:] = ens.beads.q[1:-1]
        self.rforces = ens.forces.copy(self.rbeads, self.dcell)

    def __call__(self, x):
        """Returns the potential for all beads and the gradient."""

        # Bead positions
        # Touch positions only if they have changed (to avoid triggering forces)
        # I need both dbeads and rbeads because of the endpoint tangents.
        if (self.rbeads.q[:, self.fixmask] != x).any():
            self.rbeads.q[:, self.fixmask] = x
        rbq = self.rbeads.q[:, self.fixmask]
        # We need all beads, but only reduced number is being updated.
        self.dbeads.q[1:-1, self.fixmask] = rbq

        # Forces
        rbf = dstrip(self.rforces.f).copy()[:, self.fixmask]
        be = dstrip(self.rforces.pots).copy()

        # Return potential energy of the whole string and full force gradient
        e = be.sum()
        g = -rbf
        return e, g


class StringClimbGrMapper(object):
    """Creation of the multi-dimensional function that will be minimized.
    Functional analog of a GradientMapper in geop.py

    Constructs climbing forces for a single node.
    Node index is determined inside bind()
    Fixed atoms are excluded via boolean mask. 1 = moving, 0 = fixed.
    """

    def __init__(self):
        pass

    def bind(self, ens):
        if scipy is None:
            info(" @NM: scipy import failed", vrb.low)
            raise scipy_exception

        """Creates reduced Beads object in order to calculate forces
        only for the climbing bead, and binds it to beads
        """
        self.fixatoms = ens.fixatoms.copy()
        # A mask to exclude fixed atoms from 3N-arrays
        self.fixmask = np.ones(3 * ens.beads.natoms, dtype=bool)
        if len(ens.fixatoms) > 0:
            self.fixmask[3 * ens.fixatoms] = 0
            self.fixmask[3 * ens.fixatoms + 1] = 0
            self.fixmask[3 * ens.fixatoms + 2] = 0

        # Reduced Beads object is needed to calculate only required beads.
        self.rbeads = Beads(ens.beads.natoms, 1)
        self.rcell = ens.cell.copy()
        # Coords of the bead before and after the climbing one.
        self.q_prev = np.zeros(3 * (ens.beads.natoms - len(ens.fixatoms)))
        self.q_next = np.zeros(3 * (ens.beads.natoms - len(ens.fixatoms)))
        # Make reduced forces dependent on reduced beads
        self.rforces = ens.forces.copy(self.rbeads, self.rcell)

    def __call__(self, x):
        """Returns climbing force for a climbing image."""
        if self.q_prev is None or self.q_next is None:
            raise RuntimeError("StringClimbGrMapper.q_prev or q_next is None.")

        info(
            "@STRING_CLIMB: StringClimbGrMapper.x.shape: %s" % str(x.shape),
            vrb.debug,
        )
        # Touch positions only if they have changed (to avoid triggering forces)
        if (self.rbeads.q[:, self.fixmask] != x).any():
            self.rbeads.q[:, self.fixmask] = x
        # Compared to stringgradientMapper, I need to flatten here,
        # otherwise it's (1, 3N) instead of just 3N.
        rbq = self.rbeads.q[:, self.fixmask].flatten()

        # Reduced forces
        # We don't need copy() because flatten() implies copying.
        rbf = dstrip(self.rforces.f)[:, self.fixmask].flatten()

        # I think here it's better to use plain tangents.
        # Then we don't need energies of the neighboring beads.
        d1 = rbq - self.q_prev  # tau minus
        d2 = self.q_next - rbq  # tau plus
        tau = d1 / npnorm(d1) + d2 / npnorm(d2)
        tau /= npnorm(tau)

        # Inverting the force component along the path
        rbf -= 2 * np.dot(rbf, tau) * tau

        # Return the potential energy and gradient of the climbing bead
        e = self.rforces.pot
        g = -rbf
        info("@STRING_CLIMB: StringClimbGrMapper finished.", vrb.debug)
        return e, g


class StringMover(Motion):
    """MEP optimization routine for the String method.
    The algorithm is close to the one described in
    W. E, W. Ren, E. Vanden-Eijnden, PRB 66, 052301 (2002)
    "String method for the study of rare events".

    Attributes:
        mode: minimizer to use for string optimization.
        biggest_step: maximum step size for BFGS family algorithms
        old_force: force from previous iteration (not projected)
        old_direction: direction from previous iteration
        old_stringpotential: bead potentials from previous iteration
        old_stringgradient: gradient for the string from previous iteration
        hessian_bfgs: Hessian for (damped) BFGS
        tolerances:
            energy: tolerance on change in energy for exiting minimization
            force: tolerance on force/change in force for exiting minimization
            position: tolerance and change in position for exiting minimization
        corrections_lbfgs: number of corrections to store for L-BFGS
        qlist_lbfgs: list of previous positions (x_n+1 - x_n) for L-BFGS
        glist_lbfgs: list of previous gradients (g_n+1 - g_n) for L-BFGS
        endpoints: flag for minimizing end images of the string *** NOT YET IMPLEMENTED ***
        use_climb: flag for climbing image routine
        stage: flag denoting current stage: "endpoints", "string" or "climb", or "converged"
    """

    def __init__(
        self,
        fixcom=False,
        fixatoms=None,
        mode="damped_bfgs",
        biggest_step=0.5,
        old_coord=np.zeros(0, float),
        full_force=np.zeros(0, float),
        full_pots=np.zeros(0, float),
        old_stringpotential=np.zeros(0, float),
        old_stringgradient=np.zeros(0, float),
        old_direction=np.zeros(0, float),
        hessian_bfgs=np.eye(0),
        tolerances={"energy": 1e-5, "force": 1e-5, "position": 1e-5},
        corrections_lbfgs=5,
        qlist_lbfgs=np.zeros(0, float),
        glist_lbfgs=np.zeros(0, float),
        tr_trm=[0.1],
        dtmax_fire=1.0,
        v_fire=np.zeros(0, float),
        alpha_fire=0.1,
        N_down_fire=0,
        N_up_fire=0,
        dt_fire=0.1,
        endpoints=(False, "bfgs"),
        scale_lbfgs=2,
        stage="string",
        use_climb=False,
        climb_bead=-1,
    ):
        """Initialises StringMover.

        Args:
           fixcom: An optional boolean which decides whether the centre of mass
              motion will be constrained or not. Defaults to False.
        """

        super(StringMover, self).__init__(fixcom=fixcom, fixatoms=fixatoms)

        # Optimization options
        self.tolerances = tolerances
        self.mode = mode
        self.big_step = biggest_step
        if self.big_step > 0.5:
            info(
                " @DampedBFGS: WARNING: big_step is quite big. "
                "Damped_BFGS handles it differently from other algorithms, "
                "see ipi/inputs/motion/neb.py:'biggest_step' to get an idea.\n"
                "Current big_step value: %f." % self.big_step,
                vrb.medium,
            )
        self.old_x = old_coord  # always in "masked" dimension
        self.stringpot = old_stringpotential
        self.stringgrad = old_stringgradient
        self.d = old_direction
        self.full_f = full_force  # physical forces of ALL beads even in climb stage
        self.full_v = full_pots  # potentials of ALL beads
        self.corrections = corrections_lbfgs
        self.qlist = qlist_lbfgs
        self.glist = glist_lbfgs
        self.scale = scale_lbfgs
        self.tr_trm = tr_trm
        self.endpoints = endpoints
        self.use_climb = use_climb
        self.cl_indx = climb_bead
        self.stage = stage
        # damped_bfgs
        self.hessian = hessian_bfgs
        # fire
        self.v = v_fire  # velocity
        self.a = alpha_fire  # alpha. velocity mixing factor
        self.N_dn = N_down_fire  # consecutive steps in downhill dierction
        self.N_up = N_up_fire  # consecutive steps in uphill dierction
        self.dt_fire = dt_fire
        self.dtmax = dtmax_fire  # maximum dt for FIRE

        self.stringgm = StringGradientMapper()
        self.climbgm = StringClimbGrMapper()

    def bind(self, ens, beads, nm, cell, bforce, prng, omaker):
        if scipy is None:
            info(" @NM: scipy import failed", vrb.low)
            raise scipy_exception

        super(StringMover, self).bind(ens, beads, nm, cell, bforce, prng, omaker)
        # We check dimensionality of the Hessian in the beginning
        # and reduce it if needed, because someone may want to provide
        # existing Hessian of the full system.
        if self.mode in ["damped_bfgs", "bfgstrm"]:
            n_activedim = beads.q[0].size - len(self.fixatoms) * 3
            if self.stage == "string":
                if self.hessian.size == (n_activedim * (beads.nbeads - 2)) ** 2:
                    # Desired dimension
                    info(
                        " @STRINGMover: String Hessian is treated as reduced "
                        "and masked already, according to its size.",
                        vrb.high,
                    )
                elif self.hessian.size == beads.q.size * beads.q.size:
                    info(
                        " @STRINGMover: String Hessian has full-dimensional size, "
                        "will be reduced and masked by fixatoms.",
                        vrb.high,
                    )
                    # First trim the endpoints, then mask fixed atoms
                    self.hessian = (
                        self.hessian[
                            self.beads.natoms : -self.beads.natoms,
                            self.beads.natoms : -self.beads.natoms,
                        ]
                    )[
                        np.ix_(
                            np.tile(self.stringgm.fixmask, self.beads.nbeads - 2),
                            np.tile(self.stringgm.fixmask, self.beads.nbeads - 2),
                        )
                    ]
                elif self.hessian.size == 0:
                    info(
                        " @STRINGMover: No Hessian provided, starting from unity matrix",
                        vrb.high,
                    )
                    self.hessian = np.eye(
                        (self.beads.nbeads - 2) * n_activedim, dtype=float
                    )
                else:
                    raise ValueError("Hessian size does not match system size.")
            # In climb stage, hessian has only 1 bead.
            # Not sure whether fixatoms should be handled here or not, let's see...
            elif self.stage == "climb":
                if self.hessian.size == n_activedim * n_activedim:
                    # desired dimension
                    info(
                        " @STRINGMover: Hessian of the climber is treated as masked one,"
                        " according to its size.",
                        vrb.high,
                    )
                elif self.hessian.size == (beads.q[0].size * beads.q[0].size):
                    info(
                        " @STRINGMover: Hessian of the climber has size (3natoms x 3natoms), "
                        "will be masked by fixatoms.",
                        vrb.high,
                    )
                    self.hessian = self.hessian[
                        np.ix_(self.climbgm.fixmask, self.climbgm.fixmask)
                    ]
                if self.hessian.size == 0:
                    info(
                        " @STRINGMover: No Hessian provided for climbing, "
                        "starting from unity matrix",
                        vrb.high,
                    )
                    self.hessian = np.eye(n_activedim, dtype=float)
                else:
                    raise ValueError("Hessian size does not match system size.")

        if len(self.fixatoms) == len(self.beads[0]):
            softexit.trigger(
                status="bad",
                message="WARNING: all atoms are fixed, geometry won't change. Exiting simulation.",
            )

        self.stringgm.bind(self)
        self.climbgm.bind(self)

    def init_climb(self):
        """Operations to resize hessian for the climbing bead(s),
        get neighboring beads and other things necessary for climbing.
        Returns:
          cl_indx - index of the bead for climbing
        """
        cl_indx = np.argmax(self.forces.pots)
        info(
            " @STRING: Initializing climbing. Climbing bead: %i." % cl_indx,
            vrb.medium,
        )
        if cl_indx in [0, self.beads.nbeads - 1]:
            softexit.trigger(
                status="bad", message="ERROR: climbing bead is the endpoint."
            )
        self.climbgm.rbeads.q[:] = self.beads.q[cl_indx]
        self.climbgm.q_prev[:] = self.beads.q[cl_indx - 1, self.climbgm.fixmask]
        self.climbgm.q_next[:] = self.beads.q[cl_indx + 1, self.climbgm.fixmask]
        info(
            "q_prev.shape: %s,    q_next.shape: %s"
            % (str(self.climbgm.q_prev.shape), str(self.climbgm.q_next.shape)),
            vrb.debug,
        )

        info(" @STRING_CLIMB: call StringClimbGrMapper for the first time.", vrb.debug)
        self.stringpot, self.stringgrad = self.climbgm(
            self.beads.q[cl_indx, self.climbgm.fixmask]
        )

        if self.mode in ["damped_bfgs", "bfgstrm"]:
            # Initialize BFGS Hessian for a single bead
            n_activedim = self.beads.q[0].size - len(self.fixatoms) * 3
            if self.hessian.shape != (n_activedim, n_activedim):
                self.hessian = np.eye(n_activedim)
        elif self.mode == "fire":
            # Initialize FIRE parameters
            self.a = 0.1
            self.v = -self.a * self.stringgrad
            self.N_dn = 0
            self.N_up = 0
            self.dt_fire = 0.1

        self.old_x = dstrip(self.beads.q[cl_indx, self.climbgm.fixmask]).copy()
        self.stage = "climb"
        return cl_indx

    def path_step_textbook2002(self, step=None):
        """A variant of the step which is responsible for moving N-2 intermediate beads.
        It uses StringGradientMapper to construct N-2 beads object in order to avoid
        recalculation of the static endpoints at every step.

        'textbook' variant does optimizer step and spline resampling consecutively,
        as suggested in W. E, W. Ren, E. Vanden-Eijnden, PRB 66, 052301 (2002)
        https://doi.org/10.1103/PhysRevB.66.052301

        The main problem here is that any optimizer from 'mintools.py' needs to
        receive forces at resampled positions, but returns forces from NOT resampled
        ones, therefore it's necessary to have 2 force calls per step. We provide an
        alternative step 'path_step_single_f_call', which uses different optimization
        routines with spline resampling inside, but it's not strictly the textbook
        String method anymore.
        """
        self.ptime = self.ttime = 0
        self.qtime = -time.time()
        # Shortcuts for prettier expresions
        nbeads = self.beads.nbeads
        natoms = self.beads.natoms
        fixmask = self.stringgm.fixmask

        n_activedim = self.beads.q[0].size - len(self.fixatoms) * 3

        # Store 'old' resampled positions
        self.old_x[:] = self.beads.q[1:-1, fixmask]

        # 'Self' stringpot will be updated later on,
        # therefore we store the copy to check convergence in the end.
        old_stringpot = self.stringpot.copy()

        if self.mode in ["damped_bfgs", "bfgstrm"]:
            # All BFGS-family algorithms would have similar structure.
            # For simplicity, we only use the two which use direct Hessian
            # not the inverse one.
            if step == 0:  # TODO add a condition when after the endpoints.
                # With multiple stages, the size of the hessian is different
                # at each stage, therefore we check.
                if self.hessian.shape != (
                    (nbeads - 2) * n_activedim,
                    (nbeads - 2) * n_activedim,
                ):
                    print("Dimensions of the Hessian and of the beads:")
                    print((self.hessian.shape, self.beads.q.shape))
                    softexit.trigger(
                        status="bad",
                        message="Wrong Hessian size in String step.",
                    )

            if self.mode == "damped_bfgs":
                info(" @STRING: before Damped_BFGS() call", vrb.debug)
                info("self.old_x.shape: %s" % str(self.old_x.shape), vrb.debug)
                info(
                    "self.stringgrad.shape: %s" % str(self.stringgrad.shape), vrb.debug
                )
                info("self.hessian.shape: %s" % str(self.hessian.shape), vrb.debug)
                Damped_BFGS(
                    x0=self.old_x.copy(),
                    fdf=self.stringgm,
                    fdf0=(self.stringpot, self.stringgrad),
                    hessian=self.hessian,
                    big_step=self.big_step,
                )
                info(" @STRING: after Damped_BFGS() call", vrb.debug)

            elif self.mode == "bfgstrm":
                info(" @STRING: before BFGSTRM() call", vrb.debug)
                info("self.old_x.shape: %s" % str(self.old_x.shape), vrb.debug)
                info(
                    "self.stringgrad.shape: %s" % str(self.stringgrad.shape), vrb.debug
                )
                info("self.hessian.shape: %s" % str(self.hessian.shape), vrb.debug)
                info("self.tr_trm: %s" % self.tr_trm, vrb.debug)
                BFGSTRM(
                    x0=self.old_x.copy(),
                    u0=self.stringpot,
                    # BFGSTRM expects force instead of gradient, therefore minus
                    f0=-self.stringgrad,
                    h0=self.hessian,
                    tr=self.tr_trm,
                    mapper=self.stringgm,
                    big_step=self.big_step,
                )
                info(" @STRING: after BFGSTRM() call", vrb.debug)

        elif self.mode == "fire":
            # Only initialize velocity for fresh start, not for RESTART
            if step == 0 and self.v.size == 0:
                # Mapper has already been called in the step()
                self.v = -self.a * self.stringgrad

            info(" @STRING: calling FIRE", vrb.debug)
            info(" @FIRE velocity: %s" % str(npnorm(self.v)), vrb.debug)
            info(" @FIRE alpha: %s" % str(self.a), vrb.debug)
            info(" @FIRE N down: %s" % str(self.N_dn), vrb.debug)
            info(" @FIRE N up: %s" % str(self.N_up), vrb.debug)
            info(" @FIRE dt: %s" % str(self.dt_fire), vrb.debug)
            self.v, self.a, self.N_dn, self.N_up, self.dt_fire = FIRE(
                x0=self.old_x.copy(),
                fdf=self.stringgm,
                fdf0=(self.stringpot, self.stringgrad),
                v=self.v,
                a=self.a,
                N_dn=self.N_dn,
                N_up=self.N_up,
                dt=self.dt_fire,
                dtmax=self.dtmax,
            )
            info(" @STRING: after FIRE call", vrb.debug)

        elif self.mode == "euler":
            # Dummy step "dq/dt = f" with small dt.
            info(" @EULER: making a step.", vrb.debug)
            dt = 1e-3
            dq = -self.stringgrad * dt
            if np.any(np.abs(dq) > self.big_step):
                info(" @EULER: WARNING: hit the big_step.")
                dq *= self.big_step / np.amax(np.abs(dq))
            info(" @EULER: maxdispl is %g" % np.amax(np.abs(dq)))
            self.beads.q[1:-1, fixmask] += dq.reshape((nbeads - 2, -1))
            # info(" @EULER: calling the mapper.", vrb.debug)
            self.stringgm.dbeads.q[:] = self.beads.q[:]
            # self.stringpot, self.stringgrad = self.stringgm(self.beads.q[1:-1, fixmask])

        # TODO: Routines for L-BFGS, SD, CG
        else:
            print("Error: mode %s is not supported." % self.mode)
            softexit.trigger(
                status="bad",
                message="Try 'damped_bfgs', 'bfgstrm' or 'fire'. Other algorithms are not implemented for string optimization.",
            )

        # Update 'main' beads positions
        self.beads.q[:] = self.stringgm.dbeads.q
        info(" @STRING: bead positions transferred from mapper.", vrb.debug)

        # Resample the beads after an optimizer step.
        info(" @STRING: resampling beads uniformly.", vrb.debug)
        self.beads.q[:, fixmask] = spline_resample(
            self.beads.q[:, fixmask],
            nbeads,
        )

        # Here the forces are called on resampled positions.
        info(" @STRING: calling forces on resampled positions.", vrb.debug)
        self.stringpot, self.stringgrad = self.stringgm(self.beads.q[1:-1, fixmask])
        # This transfers forces from the mapper beads to the "main" beads,
        # so that recalculation won't be triggered after the step.
        # We need to keep forces up to date to be able to output
        # full potentials and full forces.
        tmp_v = self.full_v
        tmp_v[1:-1] = self.stringgm.rforces.pots
        tmp_f = self.full_f.copy()
        tmp_f[1:-1] = self.stringgm.rforces.f
        self.forces.transfer_forces_manual(
            new_q=[self.beads.q],
            new_v=[tmp_v],
            new_forces=[tmp_f],
        )
        # Full-dimensional forces are needed in reduced-beads modes,
        # i.e. in "endpoints" and "climb".
        self.full_f = dstrip(self.forces.f)
        self.full_v = dstrip(self.forces.pots)

        # Calculate the force component perpendicular to the spline to check convergence.
        tangent = spline_derv(self.beads.q[:, fixmask], nbeads)[1:-1]
        f_perp = (
            self.full_f[1:-1, fixmask]
            - tangent
            * np.sum(self.full_f[1:-1, fixmask] * tangent, axis=1)[:, np.newaxis]
        )

        # Check convergence at the resampled positions.
        # dx = current resampled position - previous resampled position.
        dx = np.amax(np.abs(self.beads.q[1:-1, fixmask] - self.old_x))
        de = np.amax(np.abs(self.stringpot - old_stringpot)) / (nbeads * natoms)
        if (
            (de <= self.tolerances["energy"])
            and (dx <= self.tolerances["position"])
            and (np.amax(np.abs(f_perp)) <= self.tolerances["force"])
        ):
            info(
                decor
                + " @STRING: path optimization converged. Step: %i\n" % step
                + decor,
                vrb.medium,
            )

            # Set climbing stage indicator after convergence of the path
            if self.use_climb:
                self.stage = "climb"
            else:
                self.stage = "converged"
                softexit.trigger(
                    status="success",
                    message="String MEP finished successfully at STEP %i." % step,
                )

        else:
            info(
                " @STRING: Not converged, deltaEnergy = %.8f, tol = %.8f per atom"
                % (de, self.tolerances["energy"]),
                vrb.debug,
            )
            info(
                " @STRING: Not converged, f_perp_2 = %.8f, tol = %f"
                % (np.amax(np.abs(f_perp)), self.tolerances["force"]),
                vrb.debug,
            )
            info(
                " @STRING: Not converged, deltaX = %.8f, tol = %.8f"
                % (dx, self.tolerances["position"]),
                vrb.debug,
            )

        self.qtime += time.time()

    def path_step_single_f_call(self, step=None):
        """A variant of the step which is responsible for moving N-2 intermediate beads.
        It uses StringGradientMapper to construct N-2 beads object in order to avoid
        recalculation of the static endpoints every step.

        'single_f_call' variant calls forces only once per step, because it
        uses custom optimization algorithms which resample the spline in the middle of
        an optimization step. Thus, it's not the textbook algorithm anymore, but it is
        TWICE more efficient in number of force calls, so we provide both variants.

        11.02.2022: doesn't work yet.
        """
        softexit.trigger(
            status="bad",
            message="Optimizers for path_step_single_f_call() are not yet implemented.",
        )
        # self.ptime = self.ttime = 0
        # self.qtime = -time.time()
        # # Shortcuts for prettier expresions
        # nbeads = self.beads.nbeads
        # natoms = self.beads.natoms
        # fixmask = self.stringgm.fixmask

        # n_activedim = self.beads.q[0].size - len(self.fixatoms) * 3

        # # Resample beads in the beginning of every step.
        # info(" @STRING: resampling beads uniformly.", vrb.debug)
        # self.beads.q[:, fixmask] = spline_resample(
        #     self.beads.q[:, fixmask],
        #     nbeads,
        # )
        # # Store 'old' resampled positions
        # self.old_x[:] = self.beads.q[1:-1, fixmask]

        # # 'Self' instances will be updated in the optimizer,
        # # therefore we store the copy of the potential.
        # # old_stringpot is used later as a convergence criterion.
        # old_stringpot = self.stringpot.copy()

        # if self.mode in ["damped_bfgs", "bfgstrm"]:
        #     # All BFGS-family algorithms would have similar structure.
        #     # For simplicity, we only use the two which use direct Hessian
        #     # not the inverse one.
        #     if step == 0:  # TODO add a condition when after the endpoints.
        #         # With multiple stages, the size of the hessian is different
        #         # at each stage, therefore we check.
        #         if self.hessian.shape != (
        #             (nbeads - 2) * n_activedim,
        #             (nbeads - 2) * n_activedim,
        #         ):
        #             print("Dimensions of the Hessian and of the beads:")
        #             print((self.hessian.shape, self.beads.q.shape))
        #             softexit.trigger(
        #                 status="bad",
        #                 message="Wrong Hessian size in String step.",
        #             )

        #     if self.mode == "damped_bfgs":
        #         info(" @STRING: before Damped_BFGS() call", vrb.debug)
        #         info("self.old_x.shape: %s" % str(self.old_x.shape), vrb.debug)
        #         info(
        #             "self.stringgrad.shape: %s" % str(self.stringgrad.shape), vrb.debug
        #         )
        #         info("self.hessian.shape: %s" % str(self.hessian.shape), vrb.debug)
        #         Damped_BFGS(
        #             x0=self.old_x.copy(),
        #             fdf=self.stringgm,
        #             fdf0=(self.stringpot, self.stringgrad),
        #             hessian=self.hessian,
        #             big_step=self.big_step,
        #         )
        #         info(" @STRING: after Damped_BFGS() call", vrb.debug)

        #     elif self.mode == "bfgstrm":
        #         info(" @STRING: before BFGSTRM() call", vrb.debug)
        #         info("self.old_x.shape: %s" % str(self.old_x.shape), vrb.debug)
        #         info(
        #             "self.stringgrad.shape: %s" % str(self.stringgrad.shape), vrb.debug
        #         )
        #         info("self.hessian.shape: %s" % str(self.hessian.shape), vrb.debug)
        #         BFGSTRM(
        #             x0=self.old_x.copy(),
        #             u0=self.stringpot,
        #             # BFGSTRM expects force instead of gradient, therefore minus
        #             f0=-self.stringgrad,
        #             h0=self.hessian,
        #             tr=self.tr_trm,
        #             mapper=self.stringgm,
        #             big_step=self.big_step,
        #         )
        #         info(" @STRING: after BFGSTRM() call", vrb.debug)

        # elif self.mode == "fire":
        #     # Only initialize velocity for fresh start, not for RESTART
        #     if step == 0 and self.v.size == 0:
        #         # Mapper has already been called in the step()
        #         self.v = -self.a * self.stringgrad

        #     info(" @STRING: calling FIRE", vrb.debug)
        #     info(" @FIRE velocity: %s" % str(npnorm(self.v)), vrb.debug)
        #     info(" @FIRE alpha: %s" % str(self.a), vrb.debug)
        #     info(" @FIRE N down: %s" % str(self.N_dn), vrb.debug)
        #     info(" @FIRE N up: %s" % str(self.N_up), vrb.debug)
        #     info(" @FIRE dt: %s" % str(self.dt_fire), vrb.debug)
        #     self.v, self.a, self.N_dn, self.N_up, self.dt_fire = FIRE(
        #         x0=self.old_x.copy(),
        #         fdf=self.stringgm,
        #         fdf0=(self.stringpot, self.stringgrad),
        #         v=self.v,
        #         a=self.a,
        #         N_dn=self.N_dn,
        #         N_up=self.N_up,
        #         dt=self.dt_fire,
        #         dtmax=self.dtmax,
        #     )
        #     info(" @STRING: after FIRE call", vrb.debug)

        # # TODO: Routines for L-BFGS, SD, CG
        # else:
        #     print("Error: mode %s is not supported." % self.mode)
        #     softexit.trigger(
        #         status="bad",
        #         message="Try 'damped_bfgs', 'bfgstrm' or 'fire'. Other algorithms are not implemented for string optimization.",
        #     )

        # # Update positions
        # self.beads.q[:] = self.stringgm.dbeads.q
        # info(" @STRING: bead positions transferred from mapper.", vrb.debug)

        # # Recalculation won't be triggered because the position is the same.
        # self.stringpot, self.stringgrad = self.stringgm(self.beads.q[1:-1, fixmask])
        # info(" @STRING: forces transferred from mapper.", vrb.debug)

        # # This transfers forces from the mapper to the "main" beads,
        # # so that recalculation won't be triggered after the step.
        # # We need to keep forces up to date to be able to output
        # # full potentials and full forces.
        # tmp_f = self.full_f.copy()
        # tmp_f[1:-1] = self.stringgm.rforces.f
        # tmp_v = self.full_v
        # tmp_v[1:-1] = self.stringgm.rforces.pots
        # self.forces.transfer_forces_manual(
        #     new_q=[
        #         self.beads.q,
        #     ],
        #     new_v=[
        #         tmp_v,
        #     ],
        #     new_forces=[
        #         tmp_f,
        #     ],
        # )

        # # Full-dimensional forces are needed in reduced-beads modes,
        # # i.e. in "endpoints" and "climb".
        # self.full_f = dstrip(self.forces.f)
        # self.full_v = dstrip(self.forces.pots)

        # # We need the force component perpendicular to the spline to check convergence.
        # tangent = spline_derv(self.beads.q[:, fixmask], nbeads)[1:-1]
        # f_perp = (
        #     self.full_f[1:-1, fixmask]
        #     - tangent
        #     * np.sum(self.full_f[1:-1, fixmask] * tangent, axis=1)[:, np.newaxis]
        # )

        # self.qtime += time.time()

        # # Check preliminary convergence: only the perpendicular component of the force
        # # at the NOT resampled positions - because this costs nothing.
        # if np.amax(np.abs(f_perp)) <= self.tolerances["force"]:
        #     info(
        #         " @STRING: The force before resampling has converged. Step: %i\n"
        #         % step,
        #         vrb.medium,
        #     )
        #     info(" @STRING: resampling beads and checking again.", vrb.debug)
        #     self.beads.q[:, fixmask] = spline_resample(
        #         self.beads.q[:, fixmask],
        #         nbeads,
        #     )
        #     # dx = current resampled position - previous resampled position.
        #     dx = np.amax(np.abs(self.beads.q[1:-1, fixmask] - self.old_x))

        #     # Here we check forces on the resampled positions, which costs
        #     # an additional force evaluation. We again call it via the mapper
        #     # to omit the endpoints.
        #     self.stringpot, self.stringgrad = self.stringgm(self.beads.q[1:-1, fixmask])
        #     tangent = spline_derv(self.beads.q[:, fixmask], nbeads)[1:-1]
        #     f_perp = (
        #         self.stringgrad
        #         - tangent * np.sum(self.stringgrad * tangent, axis=1)[:, np.newaxis]
        #     )
        #     if (
        #         (
        #             np.amax(np.abs(self.stringpot - old_stringpot)) / (nbeads * natoms)
        #             <= self.tolerances["energy"]
        #         )
        #         and (np.amax(np.abs(f_perp)) <= self.tolerances["force"])
        #         and (dx <= self.tolerances["position"])
        #     ):
        #         info(
        #             decor
        #             + " @STRING: path optimization converged. Step: %i\n" % step
        #             + decor,
        #             vrb.medium,
        #         )

        #         # Set climbing stage indicator after convergence of the path
        #         if self.use_climb:
        #             self.stage = "climb"
        #         else:
        #             self.stage = "converged"
        #             softexit.trigger(
        #                 status="success",
        #                 message="String MEP finished successfully at STEP %i." % step,
        #             )

        #     else:  # the 2nd check didn't pass
        #         info(
        #             " @STRING: Not converged, deltaEnergy = %.8f, tol = %.8f per atom"
        #             % (
        #                 np.amax(np.abs(self.stringpot - old_stringpot))
        #                 / (nbeads * natoms),
        #                 self.tolerances["energy"],
        #             ),
        #             vrb.debug,
        #         )
        #         info(
        #             " @STRING: Not converged, f_perp_2 = %.8f, tol = %f"
        #             % (np.amax(np.abs(self.stringgrad)), self.tolerances["force"]),
        #             vrb.debug,
        #         )
        #         info(
        #             " @STRING: Not converged, deltaX = %.8f, tol = %.8f"
        #             % (dx, self.tolerances["position"]),
        #             vrb.debug,
        #         )
        # else:  # the 1st check didn't pass
        #     info(
        #         " @STRING: Not converged, f_perp_1 = %.8f, tol = %f"
        #         % (np.amax(np.abs(self.stringgrad)), self.tolerances["force"]),
        #         vrb.debug,
        #     )

    def climb_step(self, step=None):
        """Climbing image optimization. It uses StringClimbGrMapper
        which constructs single-bead object to calculate forces
        only on the "climbing" bead.
        """
        self.ptime = self.ttime = 0
        self.qtime = -time.time()

        # Self instances will be updated in the optimizer, so we store the copies.
        # old_stringpot is used later as a convergence criterion.
        old_stringpot = self.stringpot.copy()

        if self.mode == "damped_bfgs":
            info(" @STRING_CLIMB: before Damped_BFGS() call", vrb.debug)
            print("self.old_x.shape: %s" % str(self.old_x.shape))
            print("self.stringgrad.shape: %s" % str(self.stringgrad.shape))
            print("self.hessian.shape: %s" % str(self.hessian.shape))
            Damped_BFGS(
                x0=self.old_x.copy(),
                fdf=self.climbgm,
                fdf0=(self.stringpot, self.stringgrad),
                hessian=self.hessian,
                big_step=self.big_step,
            )
            info(" @STRING_CLIMB: after Damped_BFGS() call", vrb.debug)

        elif self.mode == "bfgstrm":
            info(" @STRING_CLIMB: before BFGSTRM() call", vrb.debug)
            print("self.old_x.shape: %s" % str(self.old_x.shape))
            print("self.stringgrad.shape: %s" % str(self.stringgrad.shape))
            print("self.hessian.shape: %s" % str(self.hessian.shape))
            BFGSTRM(
                x0=self.old_x.copy(),
                u0=self.stringpot,
                # BFGSTRM expects force instead of gradient, therefore minus
                f0=-self.stringgrad,
                h0=self.hessian,
                tr=self.tr_trm,
                mapper=self.climbgm,
                big_step=self.big_step,
            )
            info(" @STRING_CLIMB: after BFGSTRM() call", vrb.debug)

        elif self.mode == "fire":
            info(" @STRING: before FIRE() call", vrb.debug)
            info(" @FIRE velocity: %s" % str(npnorm(self.v)), vrb.debug)
            info(" @FIRE alpha: %s" % str(self.a), vrb.debug)
            info(" @FIRE N down: %s" % str(self.N_dn), vrb.debug)
            info(" @FIRE N up: %s" % str(self.N_up), vrb.debug)
            info(" @FIRE dt: %s" % str(self.dt_fire), vrb.debug)
            self.v, self.a, self.N_dn, self.N_up, self.dt_fire = FIRE(
                x0=self.old_x.copy(),
                fdf=self.climbgm,
                fdf0=(self.stringpot, self.stringgrad),
                v=self.v,
                a=self.a,
                N_dn=self.N_dn,
                N_up=self.N_up,
                dt=self.dt_fire,
                dtmax=self.dtmax,
            )
            info(" @STRING_CLIMB: after FIRE() call", vrb.debug)

        # TODO: Routines for L-BFGS, SD, CG, ...
        else:
            softexit.trigger(
                status="bad",
                message="Try damped_bfgs or fire, other algorithms are not implemented for climbing image optimization.",
            )

        # Update positions
        self.beads.q[self.cl_indx] = self.climbgm.rbeads.q
        info(" @STRING_CLIMB: climbing bead position is updated.", vrb.debug)

        self.stringpot, self.stringgrad = self.climbgm(
            self.beads.q[self.cl_indx, self.climbgm.fixmask]
        )

        # Use to determine converged minimization
        # max displacement
        dx = np.amax(
            np.abs(self.beads.q[self.cl_indx, self.climbgm.fixmask] - self.old_x)
        )

        # Store old positions
        self.old_x[:] = self.beads.q[self.cl_indx, self.climbgm.fixmask]

        # This transfers forces from the climbing mapper to the "main" beads,
        # so that recalculation won't be triggered after the step.
        tmp_f = self.full_f.copy()
        tmp_v = self.full_v.copy()
        tmp_f[self.cl_indx, self.climbgm.fixmask] = self.stringgrad
        tmp_v[self.cl_indx] = self.stringpot
        self.forces.transfer_forces_manual(
            new_q=[
                self.beads.q,
            ],
            new_v=[
                tmp_v,
            ],
            new_forces=[
                tmp_f,
            ],
        )
        self.full_f = dstrip(self.forces.f)
        self.full_v = dstrip(self.forces.pots)

        self.qtime += time.time()

        # Check convergence criteria
        if (
            (
                np.amax(np.abs(self.stringpot - old_stringpot)) / self.beads.natoms
                <= self.tolerances["energy"]
            )
            and (np.amax(np.abs(self.stringgrad)) <= self.tolerances["force"])
            and (dx <= self.tolerances["position"])
        ):
            info(
                decor
                + " @STRING_CLIMB: optimization converged. Step: %i\n" % step
                + decor,
                vrb.medium,
            )
            self.stage = "converged"
            softexit.trigger(
                status="success", message="STRING_CLIMB finished successfully."
            )

        else:
            info(
                " @STRING_CLIMB: Not converged, deltaEnergy = %.8f, tol = %.8f per atom"
                % (
                    np.amax(np.abs(self.stringpot - old_stringpot)) / self.beads.natoms,
                    self.tolerances["energy"],
                ),
                vrb.debug,
            )
            info(
                " @STRING_CLIMB: Not converged, climbgrad = %.8f, tol = %f"
                % (np.amax(np.abs(self.stringgrad)), self.tolerances["force"]),
                vrb.debug,
            )
            info(
                " @STRING_CLIMB: Not converged, deltaX = %.8f, tol = %.8f"
                % (dx, self.tolerances["position"]),
                vrb.debug,
            )

    def step(self, step=None):
        """Does one simulation time step. Depending on 'stage', calls one of the
        methods to deal either with the endpoints, or the path optimization,
        or the climbing-image optimization.
        """

        info(" @STRING STEP %d, stage: %s" % (step, self.stage), vrb.debug)

        # Check if we restart a converged calculation (by mistake)
        if self.stage == "converged":
            softexit.trigger(
                status="success",
                message="String MEP has already converged. Exiting simulation.",
            )

        # First, optimization of endpoints, if required
        if self.endpoints["optimize"] and self.stage == "endpoints":
            # TODO self.endpoints_step(step)
            softexit.trigger(
                status="bad",
                message="Optimization of the endpoints of a string is not implemented yet.",
            )

        # Endpoints are optimized, or their optimization is not required
        elif self.stage == "string":
            if step == 0:  # TODO add a condition when after the endpoints.
                # Make sure the beads are equidistant before we start
                info(" @STRING: resampling beads uniformly at step 0.", vrb.debug)
                self.beads.q[:, self.stringgm.fixmask] = spline_resample(
                    self.beads.q[:, self.stringgm.fixmask],
                    self.beads.nbeads,
                )

                # Store full forces and pots
                info(" @STRING: calling the full force at step 0.", vrb.debug)
                self.full_f = dstrip(self.stringgm.dforces.f).copy()
                self.full_v = dstrip(self.stringgm.dforces.pots).copy()
                # Initialize the mapper forces for the first time
                # Since we just called full forces, it makes sense to transfer them to the
                # reduced beads (rbeads). No fixmask here, rforces are not the same
                # as stringgrad!
                print("self.stringgm.rforces.shape, self.full_f.shape:")
                print(self.stringgm.rforces.f.shape, self.full_f.shape)
                self.stringgm.rforces.transfer_forces_manual(
                    new_q=[self.beads.q[1:-1]],
                    new_v=[self.full_v[1:-1]],
                    new_forces=[self.full_f[1:-1]],
                )
                info(" @STRING: calling StringGradientMapper at step 0.", vrb.debug)
                self.stringpot, self.stringgrad = self.stringgm(
                    self.beads.q[1:-1, self.stringgm.fixmask]
                )
                info(
                    " @STRING: StringGradientMapper returned stringpot and stringgrad.",
                    vrb.debug,
                )
                # Store old bead positions at step 0
                self.old_x = dstrip(self.beads.q[1:-1, self.stringgm.fixmask]).copy()

            self.path_step_textbook2002(step)

        elif self.stage == "climb":
            # We need to initialize climbing once
            if np.all(self.climbgm.q_prev == 0.0) or np.all(self.climbgm.q_next == 0.0):
                self.cl_indx = self.init_climb()
            self.climb_step(step)


def spline_resample(q, nbeads):
    """Resamples the intermediate points along the spline so that
    all points are equidistant by the spline arc length.

    Arguments:
        q - beads.q[:]
        nbeads - number of beads
    Returns:
        new_q - resampled coordinates
    """
    if nbeads <= 2:
        softexit.trigger(status="bad", message="Nbeads < 3 in string optimization.")

    # First, we calculate the current parameterization of the path
    # according to 3N-D Euclidean distances between adjacent beads.
    t = [0.0]
    current_t = 0.0
    print("Cartesian 3N-D distances between beads:")
    for i in range(1, nbeads):
        dist = npnorm(q[i] - q[i - 1])
        print("\tfrom %3d to %3d : %.6f bohr" % (i - 1, i, dist))
        if dist <= 1e-2:
            info("Warning: two adjacent beads are very close.", vrb.medium)
        current_t += dist
        t.append(current_t)

    t = np.array(t)
    t /= t[-1]
    info("Spline parameter 't' before resampling:", vrb.medium)
    info(t, vrb.medium)
    info("t[i] - t[i-1]:", vrb.medium)
    info((t - np.roll(t, 1)), vrb.medium)

    # New parameterization, equidistant by spline parameter
    new_t = np.arange(nbeads) / (nbeads - 1)
    info("New t:", vrb.debug)
    info(new_t, vrb.debug)

    # This will reshape q to a list of 3*Natoms sets consisting of
    # nbeads numbers each, describing trajectory of one Cartesian component
    # of an atom over the path.
    atoms = q.T

    new_atoms = []
    for component in atoms:
        # Interpolate the trajectory of one cartesian component  by a spline.
        # bc_type='natural' imposes zero second derivative at the endpoints.
        spl = make_interp_spline(t, component, bc_type="natural")
        # Sample new parameter values uniformly from the old range.
        # Resample positions of this atom.
        new_atoms.append(splev(new_t, spl))

    # Reshape the path back to the beads shape
    new_q = np.array(new_atoms).T

    return new_q


def spline_derv(q, nbeads):
    """Calculates normalized tangent vectors at the knots of the spline.

    Arguments:
        q - beads.q[:], shape (nbeads, 3N)
        nbeads - number of beads
    Returns:
        tangent - full-dimensional tangent of the spline, shape (nbeads, 3N).
    """
    if nbeads <= 2:
        softexit.trigger(status="bad", message="Nbeads < 3 in string optimization.")

    # First, we calculate current parameterization of the path
    # according to 3N-D Euclidean distances between adjacent beads.
    t = [0.0]
    current_t = 0.0
    for i in range(1, nbeads):
        dist = npnorm(q[i] - q[i - 1])
        if dist <= 1e-2:
            info("Warning: two adjacent beads are very close.", vrb.medium)
            print("dist[%i] = %.6f" % (i, dist))
        current_t += dist
        t.append(current_t)

    t = np.array(t)
    t /= t[-1]
    # This will reshape q to a list of 3*Natoms sets consisting of
    # nbeads numbers each, describing trajectory of one Cartesian component
    # of an atom over the path.
    atoms = q.T

    derv = []
    for component in atoms:
        # Interpolate the trajectory of one atom by a spline.
        # s=0 imposes exact interpolation at the knots.
        spl = make_interp_spline(t, component, bc_type="natural")
        # Evaluate derivatives
        derv.append(splev(t, spl, der=1))

    # Reshape the array back to the beads shape
    tangent = np.array(derv).T
    # Normalize so that at each bead it's a unitary vector
    tangent /= npnorm(tangent, axis=1)[:, np.newaxis]
    # print("tangent:")
    # print(tangent.reshape((-1, 3)))
    return tangent


def print_distances(q, nbeads):
    """Prints distances by the 3N-dimensional Euclidean metric.

    Arguments:
        q - beads.q[:]
        nbeads - number of beads
    """
    print("Cartesian 3N-D distances between beads:")
    for i in range(1, nbeads):
        dist = npnorm(q[i] - q[i - 1])
        print("\tfrom %3d to %3d : %.6f bohr" % (i - 1, i, dist))
        if dist <= 1e-2:
            info("Warning: two adjacent beads are very close.", vrb.medium)


# Damped BFGS to use in String ONLY. It resamples the spline
# right after the quasi-Newton step.
# Has no line search and no TRM.
# TODO: It's a placeholder to start from.
# def Damped_BFGS_String(x0, fdf, fdf0, hessian, big_step):
#    """BFGS, damped as described in Nocedal, Wright (2nd ed.) Procedure 18.2
#    The purpose is SOLELY using it for String optimization, because it reparameterizes
#    the string right after quasi-Newton step. The reason for such complication is
#    to avoid calculating forces twice - on resampled and non-resampled positions.
#
#    Written for a DIRECT Hessian B, not for the inverse H.
#
#    I use flattened vectors inside this function, restoring the shape only when needed.
#    I always keep x0 in the original shape.
#
#    Does one step.
#
#    Arguments:
#      x0: initial point
#      fdf: function and gradient (mapper)
#      fdf0: initial function and gradient values
#      hessian: approximate Hessian for the BFGS algorithm
#      big_step: limit on step length. It is defined differently
#                  compared to other optimization algorithms, take care.
#
#    Returns:
#      quality: minus cosine of the (gradient, dx) angle.
#               Needed for the step length adjustment.
#    """
#
#    info(" @DampedBFGS_String: Started.", vrb.debug)
#    _, g0 = fdf0
#    g0 = g0.flatten()
#
#    # Nocedal's notation
#    B = hessian
#
#    # Calculate direction
#    # When the inverse itself is not needed, people recommend solve(), not inv().
#    info(" @DampedBFGS_String: sk = np.linalg.solve(B, -g0) ...", vrb.debug)
#    info(
#        "              The code randomly crashes here with some versions of Numpy "
#        "based on OpenBLAS.\n"
#        "              If this happens, use Numpy based on MKL, e.g. from Anaconda.",
#        vrb.debug,
#    )
#    info("Operands:", vrb.debug)
#    info("%s,  %s" % (type(B), str(B.shape)), vrb.debug)
#    info("%s,  %s" % (type(g0), str(g0.shape)), vrb.debug)
#    sk = np.linalg.solve(B, -g0)
#    info(" @DampedBFGS_String: Calculated direction.", vrb.debug)
#
#    # Cosine of the (f, dx) angle
#    quality = -np.dot(sk / np.linalg.norm(sk), g0 / np.linalg.norm(g0))
#    info(" @DampedBFGS_String: Direction quality: %.4f." % quality, vrb.debug)
#
#    # I use maximal cartesian atomic displacement as a measure of step length
#    maxdispl = np.amax(np.linalg.norm(sk.reshape(-1, 3), axis=1))
#    info(" @DampedBFGS_String: big_step = %.6f" % big_step, vrb.debug)
#    if maxdispl > big_step:
#        info(
#            " @DampedBFGS_String: maxdispl before scaling: %.6f bohr" % maxdispl,
#            vrb.debug,
#        )
#        sk *= big_step / maxdispl
#
#    info(
#        " @DampedBFGS_String: maxdispl:                %.6f bohr"
#        % (np.amax(np.linalg.norm(sk.reshape(-1, 3), axis=1))),
#        vrb.debug,
#    )
#
#    # Force call
#    _, g = fdf(x0 + sk.reshape(x0.shape))
#    g = g.flatten()
#    # coordinates CHECKED
#
#    # Update hessian
#    yk = g - g0
#    skyk = np.dot(sk, yk)
#
#    # Equation 18.15 in Nocedal
#    theta = 1.0
#    Bsk = np.dot(B, sk)
#    sBs = np.dot(sk, Bsk)
#    # Damped update if rhok isn't sufficiently positive
#    if skyk < 0.2 * sBs:
#        theta = (0.8 * sBs) / (sBs - skyk)
#        info(
#            " @DampedBFGS_String: damping update of the Hessian; "
#            "(direction dot d_gradient) is small. "
#            "theta := %.6f" % theta,
#            vrb.debug,
#        )
#        yk = theta * yk + (1 - theta) * Bsk
#        skyk = np.dot(sk, yk)
#    else:
#        info(" @DampedBFGS_String: Update of the Hessian, no damping.", vrb.debug)
#
#    info(
#        " @DampedBFGS_String: (s_k dot y_k) before reciprocating: %e" % skyk,
#        vrb.debug,
#    )
#    try:
#        rhok = 1.0 / skyk
#    except:
#        warning(" @DampedBFGS_String: caught ZeroDivisionError in 1/skyk.", vrb.high)
#        rhok = 1e5
#
#    # Compute BFGS term (eq. 18.16 in Nocedal)
#    B += np.outer(yk, yk) * rhok - np.outer(Bsk, Bsk) / sBs
#
#    print("det(Hessian): %f" % np.linalg.det(B))
#    print("Condition number: %f" % np.linalg.cond(B))
#    eigvals = np.real(np.linalg.eigvals(B))
#    print("Positive definite: %r" % np.all(eigvals > 0))
#    print("Eigenvalues of the Hessian:")
#    with np.printoptions(threshold=10000, linewidth=100000):
#        print(np.sort(eigvals))
#
#    # If small numbers are found on the diagonal of the Hessian,
#    # add small positive numbers. 1 Ha/Bohr^2 is ~97.2 eV/ang^2.
#    if np.any(eigvals < 1e-1):
#        info(
#            " @DampedBFGS_String: stabilizing the diagonal of the Hessian.",
#            vrb.debug,
#        )
#        B += 1e-2 * np.eye(len(B))
#
#    return quality


# FIRE to use in String ONLY. It resamples the spline
# right after FIRE displaces atoms.
# TODO: It's a placeholder to start from.
# def FIRE_String(
#    x0,
#    fdf,
#    fdf0,
#    v=None,
#    a=0.1,
#    N_dn=0,
#    N_up=0,
#    dt=0.1,
#    maxstep=0.5,
#    dtmax=1.0,
#    dtmin=1e-5,
#    Ndelay=5,
#    Nmax=2000,
#    finc=1.1,
#    fdec=0.5,
#    astart=0.1,
#    fa=0.99,
# ):
#    """FIRE algorithm based on
#    Bitzek et al, Phys. Rev. Lett. 97, 170201 (2006) and
#    Guénolé, J. et al.  Comp. Mat. Sci. 175, 109584 (2020).
#    Semi-implicit Euler integration used.
#    Done by Guoyuan Liu <liuthepro@outlook.com>, May 2021.
#
#    FIRE does not rely on energy, therefore it is suitable for NEB calculation, where
#    the energy is not conservative. Basic principle: accelerate towards force gradient
#    (downhill direction) and stop immediately when going uphill.
#    Try adjusting dt, dtmax, dtmin for optimal performance.
#
#    Arguments:
#        x0: initial beads positions
#        fdf: energy and function mapper. call fdf(x) to update beads position and froces
#        fdf0: initial value of energy and gradient
#        v: current velocity
#        a: velocity mixing factor, in the paper it is called alpha
#        fa: a decrement factor
#        astart: initial a value
#        N_dn: number of steps since last downhill direction
#        N_up: number of steps since last uphill direction
#        dt: time interval
#        dtmax: max dt (increase when uphill)
#        dtmin: min dt (decrease when downhill)
#        finc: dt increment factor
#        fdec: dt decrement factor
#        Ndelay: min steps required to be in one direction before adjust dt and a
#        Nmax: max consecutive steps in uphill direction before trigger exit
#
#    Returns:
#        v, a, N, dt since they are dynamically adjusted
#    """
#    info(" @FIRE_String being called", vrb.debug)
#    _, g0 = fdf0
#    force = -g0
#
#    p = np.vdot(force, v)
#    # downhill
#    if p > 0.0:
#        N_dn += 1
#        N_up = 0
#        if N_dn > Ndelay:
#            dt = min(dt * finc, dtmax)
#            a = a * fa
#    # uphill
#    else:
#        N_dn = 0
#        N_up += 1
#        if N_up > Nmax:
#            softexit.trigger("@FIRE is stuck for %d steps. We stop here." % N_up)
#        dt = max(dt * fdec, dtmin)
#        a = astart
#        # correct uphill motion
#        x0 -= 0.5 * dt * v
#        # stop moving in uphill direction
#        v = np.zeros(v.shape)
#
#    # accelerate
#    v += dt * force
#    # change velocity direction with inertia
#    if p > 0.0:
#        f_unit = force / np.linalg.norm(force)
#        v = (1 - a) * v + a * np.linalg.norm(v) * f_unit
#    # update posistion
#    dx = dt * v
#    # check max dx
#    normdx = np.linalg.norm(dx)
#    if normdx > maxstep:
#        dx = maxstep * dx / normdx
#    x0 += dx
#
#    info(" @FIRE_String calling a gradient mapper to update position", vrb.debug)
#    fdf(x0)
#
#    return v, a, N_dn, N_up, dt

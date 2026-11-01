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
from ipi.utils.array_backend import xp, xp_size, to_numpy
from ipi.utils.mintools import Damped_BFGS, BFGSTRM, FIRE
from ipi.utils.messages import verbosity as vrb, info  # , warning
from ipi.engine.beads import Beads

try:
    import scipy
    from scipy.interpolate import make_interp_spline, splev
except Exception as e:
    scipy = None
    scipy_exception = e

decor = 60 * "=" + "\n"

__all__ = [
    "StringGradientMapper",
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
        self.dbeads = ens.beads.clone()
        self.dcell = ens.cell.clone()
        self.dforces = ens.forces.clone(self.dbeads, self.dcell)
        self.fixatoms_dof = ens.fixatoms_dof.copy()

        # Mask to exclude fixed atoms from 3N-arrays
        self.fixatoms_mask = np.ones(3 * ens.beads.natoms, dtype=bool)
        if len(ens.fixatoms_dof) > 0:
            self.fixatoms_mask[self.fixatoms_dof] = 0

        # Create reduced bead and force object
        self.rbeads = Beads(ens.beads.natoms, ens.beads.nbeads - 2)
        self.rbeads.q_[:] = ens.beads.q[1:-1]
        self.rforces = ens.forces.clone(self.rbeads, self.dcell)

    def __call__(self, x):
        """Returns the potential for all beads and the gradient."""

        # Bead positions
        # Touch positions only if they have changed (to avoid triggering forces)
        # I need both dbeads and rbeads because of the endpoint tangents.
        if (self.rbeads.q[:, self.fixatoms_mask] != x).any():
            self.rbeads.q_[:, self.fixatoms_mask] = x
        rbq = self.rbeads.q[:, self.fixatoms_mask]
        # We need all beads, but only reduced number is being updated.
        self.dbeads.q_[1:-1, self.fixatoms_mask] = rbq

        # Forces and energies (boolean-mask indexing returns a copy)
        rbf = self.rforces.f[:, self.fixatoms_mask]
        be = xp.asarray(self.rforces.pots, copy=True)

        # Calculate the force component perpendicular to the spline.
        tangent = spline_derv(self.dbeads.q[:, self.fixatoms_mask], self.dbeads.nbeads)[
            1:-1
        ]
        rbf_perp = rbf - tangent * (rbf * tangent).sum(axis=1)[:, None]

        # Return potential energy of the whole string and minus f_perpendicular
        e = be.sum()
        g = -rbf_perp
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
        self.fixatoms_dof = ens.fixatoms_dof.copy()
        # A mask to exclude fixed atoms from 3N-arrays
        self.fixatoms_mask = np.ones(3 * ens.beads.natoms, dtype=bool)
        if len(ens.fixatoms_dof) > 0:
            self.fixatoms_mask[ens.fixatoms_dof] = 0

        # Reduced Beads object is needed to calculate only required beads.
        self.rbeads = Beads(ens.beads.natoms, 1)
        self.rcell = ens.cell.clone()
        # Coords of the bead before and after the climbing one.
        self.q_prev = xp.zeros(3 * ens.beads.natoms - len(ens.fixatoms_dof))
        self.q_next = xp.zeros(3 * ens.beads.natoms - len(ens.fixatoms_dof))
        # Make reduced forces dependent on reduced beads
        self.rforces = ens.forces.clone(self.rbeads, self.rcell)

    def __call__(self, x):
        """Returns climbing force for a climbing image."""
        if self.q_prev is None or self.q_next is None:
            raise RuntimeError("StringClimbGrMapper.q_prev or q_next is None.")

        info(
            "@STRING_CLIMB: StringClimbGrMapper.x.shape: %s" % str(x.shape),
            vrb.debug,
        )
        # Touch positions only if they have changed (to avoid triggering forces)
        if (self.rbeads.q[:, self.fixatoms_mask] != x).any():
            self.rbeads.q_[:, self.fixatoms_mask] = x
        # Compared to stringgradientMapper, I need to flatten here,
        # otherwise it's (1, 3N) instead of just 3N.
        rbq = self.rbeads.q[:, self.fixatoms_mask].reshape(-1)

        # Reduced forces
        # We don't need copy() because boolean-mask indexing implies copying.
        rbf = self.rforces.f[:, self.fixatoms_mask].reshape(-1)

        # I think here it's better to use plain tangents.
        # Then we don't need energies of the neighboring beads.
        norm = xp.linalg.vector_norm
        d1 = rbq - self.q_prev  # tau minus
        d2 = self.q_next - rbq  # tau plus
        tau = d1 / norm(d1) + d2 / norm(d2)
        tau /= norm(tau)

        # Inverting the force component along the path
        rbf -= 2 * (rbf @ tau) * tau

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

    jax_support = False

    def __init__(
        self,
        fixcom=False,
        fixatoms_dof=None,
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

        super(StringMover, self).__init__(fixcom=fixcom, fixatoms_dof=fixatoms_dof)

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
        self.old_x = xp.asarray(old_coord)  # always in "masked" dimension
        self.stringpot = xp.asarray(old_stringpotential)
        self.stringgrad = xp.asarray(old_stringgradient)
        self.d = xp.asarray(old_direction)
        # physical forces of ALL beads even in climb stage
        self.full_f = xp.asarray(full_force)
        self.full_v = xp.asarray(full_pots)  # potentials of ALL beads
        self.corrections = corrections_lbfgs
        self.qlist = xp.asarray(qlist_lbfgs)
        self.glist = xp.asarray(glist_lbfgs)
        self.scale = scale_lbfgs
        self.tr_trm = tr_trm  # trust radius, kept on the host
        self.endpoints = endpoints
        self.use_climb = use_climb
        self.cl_indx = climb_bead
        self.stage = stage
        # damped_bfgs
        self.hessian = xp.asarray(hessian_bfgs)
        # fire
        self.v = xp.asarray(v_fire)  # velocity
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

        if len(self.fixatoms_dof) == 3 * len(self.beads[0]):
            raise ValueError(
                "All atoms are fixed, geometry won't change. Exiting simulation."
            )

        # the mappers hold the mask of the moving coordinates, used below
        self.stringgm.bind(self)
        self.climbgm.bind(self)

        # We check dimensionality of the Hessian in the beginning
        # and reduce it if needed, because someone may want to provide
        # existing Hessian of the full system.
        if self.mode in ["damped_bfgs", "bfgstrm"]:
            n_activedim = xp_size(beads.q[0]) - len(self.fixatoms_dof)
            if self.stage == "string":
                if xp_size(self.hessian) == (n_activedim * (beads.nbeads - 2)) ** 2:
                    # Desired dimension
                    info(
                        " @STRINGMover: String Hessian is treated as reduced "
                        "and masked already, according to its size.",
                        vrb.high,
                    )
                elif xp_size(self.hessian) == xp_size(beads.q) ** 2:
                    info(
                        " @STRINGMover: String Hessian has full-dimensional size, "
                        "will be reduced and masked by fixatoms.",
                        vrb.high,
                    )
                    # First trim the endpoints, then mask fixed atoms
                    self.hessian = (
                        self.hessian[
                            3 * self.beads.natoms : -3 * self.beads.natoms,
                            3 * self.beads.natoms : -3 * self.beads.natoms,
                        ]
                    )[
                        np.ix_(
                            np.tile(self.stringgm.fixatoms_mask, self.beads.nbeads - 2),
                            np.tile(self.stringgm.fixatoms_mask, self.beads.nbeads - 2),
                        )
                    ]
                elif xp_size(self.hessian) == 0:
                    info(
                        " @STRINGMover: No Hessian provided, starting from unity matrix",
                        vrb.high,
                    )
                    self.hessian = xp.eye((self.beads.nbeads - 2) * n_activedim)
                else:
                    raise ValueError("Hessian size does not match system size.")
            # In climb stage, hessian has only 1 bead.
            # Not sure whether fixatoms should be handled here or not, let's see...
            elif self.stage == "climb":
                if xp_size(self.hessian) == n_activedim * n_activedim:
                    # desired dimension
                    info(
                        " @STRINGMover: Hessian of the climber is treated as masked one,"
                        " according to its size.",
                        vrb.high,
                    )
                elif xp_size(self.hessian) == xp_size(beads.q[0]) ** 2:
                    info(
                        " @STRINGMover: Hessian of the climber has size (3natoms x 3natoms), "
                        "will be masked by fixatoms.",
                        vrb.high,
                    )
                    self.hessian = self.hessian[
                        np.ix_(self.climbgm.fixatoms_mask, self.climbgm.fixatoms_mask)
                    ]
                elif xp_size(self.hessian) == 0:
                    info(
                        " @STRINGMover: No Hessian provided for climbing, "
                        "starting from unity matrix",
                        vrb.high,
                    )
                    self.hessian = xp.eye(n_activedim)
                else:
                    raise ValueError("Hessian size does not match system size.")

    def init_climb(self):
        """Operations to resize hessian for the climbing bead(s),
        get neighboring beads and other things necessary for climbing.
        Returns:
          cl_indx - index of the bead for climbing
        """
        cl_indx = int(xp.argmax(self.forces.pots))
        info(
            " @STRING: Initializing climbing. Climbing bead: %i." % cl_indx,
            vrb.medium,
        )
        if cl_indx in [0, self.beads.nbeads - 1]:
            self.finish(status="bad", message="ERROR: climbing bead is the endpoint.")
            return None
        self.climbgm.rbeads.q_[:] = self.beads.q[cl_indx]
        self.climbgm.q_prev[:] = self.beads.q[cl_indx - 1, self.climbgm.fixatoms_mask]
        self.climbgm.q_next[:] = self.beads.q[cl_indx + 1, self.climbgm.fixatoms_mask]
        info(
            "q_prev.shape: %s,    q_next.shape: %s"
            % (str(self.climbgm.q_prev.shape), str(self.climbgm.q_next.shape)),
            vrb.debug,
        )

        info(" @STRING_CLIMB: call StringClimbGrMapper for the first time.", vrb.debug)
        self.stringpot, self.stringgrad = self.climbgm(
            self.beads.q[cl_indx, self.climbgm.fixatoms_mask]
        )

        if self.mode in ["damped_bfgs", "bfgstrm"]:
            # Initialize BFGS Hessian for a single bead
            n_activedim = xp_size(self.beads.q[0]) - len(self.fixatoms_dof)
            if self.hessian.shape != (n_activedim, n_activedim):
                self.hessian = xp.eye(n_activedim)
        elif self.mode == "fire":
            # Initialize FIRE parameters
            self.a = 0.1
            self.v = -self.a * self.stringgrad
            self.N_dn = 0
            self.N_up = 0
            self.dt_fire = 0.1

        self.old_x = xp.asarray(
            self.beads.q[cl_indx, self.climbgm.fixatoms_mask], copy=True
        )
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
        ones, therefore it's necessary to have 2 force calls per step.
        """
        self.ptime = self.ttime = 0
        self.qtime = -time.time()
        # Shortcuts for prettier expresions
        nbeads = self.beads.nbeads
        natoms = self.beads.natoms
        fixatoms_mask = self.stringgm.fixatoms_mask

        n_activedim = xp_size(self.beads.q[0]) - len(self.fixatoms_dof)

        # Store 'old' resampled positions
        self.old_x[:] = self.beads.q[1:-1, fixatoms_mask]

        # 'Self' stringpot will be updated later on,
        # therefore we store the copy to check convergence in the end.
        old_stringpot = xp.asarray(self.stringpot, copy=True)

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
                    raise ValueError(
                        f"Wrong Hessian size in String step: {self.hessian.shape} "
                        f"for beads of shape {self.beads.q.shape}."
                    )

            if self.mode == "damped_bfgs":
                info(" @STRING: before Damped_BFGS() call", vrb.debug)
                info("self.old_x.shape: %s" % str(self.old_x.shape), vrb.debug)
                info(
                    "self.stringgrad.shape: %s" % str(self.stringgrad.shape), vrb.debug
                )
                info("self.hessian.shape: %s" % str(self.hessian.shape), vrb.debug)
                Damped_BFGS(
                    x0=xp.asarray(self.old_x, copy=True),
                    fdf=self.stringgm,
                    fdf0=(self.stringpot, self.stringgrad),
                    hessian=self.hessian,
                    big_step=self.big_step,
                    mask=fixatoms_mask,
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
                    x0=xp.asarray(self.old_x, copy=True),
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
            if step == 0 and xp_size(self.v) == 0:
                # Mapper has already been called in the step()
                self.v = -self.a * self.stringgrad

            info(" @STRING: calling FIRE", vrb.debug)
            info(
                " @FIRE velocity: %s" % str(float(xp.linalg.vector_norm(self.v))),
                vrb.debug,
            )
            info(" @FIRE alpha: %s" % str(self.a), vrb.debug)
            info(" @FIRE N down: %s" % str(self.N_dn), vrb.debug)
            info(" @FIRE N up: %s" % str(self.N_up), vrb.debug)
            info(" @FIRE dt: %s" % str(self.dt_fire), vrb.debug)
            self.v, self.a, self.N_dn, self.N_up, self.dt_fire = FIRE(
                x0=xp.asarray(self.old_x, copy=True),
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
            if xp.any(xp.abs(dq) > self.big_step):
                info(" @EULER: WARNING: hit the big_step.")
                dq *= self.big_step / xp.max(xp.abs(dq))
            info(" @EULER: maxdispl is %g" % xp.max(xp.abs(dq)))
            self.beads.q_[1:-1, fixatoms_mask] += dq.reshape((nbeads - 2, -1))
            # info(" @EULER: calling the mapper.", vrb.debug)
            self.stringgm.dbeads.q_[:] = self.beads.q[:]
            # self.stringpot, self.stringgrad = self.stringgm(self.beads.q[1:-1, fixmask])

        # TODO: Routines for L-BFGS, SD, CG
        else:
            raise NotImplementedError(
                f"Mode '{self.mode}' is not implemented for string optimization: "
                "try 'damped_bfgs', 'bfgstrm' or 'fire'."
            )

        # Update 'main' beads positions
        self.beads.q_[:] = self.stringgm.dbeads.q
        info(" @STRING: bead positions transferred from mapper.", vrb.debug)

        # Resample the beads after an optimizer step.
        info(" @STRING: resampling beads uniformly.", vrb.debug)
        self.beads.q_[:, fixatoms_mask] = spline_resample(
            self.beads.q[:, fixatoms_mask],
            nbeads,
        )

        # Here the forces are called on resampled positions.
        info(" @STRING: calling forces on resampled positions.", vrb.debug)
        self.stringpot, self.stringgrad = self.stringgm(
            self.beads.q[1:-1, fixatoms_mask]
        )
        # This transfers forces from the mapper beads to the "main" beads,
        # so that recalculation won't be triggered after the step.
        # We need to keep forces up to date to be able to output
        # full potentials and full forces.
        tmp_v = xp.asarray(self.full_v, copy=True)
        tmp_v[1:-1] = self.stringgm.rforces.pots
        tmp_f = xp.asarray(self.full_f, copy=True)
        tmp_f[1:-1] = self.stringgm.rforces.f
        # Forcefield results are stored as numpy arrays
        self.forces.transfer_forces_manual(
            new_q=[self.beads.q],
            new_v=[to_numpy(tmp_v)],
            new_forces=[to_numpy(tmp_f)],
        )
        # Full-dimensional forces are needed in reduced-beads modes,
        # i.e. in "endpoints" and "climb".
        self.full_f = self.forces.f
        self.full_v = self.forces.pots

        # Calculate the force component perpendicular to the spline to check convergence.
        tangent = spline_derv(self.beads.q[:, fixatoms_mask], nbeads)[1:-1]
        f_perp = (
            self.full_f[1:-1, fixatoms_mask]
            - tangent
            * (self.full_f[1:-1, fixatoms_mask] * tangent).sum(axis=1)[:, None]
        )

        # Check convergence at the resampled positions.
        # dx = current resampled position - previous resampled position.
        dx = xp.max(xp.abs(self.beads.q[1:-1, fixatoms_mask] - self.old_x))
        de = xp.max(xp.abs(self.stringpot - old_stringpot)) / (nbeads * natoms)
        if (
            (de <= self.tolerances["energy"])
            and (dx <= self.tolerances["position"])
            and (xp.max(xp.abs(f_perp)) <= self.tolerances["force"])
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
                self.finish(
                    status="success",
                    message="String MEP finished successfully at STEP %i." % step,
                )
                return

        else:
            info(
                " @STRING: Not converged, deltaEnergy = %.8f, tol = %.8f per atom"
                % (de.item(), self.tolerances["energy"]),
                vrb.debug,
            )
            info(
                " @STRING: Not converged, f_perp_2 = %.8f, tol = %f"
                % (xp.max(xp.abs(f_perp)).item(), self.tolerances["force"]),
                vrb.debug,
            )
            info(
                " @STRING: Not converged, deltaX = %.8f, tol = %.8f"
                % (dx.item(), self.tolerances["position"]),
                vrb.debug,
            )

        self.qtime += time.time()

    def climb_step(self, step=None):
        """Climbing image optimization. It uses StringClimbGrMapper
        which constructs single-bead object to calculate forces
        only on the "climbing" bead.
        """
        self.ptime = self.ttime = 0
        self.qtime = -time.time()

        # Self instances will be updated in the optimizer, so we store the copies.
        # old_stringpot is used later as a convergence criterion.
        old_stringpot = xp.asarray(self.stringpot, copy=True)

        if self.mode == "damped_bfgs":
            info(" @STRING_CLIMB: before Damped_BFGS() call", vrb.debug)
            info("self.old_x.shape: %s" % str(self.old_x.shape), vrb.debug)
            info("self.stringgrad.shape: %s" % str(self.stringgrad.shape), vrb.debug)
            info("self.hessian.shape: %s" % str(self.hessian.shape), vrb.debug)
            Damped_BFGS(
                x0=xp.asarray(self.old_x, copy=True),
                fdf=self.climbgm,
                fdf0=(self.stringpot, self.stringgrad),
                hessian=self.hessian,
                big_step=self.big_step,
                mask=self.climbgm.fixatoms_mask,
            )
            info(" @STRING_CLIMB: after Damped_BFGS() call", vrb.debug)

        elif self.mode == "bfgstrm":
            info(" @STRING_CLIMB: before BFGSTRM() call", vrb.debug)
            info("self.old_x.shape: %s" % str(self.old_x.shape), vrb.debug)
            info("self.stringgrad.shape: %s" % str(self.stringgrad.shape), vrb.debug)
            info("self.hessian.shape: %s" % str(self.hessian.shape), vrb.debug)
            BFGSTRM(
                x0=xp.asarray(self.old_x, copy=True),
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
            info(
                " @FIRE velocity: %s" % str(float(xp.linalg.vector_norm(self.v))),
                vrb.debug,
            )
            info(" @FIRE alpha: %s" % str(self.a), vrb.debug)
            info(" @FIRE N down: %s" % str(self.N_dn), vrb.debug)
            info(" @FIRE N up: %s" % str(self.N_up), vrb.debug)
            info(" @FIRE dt: %s" % str(self.dt_fire), vrb.debug)
            self.v, self.a, self.N_dn, self.N_up, self.dt_fire = FIRE(
                x0=xp.asarray(self.old_x, copy=True),
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
            self.finish(
                status="bad",
                message="Try damped_bfgs or fire, other algorithms are not implemented for climbing image optimization.",
            )
            return

        # Update positions
        self.beads.q_[self.cl_indx] = self.climbgm.rbeads.q
        info(" @STRING_CLIMB: climbing bead position is updated.", vrb.debug)

        self.stringpot, self.stringgrad = self.climbgm(
            self.beads.q[self.cl_indx, self.climbgm.fixatoms_mask]
        )

        # Use to determine converged minimization
        # max displacement
        dx = xp.max(
            xp.abs(self.beads.q[self.cl_indx, self.climbgm.fixatoms_mask] - self.old_x)
        )

        # Store old positions
        self.old_x[:] = self.beads.q[self.cl_indx, self.climbgm.fixatoms_mask]

        # This transfers forces from the climbing mapper to the "main" beads,
        # so that recalculation won't be triggered after the step.
        tmp_f = xp.asarray(self.full_f, copy=True)
        tmp_v = xp.asarray(self.full_v, copy=True)
        # the physical forces of the climbing bead, not the climbing ones
        tmp_f[self.cl_indx] = self.climbgm.rforces.f[0]
        tmp_v[self.cl_indx] = self.stringpot
        # Forcefield results are stored as numpy arrays
        self.forces.transfer_forces_manual(
            new_q=[
                self.beads.q,
            ],
            new_v=[
                to_numpy(tmp_v),
            ],
            new_forces=[
                to_numpy(tmp_f),
            ],
        )
        self.full_f = self.forces.f
        self.full_v = self.forces.pots

        self.qtime += time.time()

        # Check convergence criteria
        if (
            (
                xp.max(xp.abs(self.stringpot - old_stringpot)) / self.beads.natoms
                <= self.tolerances["energy"]
            )
            and (xp.max(xp.abs(self.stringgrad)) <= self.tolerances["force"])
            and (dx <= self.tolerances["position"])
        ):
            info(
                decor
                + " @STRING_CLIMB: optimization converged. Step: %i\n" % step
                + decor,
                vrb.medium,
            )
            self.stage = "converged"
            self.finish(status="success", message="STRING_CLIMB finished successfully.")
            return

        else:
            info(
                " @STRING_CLIMB: Not converged, deltaEnergy = %.8f, tol = %.8f per atom"
                % (
                    (
                        xp.max(xp.abs(self.stringpot - old_stringpot))
                        / self.beads.natoms
                    ).item(),
                    self.tolerances["energy"],
                ),
                vrb.debug,
            )
            info(
                " @STRING_CLIMB: Not converged, climbgrad = %.8f, tol = %f"
                % (
                    xp.max(xp.abs(self.stringgrad)).item(),
                    self.tolerances["force"],
                ),
                vrb.debug,
            )
            info(
                " @STRING_CLIMB: Not converged, deltaX = %.8f, tol = %.8f"
                % (dx.item(), self.tolerances["position"]),
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
            self.finish(
                status="success",
                message="String MEP has already converged. Exiting simulation.",
            )
            return

        # First, optimization of endpoints, if required
        if self.endpoints["optimize"] and self.stage == "endpoints":
            # TODO self.endpoints_step(step)
            raise NotImplementedError(
                "Optimization of the endpoints of a string is not implemented yet."
            )

        # Endpoints are optimized, or their optimization is not required
        elif self.stage == "string":
            if step == 0:  # TODO add a condition when after the endpoints.
                # Make sure the beads are equidistant before we start
                info(" @STRING: resampling beads uniformly at step 0.", vrb.debug)
                self.beads.q_[:, self.stringgm.fixatoms_mask] = spline_resample(
                    self.beads.q[:, self.stringgm.fixatoms_mask],
                    self.beads.nbeads,
                )

                # Store full forces and pots
                info(" @STRING: calling the full force at step 0.", vrb.debug)
                self.full_f = xp.asarray(self.stringgm.dforces.f, copy=True)
                self.full_v = xp.asarray(self.stringgm.dforces.pots, copy=True)
                # Initialize the mapper forces for the first time
                # Since we just called full forces, it makes sense to transfer them to the
                # reduced beads (rbeads). No fixmask here, rforces are not the same
                # as stringgrad! Forcefield results are stored as numpy arrays.
                info(
                    " @STRING: shapes of rforces.f and full_f: %s, %s"
                    % (self.stringgm.rforces.f.shape, self.full_f.shape),
                    vrb.debug,
                )
                self.stringgm.rforces.transfer_forces_manual(
                    new_q=[self.beads.q[1:-1]],
                    new_v=[to_numpy(self.full_v[1:-1])],
                    new_forces=[to_numpy(self.full_f[1:-1])],
                )
                info(" @STRING: calling StringGradientMapper at step 0.", vrb.debug)
                self.stringpot, self.stringgrad = self.stringgm(
                    self.beads.q[1:-1, self.stringgm.fixatoms_mask]
                )
                info(
                    " @STRING: StringGradientMapper returned stringpot and stringgrad.",
                    vrb.debug,
                )
                # Store old bead positions at step 0
                self.old_x = xp.asarray(
                    self.beads.q[1:-1, self.stringgm.fixatoms_mask], copy=True
                )

            self.path_step_textbook2002(step)

        elif self.stage == "climb":
            # We need to initialize climbing once
            if xp.all(self.climbgm.q_prev == 0.0) or xp.all(self.climbgm.q_next == 0.0):
                self.cl_indx = self.init_climb()
                if self.finished:
                    return
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
        raise ValueError("Nbeads < 3 in string optimization.")

    # The spline is computed by scipy, on the host
    q = to_numpy(q)

    # First, we calculate the current parameterization of the path
    # according to 3N-D Euclidean distances between adjacent beads.
    t = [0.0]
    current_t = 0.0
    for i in range(1, nbeads):
        dist = npnorm(q[i] - q[i - 1])
        info(
            " @STRING: distance from %3d to %3d: %.6f bohr" % (i - 1, i, dist),
            vrb.debug,
        )
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

    return xp.asarray(new_q)


def spline_derv(q, nbeads):
    """Calculates normalized tangent vectors at the knots of the spline.

    Arguments:
        q - beads.q[:], shape (nbeads, 3N)
        nbeads - number of beads
    Returns:
        tangent - full-dimensional tangent of the spline, shape (nbeads, 3N).
    """
    if nbeads <= 2:
        raise ValueError("Nbeads < 3 in string optimization.")

    # The spline is computed by scipy, on the host
    q = to_numpy(q)

    # First, we calculate current parameterization of the path
    # according to 3N-D Euclidean distances between adjacent beads.
    t = [0.0]
    current_t = 0.0
    for i in range(1, nbeads):
        dist = npnorm(q[i] - q[i - 1])
        if dist <= 1e-2:
            info(
                "Warning: beads %d and %d are very close (%.6f bohr)."
                % (i - 1, i, dist),
                vrb.medium,
            )
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
    return xp.asarray(tangent)

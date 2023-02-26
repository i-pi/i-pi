"""Holds the algorithms to perform nudged elastic band (NEB) calculations.
J. Chem. Phys. 113, 9901 (2000); https://doi.org/10.1063/1.1329672

The algorithms are first implemented by Michele Ceriotti and Benjamin Helfrecht, 2015.
Considerably reworked by Karen Fidanyan in 2021.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2021 i-PI developers
# See the "licenses" directory for full license information.

import numpy as np
from numpy.linalg import norm as npnorm
import time

from ipi.engine.motion import Motion
from ipi.utils.depend import dstrip
from ipi.utils.softexit import softexit
from ipi.utils.mintools import Damped_BFGS, FIRE
from ipi.utils.messages import verbosity, info
from ipi.engine.beads import Beads

np.set_printoptions(threshold=10000, linewidth=1000)  # Remove in cleanup

__all__ = ["NEBGradientMapper", "NEBClimbGrMapper", "NEBMover"]


class NEBGradientMapper(object):
    """Creation of the multi-dimensional function that will be minimized.
        Functional analog of a GradientMapper in geop.py

        Fixed atoms are excluded via boolean mask. 1 = moving, 0 = fixed.

    Attributes:
        kappa: spring constants
        tangent: plain or improved tangents
          "plain":    J. Chem. Phys. 113, 9901 (2000); https://doi.org/10.1063/1.1329672
          "improved": J. Chem. Phys. 113, 9978 (2000); https://doi.org/10.1063/1.1323224
    """

    def __init__(self, tangent=None):
        self.kappa = None
        self.allpots = None
        self.tangent = tangent
        info(
            "NEBGradientMapper: Using %s tangent." % tangent,
            verbosity.low,
        )

    def bind(self, ens):
        # In principle, there is no need in dforces within the Mapper,
        # BUT dbeads are needed to calculate tangents for the endpoints,
        # and dforces are needed outside the Mapper to construct the "main" forces.
        self.dbeads = ens.beads.copy()
        self.dcell = ens.cell.copy()
        self.dforces = ens.forces.copy(self.dbeads, self.dcell)
        self.fixatoms = ens.fixatoms.copy()

        # Mask to exclude fixed atoms from 3N-arrays
        self.fixatoms_mask = np.ones(3 * ens.beads.natoms, dtype=bool)
        if len(ens.fixatoms) > 0:
            self.fixatoms_mask[3 * ens.fixatoms] = 0
            self.fixatoms_mask[3 * ens.fixatoms + 1] = 0
            self.fixatoms_mask[3 * ens.fixatoms + 2] = 0

        # Create reduced bead and force object
        self.rbeads = Beads(ens.beads.natoms, ens.beads.nbeads - 2)
        self.rbeads.q[:] = ens.beads.q[1:-1]
        self.rforces = ens.forces.copy(self.rbeads, self.dcell)

    def __call__(self, x):
        """Returns the potential for all beads and the gradient."""

        # Bead positions
        # Touch positions only if they have changed (to avoid triggering forces)
        # I need both dbeads and rbeads because of the endpoint tangents.
        if (self.rbeads.q[:, self.fixatoms_mask] != x).any():
            self.rbeads.q[:, self.fixatoms_mask] = x
        rbq = self.rbeads.q[:, self.fixatoms_mask]
        # We need all beads, but only reduced number is being updated.
        self.dbeads.q[1:-1, self.fixatoms_mask] = rbq
        bq = self.dbeads.q[:, self.fixatoms_mask]

        # Bead energies (needed for improved tangents)
        # Full-dimensional, but only reduced number is being updated.
        if self.allpots is None:
            info(
                "Calculating all beads once to get potentials on the endpoints",
                verbosity.medium,
            )
            self.allpots = dstrip(self.dforces.pots).copy()
            # We want to be greedy about force calls,
            # so we transfer from full beads to the reduced ones.
            tmp_f = self.dforces.f.copy()[1:-1]
            tmp_v = self.allpots.copy()[1:-1]
            self.rforces.transfer_forces_manual(
                new_q=[self.dbeads.q[1:-1]],
                new_v=[tmp_v],
                new_forces=[tmp_f],
            )
        be = self.allpots
        be[1:-1] = dstrip(self.rforces.pots).copy()

        # Forces
        rbf = dstrip(self.rforces.f).copy()[:, self.fixatoms_mask]

        # Number of images
        nimg = self.dbeads.nbeads

        # Number of atoms
        nat = self.dbeads.natoms - len(self.fixatoms)

        # Array for spring constants
        kappa = np.zeros(nimg)

        btau = np.zeros((nimg, 3 * nat), float)
        for ii in range(1, nimg - 1):
            d1 = bq[ii] - bq[ii - 1]  # tau minus
            d2 = bq[ii + 1] - bq[ii]  # tau plus

            if self.tangent == "plain":
                # Old implementation of NEB tangents
                btau[ii] = d1 / npnorm(d1) + d2 / npnorm(d2)
                btau[ii] /= npnorm(btau[ii])

            elif self.tangent == "improved":
                # Improved tangent estimate
                # J. Chem. Phys. 113, 9978 (2000) https://doi.org/10.1063/1.1323224

                # Energy of images: (ii+1) < (ii) < (ii-1)
                if be[ii + 1] < be[ii] < be[ii - 1]:
                    btau[ii] = d1
                # Energy of images (ii-1) < (ii) < (ii+1)
                elif be[ii - 1] <= be[ii] <= be[ii + 1]:
                    btau[ii] = d2
                # Energy of image (ii) is a minimum or maximum
                else:
                    maxpot = max(abs(be[ii + 1] - be[ii]), abs(be[ii - 1] - be[ii]))
                    minpot = min(abs(be[ii + 1] - be[ii]), abs(be[ii - 1] - be[ii]))

                    if be[ii + 1] >= be[ii - 1]:
                        btau[ii] = d2 * maxpot + d1 * minpot
                    elif be[ii + 1] < be[ii - 1]:
                        btau[ii] = d2 * minpot + d1 * maxpot
                btau[ii] /= npnorm(btau[ii])

            else:
                softexit.trigger(
                    status="bad",
                    message="Error: unknown tangent kind %s." % self.tangent,
                )

        # if mode == "variablesprings":
        #    if mode == "ci":
        #        # Climbing NEB term. Choose highest energy bead
        #        # after 5 (arbitrary) iterations
        #        if step >= 5:
        #            imax = np.argmax(be)
        #            bf[imax] = bf[imax] - 2 * np.dot(bf[imax], btau[imax]) * btau[imax]
        #
        #    # Determine variable spring constants
        #    kappa = np.zeros(nimg)
        #    ei = np.zeros(nimg)
        #    emax = np.amax(be)
        #    eref = max(be[0], be[nimg])
        #    kappamax = self.spring["kappa_max"]
        #    kappamin = self.spring["kappa_min"]
        #    deltakappa = kappamax - kappamin
        #    for ii in range(1, nimg - 1):
        #    ei[ii] = max(be[ii], be[ii - 1])
        #    if ei[j] > eref:
        #        kappa[ii] = kappamax - deltakappa * ((emax - ei[ii]) / (emax - eref))
        #    else:
        #        kappa[ii] = kappamin
        #
        # else:
        #     kappa.fill(self.kappa)

        # Array of spring constants; all are equal
        kappa.fill(self.kappa)

        # Get perpendicular forces
        for ii in range(1, nimg - 1):
            rbf[ii - 1] -= np.dot(rbf[ii - 1], btau[ii]) * btau[ii]

        for ii in range(1, nimg):
            print(
                "Bead %2i, distance to previous bead:   %f"
                % (ii, npnorm(bq[ii] - bq[ii - 1]))
            )

        # Adds the spring forces
        for ii in range(1, nimg - 1):
            # Old implementation (simple tangents, single kappa)
            # rbf[ii] += (
            #     kappa[ii]
            #     * btau[ii]
            #     * np.dot(btau[ii], (bq[ii + 1] + bq[ii - 1] - 2 * bq[ii]))
            # )

            # Improved tangent implementation
            # Eq. 12 in J. Chem. Phys. 113, 9978 (2000):
            rbf[ii - 1] += (
                kappa[ii]
                * (npnorm(bq[ii + 1] - bq[ii]) - npnorm(bq[ii] - bq[ii - 1]))
                * btau[ii]
            )

        # Including spring energy into potential (experimental)
        # The reason for having such term is that many optimization algorithms
        # need energy for line search or TRM update. But it doesn't work well.
        # for ii in range(1, nimg - 1):
        #     dl = (npnorm(bq[ii + 1] - bq[ii]) - npnorm(bq[ii] - bq[ii - 1]))
        #     be[ii] += 0.5 * kappa[ii] * dl * dl

        # Return NEB total potential energy of the band and total force gradient
        e = be.sum()
        g = -rbf
        return e, g


class NEBClimbGrMapper(object):
    """Creation of the multi-dimensional function that will be minimized.
    Functional analog of a GradientMapper in geop.py

    Constructs climbing forces for a single node.
    Node index is determined inside bind()
    Fixed atoms are excluded via boolean mask. 1 = moving, 0 = fixed.
    """

    def __init__(self):
        pass

    def bind(self, ens):
        """Creates reduced Beads object in order to calculate forces
        only for the climbing bead, and binds it to beads
        """
        self.fixatoms = ens.fixatoms.copy()
        # A mask to exclude fixed atoms from 3N-arrays
        self.fixatoms_mask = np.ones(3 * ens.beads.natoms, dtype=bool)
        if len(ens.fixatoms) > 0:
            self.fixatoms_mask[3 * ens.fixatoms] = 0
            self.fixatoms_mask[3 * ens.fixatoms + 1] = 0
            self.fixatoms_mask[3 * ens.fixatoms + 2] = 0

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
            raise RuntimeError("NEBClimbGrMapper.q_prev or q_next is None.")

        info("@NEB_CLIMB: NEBClimbGrMapper.x.shape: %s" % str(x.shape), verbosity.debug)
        # Touch positions only if they have changed (to avoid triggering forces)
        if (self.rbeads.q[:, self.fixatoms_mask] != x).any():
            self.rbeads.q[:, self.fixatoms_mask] = x
        # Compared to NEBGradientMapper, I need to flatten here,
        # otherwise it's (1, 3N) instead of just 3N.
        rbq = self.rbeads.q[:, self.fixatoms_mask].flatten()

        # Reduced forces
        # We don't need copy() because flatten() implies copying.
        rbf = dstrip(self.rforces.f)[:, self.fixatoms_mask].flatten()

        # I think here it's better to use plain tangents.
        # Then we don't need energies of the neighboring beads.
        d1 = rbq - self.q_prev  # tau minus
        d2 = self.q_next - rbq  # tau plus
        tau = d1 / npnorm(d1) + d2 / npnorm(d2)
        tau /= npnorm(tau)

        # Inverting the force component along the path
        rbf -= 2 * np.dot(rbf, tau) * tau

        # Return the potential energy and gradient of the climbing bead
        e = dstrip(self.rforces.pot.copy())
        g = -rbf
        info("@NEB_CLIMB: NEBClimbGrMapper finished.", verbosity.debug)
        return e, g


# TODO
# class NEB_GPR(object):
#    """KF: I strongly suggest that someone implements and this method:
#    J. A. G. Torres et al, PRL 122, 156001 (2019),
#    https://doi.org/10.1103/PhysRevLett.122.156001
#
#    Given that the stability of NEB is so far quite bad
#    and small step size is desired, it sounds very promising, and the examples
#    shown in the SI of that paper are really practical.
#
#    Also, it would be nice to build a unified framework
#    for NEB, String and probably Instanton methods to work with GPR.
#    """


class NEBMover(Motion):
    """Nudged elastic band routine.

    Attributes:
        mode: minimizer to use for NEB
        biggest_step: maximum step size for BFGS/L-BFGS
        old_force: force from previous iteration (not NEB-projected)
        old_direction: direction from previous iteration
        old_nebpotential: bead potentials from previous iteration
        old_nebgradient: gradient wrt springs from previous iteration
        hessian_bfgs: Hessian for (damped) BFGS
        tolerances:
            energy: tolerance on change in energy for exiting minimization
            force: tolerance on force/change in force for exiting minimization
            position: tolerance and change in position for exiting minimization
        corrections_lbfgs: number of corrections to store for L-BFGS
        qlist_lbfgs: list of previous positions (x_n+1 - x_n) for L-BFGS
        glist_lbfgs: list of previous gradients (g_n+1 - g_n) for L-BFGS
        endpoints: flag for minimizing end images in NEB *** NOT YET IMPLEMENTED ***
        use_climb: flag for climbing image NEB
        stage: flag denoting current stage: "endpoints", "neb" or "climb", or "converged"
        tangents: "plain" or "improved"
        spring:
            varsprings: T/F for variable spring constants
            kappa: single spring constant if varsprings is F
            kappamax: max spring constant if varsprings is T *** NOT YET IMPLEMENTED ***
            kappamin: min spring constant if varsprings is T *** NOT YET IMPLEMENTED ***
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
        old_nebpotential=np.zeros(0, float),
        old_nebgradient=np.zeros(0, float),
        old_direction=np.zeros(0, float),
        hessian_bfgs=np.eye(0),
        tolerances={"energy": 1e-5, "force": 1e-5, "position": 1e-5},
        corrections_lbfgs=5,
        qlist_lbfgs=np.zeros(0, float),
        glist_lbfgs=np.zeros(0, float),
        dtmax_fire=1.0,
        v_fire=np.zeros(0, float),
        alpha_fire=0.1,
        N_down_fire=0,
        N_up_fire=0,
        dt_fire=0.1,
        endpoints=(False, "bfgs"),
        spring={"varsprings": False, "kappa": 1.0, "kappamax": 1.5, "kappamin": 0.5},
        tangent=None,
        scale_lbfgs=2,
        stage="neb",
        use_climb=False,
        climb_bead=-1,
    ):
        """Initialises NEBMover.

        Args:
           fixcom: An optional boolean which decides whether the centre of mass
              motion will be constrained or not. Defaults to False.
        """

        super(NEBMover, self).__init__(fixcom=fixcom, fixatoms=fixatoms)

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
                verbosity.medium,
            )
        self.old_x = old_coord  # always in "masked" dimension
        self.nebpot = old_nebpotential
        self.nebgrad = old_nebgradient
        self.d = old_direction
        self.full_f = full_force  # physical forces of ALL beads even in climb stage
        self.full_v = full_pots  # potentials of ALL beads
        self.corrections = corrections_lbfgs
        self.qlist = qlist_lbfgs
        self.glist = glist_lbfgs
        self.scale = scale_lbfgs
        self.endpoints = endpoints
        self.spring = spring
        self.tangent = tangent
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

        self.nebgm = NEBGradientMapper(tangent=self.tangent)
        self.climbgm = NEBClimbGrMapper()

    def bind(self, ens, beads, nm, cell, bforce, prng, omaker):
        super(NEBMover, self).bind(ens, beads, nm, cell, bforce, prng, omaker)
        # I check dimensionality of the hessian in the beginning
        # and reduce it if needed, because someone may want to provide
        # existing hessian of the full system.
        if self.mode == "damped_bfgs":
            n_activedim = beads.q[0].size - len(self.fixatoms) * 3
            if self.stage == "neb":
                if self.hessian.size == (n_activedim * (beads.nbeads - 2)) ** 2:
                    # Desired dimension
                    info(
                        " @NEBMover: NEB Hessian is treated as reduced "
                        "and masked already, according to its size.",
                        verbosity.high,
                    )
                elif self.hessian.size == beads.q.size * beads.q.size:
                    info(
                        " @NEBMover: NEB Hessian has full-dimensional size, "
                        "will be reduced and masked by fixatoms.",
                        verbosity.high,
                    )
                    # First trim the endpoints, then mask fixed atoms
                    self.hessian = (
                        self.hessian[
                            self.beads.natoms : -self.beads.natoms,
                            self.beads.natoms : -self.beads.natoms,
                        ]
                    )[
                        np.ix_(
                            np.tile(self.nebgm.fixatoms_mask, self.beads.nbeads - 2),
                            np.tile(self.nebgm.fixatoms_mask, self.beads.nbeads - 2),
                        )
                    ]
                elif self.hessian.size == 0:
                    info(
                        " @NEBMover: No Hessian provided, starting from unity matrix",
                        verbosity.high,
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
                        " @NEBMover: NEB climbing Hessian is treated as masked one,"
                        " according to its size.",
                        verbosity.high,
                    )
                elif self.hessian.size == (beads.q[0].size * beads.q[0].size):
                    info(
                        " @NEBMover: NEB climbing Hessian has size (3natoms x 3natoms), "
                        "will be masked by fixatoms.",
                        verbosity.high,
                    )
                    self.hessian = self.hessian[
                        np.ix_(self.climbgm.fixatoms_mask, self.climbgm.fixatoms_mask)
                    ]
                if self.hessian.size == 0:
                    info(
                        " @NEBMover: No climbing Hessian provided, "
                        "starting from unity matrix",
                        verbosity.high,
                    )
                    self.hessian = np.eye(n_activedim, dtype=float)
                else:
                    raise ValueError("Hessian size does not match system size.")

        if len(self.fixatoms) == len(self.beads[0]):
            softexit.trigger(
                status="bad",
                message="WARNING: all atoms are fixed, geometry won't change. Exiting simulation.",
            )

        self.nebgm.bind(self)
        self.climbgm.bind(self)

    def init_climb(self):
        """Operations to resize hessian for the climbing bead(s),
        get neighboring beads and other things necessary for climbing.
        Returns:
          cl_indx - index of the bead for climbing
        """
        cl_indx = np.argmax(self.forces.pots)
        info(
            " @NEB: Initializing climbing. Climbing bead: %i." % cl_indx,
            verbosity.medium,
        )
        if cl_indx in [0, self.beads.nbeads - 1]:
            softexit.trigger(
                status="bad", message="ERROR: climbing bead is the endpoint."
            )
        self.climbgm.rbeads.q[:] = self.beads.q[cl_indx]
        self.climbgm.q_prev[:] = self.beads.q[cl_indx - 1, self.climbgm.fixatoms_mask]
        self.climbgm.q_next[:] = self.beads.q[cl_indx + 1, self.climbgm.fixatoms_mask]
        info(
            "q_prev.shape: %s,    q_next.shape: %s"
            % (str(self.climbgm.q_prev.shape), str(self.climbgm.q_next.shape)),
            verbosity.debug,
        )
        n_activedim = self.beads.q[0].size - len(self.fixatoms) * 3
        if self.hessian.shape != (n_activedim, n_activedim):
            self.hessian = np.eye(n_activedim)

        info(" @NEB_CLIMB: calling NEBClimbGrMapper first time.", verbosity.debug)
        self.nebpot, self.nebgrad = self.climbgm(
            self.beads.q[cl_indx, self.climbgm.fixatoms_mask]
        )

        # If climbing is launched from geometries rather than i-pi restart,
        # these need to be initialized
        if self.full_f.shape == (0,) or self.full_v.shape == (0,):
            self.full_f = dstrip(self.nebgm.dforces.f).copy()
            self.full_v = dstrip(self.nebgm.dforces.pots).copy()
            info(
                " @NEB_CLIMB: full_f and full_v are initialized in init_climb()."
                "\n\tself.full_f.shape: %s"
                "\n\tself.full_v.shape: %s"
                % (str(self.full_f.shape), str(self.full_v.shape)),
                verbosity.debug,
            )

        # Initialize FIRE parameters
        self.v = -self.a * self.nebgrad
        self.a = 0.1
        self.N_dn = 0
        self.N_up = 0
        self.dt_fire = 0.1

        self.old_x = dstrip(self.beads.q[cl_indx, self.climbgm.fixatoms_mask]).copy()
        self.stage = "climb"
        return cl_indx

    def adjust_big_step(self, quality):
        """We try to make big_step brave enough, but if relaxation goes wrong,
        we decrease the maximal allowed step."""

        if quality > 0.8:
            self.big_step *= 1.3
            info(
                " @NEBMover: increasing big_step to %.6f bohr." % self.big_step,
                verbosity.debug,
            )
        elif quality > 0.6:
            self.big_step *= 1.1
            info(
                " @NEBMover: increasing big_step to %.6f bohr." % self.big_step,
                verbosity.debug,
            )
        else:
            self.big_step *= 0.5
            # In no case we need big_step going to zero completely
            if self.big_step <= 0.01:
                self.big_step = 0.02
            info(
                " @NEBMover: Step direction far from nebgrad/climbgrad, "
                "reducing big_step to %.6f bohr." % self.big_step,
                verbosity.debug,
            )

    def build_precon_hessian(self):
        """Preconditioner for the Hessian with Lindh for intra-bead parts
        and explicit spring terms between neighboring beads"""
        # TODO
        pass

    def step(self, step=None):
        """Does one simulation time step.
        Dimensionality is reduced in the very beginning.
        """

        info(" @NEB STEP %d, stage: %s" % (step, self.stage), verbosity.debug)

        n_activedim = self.beads.q[0].size - len(self.fixatoms) * 3

        # Check if we restarted a converged calculation (by mistake)
        if self.stage == "converged":
            softexit.trigger(
                status="success",
                message="NEB has already converged. Exiting simulation.",
            )

        # First, optimization of endpoints, if required
        if self.endpoints["optimize"] and self.stage == "endpoints":
            # TODO
            softexit.trigger(
                status="bad",
                message="Optimization of endpoints in NEB is not implemented yet.",
            )

        # Endpoints are optimized or optimization is not required
        elif self.stage == "neb":
            # Fetch spring constants
            if self.spring["varsprings"] == True:
                softexit.trigger(
                    status="bad",
                    message="Variable springs in NEB are not implemented yet.",
                )
            self.nebgm.kappa = self.spring["kappa"]

            self.ptime = self.ttime = 0
            self.qtime = -time.time()

            if self.mode == "damped_bfgs":
                # All BFGS-family algorithms would have similar structure, but currently
                # BFGS, LBFGS and BFGSTRM are not suited for NEB because they use energy,
                # which is ill-defined in NEB.
                if step == 0:  # TODO add a condition when after the endpoints.
                    # Initialize direction to the steepest descent direction
                    info(" @NEB: calling NEBGradientMapper at step 0.", verbosity.debug)
                    self.nebpot, self.nebgrad = self.nebgm(
                        self.beads.q[1:-1, self.nebgm.fixatoms_mask]
                    )
                    info(
                        " @NEB: NEBGradientMapper returned nebpot and nebgrad.",
                        verbosity.debug,
                    )

                    # Store old bead positions
                    self.old_x = dstrip(
                        self.beads.q[1:-1, self.nebgm.fixatoms_mask]
                    ).copy()

                    # At step 0, we also need to store full forces and pots.
                    self.full_f = dstrip(self.nebgm.dforces.f).copy()
                    self.full_v = dstrip(self.nebgm.dforces.pots).copy()

                    # With multiple stages, the size of the hessian is different
                    # at each stage, therefore we check.
                    if self.hessian.shape != (
                        (self.beads.nbeads - 2) * n_activedim,
                        (self.beads.nbeads - 2) * n_activedim,
                    ):
                        print("Dimensions of the Hessian and of the beads:")
                        print((self.hessian.shape, self.beads.q.shape))
                        softexit.trigger(
                            status="bad",
                            message="Hessian not initialized correctly in NEB.",
                        )

                # Self instances will be updated in the optimizer, so we store the copies.
                # old_nebpot is used later as a convergence criterion.
                old_nebpot = self.nebpot.copy()
                # old_nebgrad = self.nebgrad.copy()

                info(" @NEB: before Damped_BFGS() call", verbosity.debug)
                print("self.old_x.shape: %s" % str(self.old_x.shape))
                print("self.nebgrad.shape: %s" % str(self.nebgrad.shape))
                print("self.hessian.shape: %s" % str(self.hessian.shape))
                quality = Damped_BFGS(
                    x0=self.old_x.copy(),
                    fdf=self.nebgm,
                    fdf0=(self.nebpot, self.nebgrad),
                    hessian=self.hessian,
                    big_step=self.big_step,
                )
                info(" @NEB: after Damped_BFGS() call", verbosity.debug)

                self.adjust_big_step(quality)

                # tmp printout, remove after neb is finished
                # if step % 100 == 0:
                #     eigvals, eigvecs = np.linalg.eigh(self.hessian)
                #     idx = eigvals.argsort()
                #     np.savetxt("eigvecs.s%04d.dat" % step, eigvecs[:, idx])
                #     np.savetxt("eigvals.s%04d.dat" % step, eigvals[idx])

            elif self.mode == "fire":
                # Only initialize velocity for fresh start, not for RESTART
                if step == 0 and self.v.size == 0:
                    info(
                        " @NEB: calling NEBGradientMapper at step 0 by FIRE",
                        verbosity.debug,
                    )
                    self.nebpot, self.nebgrad = self.nebgm(
                        self.beads.q[1:-1, self.nebgm.fixatoms_mask]
                    )
                    self.old_x = dstrip(
                        self.beads.q[1:-1, self.nebgm.fixatoms_mask].copy()
                    )

                    self.v = -self.a * self.nebgrad

                    # At step 0, we also need to store full forces and pots
                    self.full_f = dstrip(self.nebgm.dforces.f).copy()
                    self.full_v = dstrip(self.nebgm.dforces.pots).copy()

                # Store potential and force gradient for convergence criterion
                old_nebpot = self.nebpot.copy()
                # old_nebgrad = self.nebgrad.copy()
                info(" @NEB: using FIRE", verbosity.debug)
                info(" @FIRE velocity: %s" % str(npnorm(self.v)), verbosity.debug)
                info(" @FIRE alpha: %s" % str(self.a), verbosity.debug)
                info(" @FIRE N down: %s" % str(self.N_dn), verbosity.debug)
                info(" @FIRE N up: %s" % str(self.N_up), verbosity.debug)
                info(" @FIRE dt: %s" % str(self.dt_fire), verbosity.debug)
                self.v, self.a, self.N_dn, self.N_up, self.dt_fire = FIRE(
                    x0=self.old_x.copy(),
                    fdf=self.nebgm,
                    fdf0=(self.nebpot, self.nebgrad),
                    v=self.v,
                    a=self.a,
                    N_dn=self.N_dn,
                    N_up=self.N_up,
                    dt=self.dt_fire,
                    dtmax=self.dtmax,
                )
                info(" @NEB: after FIRE call")

            # TODO: Routines for L-BFGS, SD, CG
            else:
                softexit.trigger(
                    status="bad",
                    message="Try 'damped_bfgs' or 'fire'. Other algorithms are not implemented for NEB.",
                )

            # Update positions
            self.beads.q[:] = self.nebgm.dbeads.q
            info(" @NEB: bead positions transferred from mapper.", verbosity.debug)

            # Recalculation won't be triggered because the position is the same.
            self.nebpot, self.nebgrad = self.nebgm(
                self.beads.q[1:-1, self.nebgm.fixatoms_mask]
            )
            info(" @NEB: NEB forces transferred from mapper.", verbosity.debug)

            # dx = current position - previous position.
            # Use to determine converged minimization
            dx = np.amax(
                np.abs(self.beads.q[1:-1, self.nebgm.fixatoms_mask] - self.old_x)
            )

            # Store old positions
            self.old_x[:] = self.beads.q[1:-1, self.nebgm.fixatoms_mask]

            # This transfers forces from the mapper to the "main" beads,
            # so that recalculation won't be triggered after the step.
            # I need to keep forces up to date, because I should be able to output
            # full potentials and full forces.
            tmp_f = self.full_f.copy()
            tmp_f[1:-1] = self.nebgm.rforces.f
            tmp_v = self.full_v
            tmp_v[1:-1] = self.nebgm.rforces.pots
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

            # Full-dimensional forces are not really needed all the time,
            # but it's easier to have it. They are needed in reduced-beads modes,
            # i.e. in "endpoints" and "climb".
            self.full_f = dstrip(self.forces.f)
            self.full_v = dstrip(self.forces.pots)

            info(
                " @NEB: remaining max force component: {}".format(
                    np.amax(np.abs(self.nebgrad))
                )
            )
            info(" @NEB: max delta x: {}".format(dx))

            self.qtime += time.time()

            # Check convergence criteria
            if (
                (
                    np.amax(np.abs(self.nebpot - old_nebpot))
                    / (self.beads.nbeads * self.beads.natoms)
                    <= self.tolerances["energy"]
                )
                and (np.amax(np.abs(self.nebgrad)) <= self.tolerances["force"])
                and (dx <= self.tolerances["position"])
            ):
                decor = 60 * "=" + "\n"
                info(
                    decor
                    + " @NEB: path optimization converged. Step: %i\n" % step
                    + decor,
                    verbosity.medium,
                )

                # Set climbing stage indicator after convergence of the path
                if self.use_climb:
                    self.stage = "climb"
                else:
                    self.stage = "converged"
                    softexit.trigger(
                        status="success",
                        message="NEB finished successfully at STEP %i." % step,
                    )

            else:
                info(
                    " @NEB: Not converged, deltaEnergy = %.8f, tol = %.8f per atom"
                    % (
                        np.amax(np.abs(self.nebpot - old_nebpot))
                        / (self.beads.nbeads * self.beads.natoms),
                        self.tolerances["energy"],
                    ),
                    verbosity.debug,
                )
                info(
                    " @NEB: Not converged, nebgrad = %.8f, tol = %f"
                    % (np.amax(np.abs(self.nebgrad)), self.tolerances["force"]),
                    verbosity.debug,
                )
                info(
                    " @NEB: Not converged, deltaX = %.8f, tol = %.8f"
                    % (dx, self.tolerances["position"]),
                    verbosity.debug,
                )

        # ============================== C L I M B ==============================

        # Climbing image optimization
        elif self.stage == "climb":
            self.ptime = self.ttime = 0
            self.qtime = -time.time()

            # We need to initialize climbing once
            if np.all(self.climbgm.q_prev == 0.0) or np.all(self.climbgm.q_next == 0.0):
                self.cl_indx = self.init_climb()

            if self.mode == "damped_bfgs":
                # BFGS-family algorithms

                # Self instances will be updated in the optimizer, so we store the copies.
                # old_nebpot is used later as a convergence criterion.
                old_nebpot = self.nebpot.copy()
                # old_nebgrad = self.nebgrad.copy()

                # if self.mode == "damped_bfgs":
                info(" @NEB_CLIMB: before Damped_BFGS() call", verbosity.debug)
                print("self.old_x.shape: %s" % str(self.old_x.shape))
                print("self.nebgrad.shape: %s" % str(self.nebgrad.shape))
                print("self.hessian.shape: %s" % str(self.hessian.shape))
                quality = Damped_BFGS(
                    x0=self.old_x.copy(),
                    fdf=self.climbgm,
                    fdf0=(self.nebpot, self.nebgrad),
                    hessian=self.hessian,
                    big_step=self.big_step,
                )
                info(" @NEB_CLIMB: after Damped_BFGS() call", verbosity.debug)

                self.adjust_big_step(quality)

            elif self.mode == "fire":
                # FIRE algorithm

                # Self instances will be updated in the optimizer, so we store the copies.
                # old_nebpot is used later as a convergence criterion.
                old_nebpot = self.nebpot.copy()
                # old_nebgrad = self.nebgrad.copy()
                info(" @NEB: using FIRE", verbosity.debug)
                info(" @FIRE velocity: %s" % str(npnorm(self.v)), verbosity.debug)
                info(" @FIRE alpha: %s" % str(self.a), verbosity.debug)
                info(" @FIRE N down: %s" % str(self.N_dn), verbosity.debug)
                info(" @FIRE N up: %s" % str(self.N_up), verbosity.debug)
                info(" @FIRE dt: %s" % str(self.dt_fire), verbosity.debug)
                self.v, self.a, self.N_dn, self.N_up, self.dt_fire = FIRE(
                    x0=self.old_x.copy(),
                    fdf=self.climbgm,
                    fdf0=(self.nebpot, self.nebgrad),
                    v=self.v,
                    a=self.a,
                    N_dn=self.N_dn,
                    N_up=self.N_up,
                    dt=self.dt_fire,
                    dtmax=self.dtmax,
                )

            # TODO: Routines for L-BFGS, SD, CG, ...
            else:
                softexit.trigger(
                    status="bad",
                    message="Try damped_bfgs or fire, other algorithms are not implemented for NEB.",
                )

            # Update positions
            self.beads.q[self.cl_indx] = self.climbgm.rbeads.q
            info(" @NEB_CLIMB: climb beads positions updated.", verbosity.debug)

            self.nebpot, self.nebgrad = self.climbgm(
                self.beads.q[self.cl_indx, self.climbgm.fixatoms_mask]
            )

            # Use to determine converged minimization
            # max movement
            dx = np.amax(
                np.abs(
                    self.beads.q[self.cl_indx, self.climbgm.fixatoms_mask] - self.old_x
                )
            )

            # Store old positions
            self.old_x[:] = self.beads.q[self.cl_indx, self.climbgm.fixatoms_mask]

            # This transfers forces from the ClimbMapper to the "main" beads,
            # so that recalculation won't be triggered after the step.
            tmp_f = self.full_f.copy()
            tmp_v = self.full_v.copy()
            tmp_f[self.cl_indx, self.climbgm.fixatoms_mask] = self.nebgrad
            tmp_v[self.cl_indx] = self.nebpot
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
                    np.amax(np.abs(self.nebpot - old_nebpot)) / self.beads.natoms
                    <= self.tolerances["energy"]
                )
                and (np.amax(np.abs(self.nebgrad)) <= self.tolerances["force"])
                and (dx <= self.tolerances["position"])
            ):
                decor = 60 * "=" + "\n"
                info(
                    decor
                    + " @NEB_CLIMB: optimization converged. Step: %i\n" % step
                    + decor,
                    verbosity.medium,
                )
                self.stage = "converged"
                softexit.trigger(
                    status="success", message="NEB_CLIMB finished successfully."
                )

            else:
                info(
                    " @NEB_CLIMB: Not converged, deltaEnergy = %.8f, tol = %.8f per atom"
                    % (
                        np.amax(np.abs(self.nebpot - old_nebpot)) / self.beads.natoms,
                        self.tolerances["energy"],
                    ),
                    verbosity.debug,
                )
                info(
                    " @NEB_CLIMB: Not converged, climbgrad = %.8f, tol = %f"
                    % (np.amax(np.abs(self.nebgrad)), self.tolerances["force"]),
                    verbosity.debug,
                )
                info(
                    " @NEB_CLIMB: Not converged, deltaX = %.8f, tol = %.8f"
                    % (dx, self.tolerances["position"]),
                    verbosity.debug,
                )

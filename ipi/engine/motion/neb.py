"""Holds the algorithms to perform nudged elastic band calculations.

Algorithms implemented by Michele Ceriotti and Benjamin Helfrecht, 2015
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

import numpy as np
from numpy.linalg import norm as npnorm
import time

from ipi.engine.motion import Motion
from ipi.utils.depend import *
from ipi.utils.softexit import softexit
from ipi.utils.mintools import L_BFGS, BFGS, BFGSTRM, Damped_BFGS, min_brent_neb
from ipi.utils.messages import verbosity, info
from ipi.engine.beads import Beads

np.set_printoptions(threshold=10000, linewidth=1000)

__all__ = ["NEBMover"]


# TODO: Do not shout :-)
# NOTE: CURRENTLY, NEB DOESN'T WORK.
#       IF SD, CG, OR BFGS OPTIONS ARE NOT DESIRED, CONSIDER ELIMINATING
#       NEBLineMapper AND THE RELEVANT BLOCKS IN NEBMover. IF THESE OPTIONS ARE DESIRED,
#       THE INFRASTRUCTURE IS PRESENT BUT MUST BE DEBUGGED AND MADE CONSISTENT WITH
#       THAT PRESENT IN NEBGradientMapper (TO REMOVE REMAINING ERRORS IN COMPUTATION).
#       CLIMBING IMAGE AND VARIABLE SPRING CONSTANTS HAVE NOT YET BEEN IMPLEMENTED, BUT
#       THE GENERIC INFRASTRUCTURE IS PRESENT (SEE COMMENTED BLOCKS IN NEBLineMapper AND
#       NEBGradientMapper). REQUIRES REARRANGEMENT AND DEBUGGING.
#
#       THIS NEB IMPLEMENTATION USES THE "IMPROVED TANGENTS" OF HENKELMAN AND JONSSON, 2000.
#       THE "OLD IMPLEMENTATION" IS PRESERVED IN COMMENTS


class NEBLineMapper(object):
    """Creation of the one-dimensional function that will be minimized

    Attributes:
        x0: initial position
        d: move direction
        kappa: spring constants
        first: flag indicating first iteration of simulation
    """

    def __init__(self):
        self.x0 = None
        self.d = None
        self.kappa = None
        self.first = True

    def bind(self, ens):
        self.dbeads = ens.beads.copy()
        self.dcell = ens.cell.copy()
        self.dforces = ens.forces.copy(self.dbeads, self.dcell)

    def set_dir(self, x0, mdir):
        self.x0 = x0.copy()
        self.d = mdir.copy() / npnorm(mdir.flatten())
        if self.x0.shape != self.d.shape:
            raise ValueError(
                "Incompatible shape of initial value and displacement direction"
            )

    def __call__(self, x):
        if self.first is True:
            self.dbeads.q = x
        else:
            self.dbeads.q = self.x0 + self.d * x

        # List of atom/bead positions
        bq = dstrip(self.dbeads.q).copy()

        # List of forces
        bf = dstrip(self.dforces.f).copy()

        # List of bead energies
        # be = dstrip(self.dforces.pots).copy()

        # Number of images
        nimg = self.dbeads.nbeads

        # Number of atoms
        nat = self.dbeads.natoms

        kappa = np.zeros(nimg)

        # Get tangents, end images are distinct, fixed, pre-relaxed configurations
        btau = np.zeros((nimg, 3 * nat), float)
        for ii in range(1, nimg - 1):
            d1 = bq[ii] - bq[ii - 1]  # tau minus
            d2 = bq[ii + 1] - bq[ii]  # tau plus

            # "Old" implementation of NEB
            btau[ii] = d1 / npnorm(d1) + d2 / npnorm(d2)
            btau[ii] *= 1.0 / npnorm(btau)

        # Energy of images: (ii+1) < (ii) < (ii-1)
        #            if (be[ii + 1] < be[ii]) and (be[ii] < be[ii - 1]):
        #                btau[ii] = d2
        #
        # Energy of images (ii-1) < (ii) < (ii+1)
        #            elif (be[ii - 1] < be[ii]) and (be[ii] < be[ii + 1]):
        #                btau[ii] = d1
        #
        # Energy of image (ii) is a minimum or maximum
        #            else:
        #                maxpot = max(be[ii + 1] - be[ii], be[ii - 1], be[ii])
        #                minpot = min(be[ii + 1] - be[ii], be[ii - 1], be[ii])
        #
        #                if be[ii + 1] < be[ii - 1]:
        #                    btau[ii] = d1 * minpot + d2 * maxpot
        #
        #                elif be[ii - 1] < be[ii + 1]:
        #                    btau[ii] = d1 * maxpot + d2 * minpot
        #
        #                else:
        #                    print "Error in NEB tangents: Energy of images are equal"
        #
        #            btau[ii] *= 1.0 / npnorm(btau)

        # if mode == "variablesprings": #TODO: input option for variable spring mode

        #        if mode == "ci":
        #
        # Climbing NEB term. Choose highest energy bead after 5 (arbitrary) iterations
        #            if step >= 5:
        #                imax = np.argmax(be)
        #                bf[imax] = bf[imax] - 2 * np.dot(bf[imax], btau[imax]) * btau[imax]
        #
        # Determine variable spring constants
        # kappa = np.zeros(nimg)
        # ei = np.zeros(nimg)
        # emax = np.amax(be)
        # eref = max(be[0], be[nimg])
        # kappamax = self.kappa_max
        # kappamin = self.kappa_min #TODO: input options for max and min spring constant
        # deltakappa = kappamax - kappamin
        # for ii in range(1, nimg - 1):
        # ei[ii] = max(be[ii], be[ii - 1])
        # if ei[j] > eref:
        # kappa[ii] = kappamax - deltakappa * ((emax - ei[ii]) / (emax - eref))
        # else:
        # kappa[ii] = kappamin
        #
        #        else:
        #            kappa.fill(self.kappa)
        #
        #
        # get perpendicular forces
        #            for ii in range(1, nimg - 1):
        #                bf[ii] = bf[ii] - np.dot(bf[ii], btau[ii]) * btau[ii]
        #
        # adds the spring forces
        #            for ii in range(1, nimg - 1):
        #                bf[ii] += kappa[ii] * btau[ii] * np.dot(btau[ii], (bq[ii + 1] + bq[ii - 1] - 2 * bq[ii]))

        kappa.fill(self.kappa)

        # get perpendicular forces
        for ii in range(1, nimg - 1):
            bf[ii] = bf[ii] - np.dot(bf[ii], btau[ii]) * btau[ii]

        # adds the spring forces
        for ii in range(1, nimg - 1):
            bf[ii] += (
                kappa[ii]
                * btau[ii]
                * np.dot(btau[ii], (bq[ii + 1] + bq[ii - 1] - 2 * bq[ii]))
            )

        # For first iteration, move in direction of the force
        if self.first is True:
            self.d = bf
            self.first = False

        force = bf
        g = -np.vdot(bf, self.d)
        g = abs(g)

        # Return NEB forces and gradient modulus
        linefunc = (force, g)
        return linefunc


class NEBGradientMapper(object):
    """ Creation of the multi-dimensional function that will be minimized.
        Functional analog of a GradientMapper in geop.py

        Fixed atoms are excluded via boolean mask. 1 = moving, 0 = fixed.

    Attributes:
        kappa: spring constants"""

    def __init__(self):
        self.kappa = None

    def bind(self, ens):
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

        # This needs to be sorted out. In principle, there is no need in dforces,
        # but dbeads are needed to calculate tangents for the endpoints.
        # # Create reduced bead and force object and evaluate forces
        # self.reduced_b = Beads(ens.beads.natoms, ens.beads.nbeads)
        # self.reduced_b.q[:] = ens.beads.q[1:-1]
        # self.reduced_forces = self.dforces.copy(self.reduced_b, self.dcell)
        # # Do I really need this?
        # self.rpots = self.reduced_forces.pots  # reduced energy
        # self.rforces = self.reduced_forces.f  # reduced gradient

    def get_tangents(bq=None, be=None):
        """ Compute tangents.
            The end images are distinct, fixed, pre-relaxed configurations.
            Arguments:
              bq - bead coordinates (nbeads, 3*natoms)
              be - bead energies (nbeads)
            Returns:
              btau - array of tangents (nbeads, 3)
        """
        btau = np.zeros((nimg, 3 * nat), float)
        for ii in range(1, nimg - 1):
            d1 = bq[ii] - bq[ii - 1]  # tau minus
            d2 = bq[ii + 1] - bq[ii]  # tau plus

            # Old implementation of NEB tangents
            # btau[ii] = d1 / npnorm(d1) + d2 / npnorm(d2)
            # btau[ii] /= npnorm(btau[ii])

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
        return btau

    def __call__(self, x):
        """ Returns the potential for all beads and the gradient.
        """

        # Bead positions
        # Touch positions only if they have changed (to avoid triggering forces)
        # I need both dbeads and reduced_b because of the endpoint tangents.
        if (self.dbeads.q[:, self.fixatoms_mask] != x).any():
            self.dbeads.q[:, self.fixatoms_mask] = x
        bq = self.dbeads.q[:, self.fixatoms_mask]
#        if (self.reduced_b.q[:, self.fixatoms_mask] != x).any():
#            self.reduced_b.q[:, self.fixatoms_mask] = x
#        rbq = self.reduced_b.q[:, self.fixatoms_mask]

        # Forces
        bf = self.dforces.f.copy()[:, self.fixatoms_mask]
        # Zeroing endpoint forces
        bf[0,:] = bf[-1,:] = 0.

        # Bead energies (needed for improved tangents)
        be = self.dforces.pots.copy()

        # Number of images
        nimg = self.dbeads.nbeads

#        # Number of atoms
#        nat = self.dbeads.natoms - len(self.fixatoms)

        # Array for spring constants
        kappa = np.zeros(nimg)

        btau = self.get_tangents(bq, be)

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
            bf[ii] -= np.dot(bf[ii], btau[ii]) * btau[ii]

        for ii in range(1, nimg):
            print(
                "Bead %2i, distance to previous:   %f"
                % (ii, npnorm(bq[ii] - bq[ii - 1]))
            )

        # Adds the spring forces
        for ii in range(1, nimg - 1):
            # Old implementation (simple tangents, single kappa)
            # bf[ii] += (
            #     kappa[ii]
            #     * btau[ii]
            #     * np.dot(btau[ii], (bq[ii + 1] + bq[ii - 1] - 2 * bq[ii]))
            # )

            # Improved tangent implementation
            # Eq. 12 in J. Chem. Phys. 113, 9978 (2000):
            bf[ii] += (
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

        # Return NEB forces (negative) and the total potential energy of the band
        e = be.sum()
        g = -bf
        return e, g


class NEBClimbGrMapper(object):
    """ Creation of the multi-dimensional function that will be minimized.
        Functional analog of a GradientMapper in geop.py

        Constructs climbing forces for a single node.
        Node index is determined inside bind()
        Fixed atoms are excluded via boolean mask. 1 = moving, 0 = fixed.
    """

    def __init__(self):
        self.q_prev = None  # Coordinates of the bead before the climbing one.
        self.q_next = None  #           -''-          after      -''-

    def bind(self, ens):
        """ Creates reduced Beads object in order to calculate forces
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
        self.reduced_b = Beads(ens.beads.natoms, 1)
        self.dcell = ens.cell.copy()
        self.reduced_forces = ens.forces.copy(self.reduced_b, self.dcell)
        self.rforces = self.reduced_forces.f  # reduced gradient

    def __call__(self, x):
        """ Returns climbing force for a climbing image.
        """
        if self.q_prev is None or self.q_next is None:
            raise RuntimeError("NEBClimbGrMapper.q_prev or .q_next is None.")

        # Touch positions only if they have changed (to avoid triggering forces)
        if (self.reduced_b.q[:, self.fixatoms_mask] != x).any():
            self.reduced_b.q[:, self.fixatoms_mask] = x
        rbq = self.reduced_b.q[:, self.fixatoms_mask]

        # Reduced forces
        rbf = self.rforces.copy()[:, self.fixatoms_mask]

        # I think here it's better to use plain tangents.
        # Then we don't need energies of the neighboring beads.
        d1 = rbq - self.q_prev  # tau minus
        d2 = self.q_next - rbq  # tau plus
        tau = d1 / npnorm(d1) + d2 / npnorm(d2)
        tau /= npnorm(tau)

        # Inverting the force component along the path
        rbf -= 2 * np.dot(rbf, tau) * tau

        # Return the potential energy and gradient of the climbing bead
        e = self.reduced_forces.pot
        g = -rbf
        return e, g

class NEBMover(Motion):
    """Nudged elastic band routine.

    Attributes:
        mode: minimizer to use for NEB
        biggest_step: maximum step size for BFGS/L-BFGS
        old_force: force from previous iteration (not NEB-projected)
        old_direction: direction from previous iteration
        old_nebpotential: bead potentials from previous iteration
        old_nebgradient: gradient wrt springs from previous iteration
        invhessian_bfgs: inverse Hessian for (damped) BFGS
        ls_options:
            tolerance: tolerance for exit of line search
            iter: maximum iterations for line search per MD step
            step: initial step size for SD/CG
            adaptive: flag for adaptive step size
        tolerances:
            energy: tolerance on change in energy for exiting minimization
            force: tolerance on force/change in force for exiting minimization
            position: tolerance and change in position for exiting minimization
        corrections_lbfgs: number of corrections to store for L-BFGS
        qlist_lbfgs: list of previous positions (x_n+1 - x_n) for L-BFGS
        glist_lbfgs: list of previous gradients (g_n+1 - g_n) for L-BFGS
        endpoints: flag for minimizing end images in NEB *** NOT YET IMPLEMENTED ***
        use_climb: flag for climbing image NEB
        stage: flag denoting current procedure: "endpoints", "neb" or "climb"
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
        old_force=np.zeros(0, float),
        old_nebpotential=np.zeros(0, float),
        old_nebgradient=np.zeros(0, float),
        old_direction=np.zeros(0, float),
        invhessian_bfgs=np.eye(0),
        ls_options={"tolerance": 1e-5, "iter": 100.0, "step": 1e-3, "adaptive": 1.0},
        tolerances={"energy": 1e-5, "force": 1e-5, "position": 1e-5},
        corrections_lbfgs=5,
        qlist_lbfgs=np.zeros(0, float),
        glist_lbfgs=np.zeros(0, float),
        tr_trm=np.zeros(0, float),
        endpoints=(False, "bfgs"),
        spring={"varsprings": False, "kappa": 1.0, "kappamax": 1.5, "kappamin": 0.5},
        scale_lbfgs=2,
        stage="neb",
        use_climb=False,
    ):
        """Initialises NEBMover.

        Args:
           fixcom: An optional boolean which decides whether the centre of mass
              motion will be constrained or not. Defaults to False.
        """

        super(NEBMover, self).__init__(fixcom=fixcom, fixatoms=fixatoms)

        # Optimization options
        self.ls_options = ls_options
        self.tolerances = tolerances
        self.mode = mode
        self.big_step = biggest_step
        if self.big_step > 0.5:
            softexit.trigger(" @DampedBFGS: ERROR: big_step is too big. "
                             "Damped_BFGS handles it differently from other algorithms, "
                             "see ipi/inputs/motion/neb.py:'biggest_step' to get an idea.\n"
                             "Current big_step value: %f." % big_step)
        self.old_x = old_coord
#        self.old_f = old_force
        self.nebpot = old_nebpotential
        self.nebgrad = old_nebgradient
        self.d = old_direction
        self.invhessian = invhessian_bfgs
        self.corrections = corrections_lbfgs
        self.qlist = qlist_lbfgs
        self.glist = glist_lbfgs
        self.tr_trm = tr_trm
        self.endpoints = endpoints
        self.spring = spring
        self.use_climb = use_climb
        self.scale = scale_lbfgs

        self.neblm = NEBLineMapper()
        self.nebgm = NEBGradientMapper()
        self.climbgm = NEBClimbGrMapper()


    def bind(self, ens, beads, nm, cell, bforce, prng, omaker):

        super(NEBMover, self).bind(ens, beads, nm, cell, bforce, prng, omaker)
#        if self.old_f.shape != beads.q.shape:
#            if self.old_f.shape == (0,):
#                self.old_f = np.zeros(beads.q.shape, float)
#            else:
#                raise ValueError("Force size does not match system size.")
#        if self.d.shape != beads.q.shape:
#            if self.d.shape == (0,):
#                self.d = np.zeros(beads.q.shape, float)
#            else:
#                raise ValueError("Direction size does not match system size.")
        # I don't reduce dimensionality of the invhessian from the beginning,
        # because someone may want to provide existing invessian of the full system.
        if (self.mode in ["damped_bfgs", "bfgs", "bfgstrm"]
            and self.invhessian.size != (beads.q.size * beads.q.size)
        ):
            if self.invhessian.size == 0:
                self.invhessian = np.eye(beads.q.size, dtype=float)
            else:
                raise ValueError("Inverse Hessian size does not match system size.")

        if len(self.fixatoms) == len(self.beads[0]):
            softexit.trigger(
                "WARNING: all atoms are fixed, geometry won't change. Exiting simulation"
            )

        self.nebgm.bind(self)
        self.climbgm.bind(self)


    def step(self, step=None):
        """ Does one simulation time step.
            Dimensionality is reduced in the very beginning.
        """

        info(" @NEB STEP %d" % step, verbosity.debug)

        # First, optimization of endpoints, if required
        if self.endpoints["optimize"] and self.stage == "endpoints":
            softexit.trigger("Optimization of endpoints in NEB is not implemented yet.")

        # Endpoints are optimized or optimization is not required
        elif self.stage == "neb" :
            # Fetch spring constants
            self.nebgm.kappa = self.spring["kappa"]

            self.ptime = self.ttime = 0
            self.qtime = -time.time()

            if self.mode in ["damped_bfgs", ]:
                # BFGS-family algorithms
                if step == 0: # need to add a condition when after the endpoints.
                    # Initialize direction to the steepest descent direction
                    info(" @NEB: calling NEBGradientMapper at step 0.", verbosity.debug)
                    self.nebpot, self.nebgrad = self.nebgm(self.beads.q[:, self.nebgm.fixatoms_mask])
                    info(" @NEB: NEBGradientMapper returned nebpot and nebgrad.",
                         verbosity.debug)

                    # Set the initial direction to the direction of NEB forces
#                    self.d = -self.nebgrad
                    self.old_x = self.beads.q.copy()[:, self.nebgm.fixatoms_mask]
#                    self.old_f = self.forces.f.copy()[:, self.nebgm.fixatoms_mask]

                    # Inverse Hessian for BFGS should be
                    # initialized in bind(), but better to check. (remove in cleanup)
                    if self.invhessian.shape != (self.beads.q.size, self.beads.q.size):
                        print((self.invhessian.shape, self.beads.q.size))
                        softexit.trigger("Invhessian not initialized correctly in NEB.")

                # Self instances will be updated in the optimizer, so we store the copies.
                # old_nebpot is used later as a convergence criterion.
                old_nebpot, old_nebgrad = (self.nebpot.copy(), self.nebgrad.copy())

                # Ensure zero displacements on the endpoints
#                self.d[0, :] = self.d[-1, :] = 0.

                # Reducing dimension of the invhessian before passing it to optimizer
                masked_invhessian = self.invhessian[
                    np.ix_(np.tile(self.nebgm.fixatoms_mask, self.beads.nbeads),
                           np.tile(self.nebgm.fixatoms_mask, self.beads.nbeads)
                    )
                ]

                if self.mode == "damped_bfgs":
                    info(" @NEB: before Damped_BFGS() call", verbosity.debug)
                    quality = Damped_BFGS(
                                x0=self.old_x,
                                fdf=self.nebgm,
                                fdf0=(self.nebpot, self.nebgrad),
                                invhessian=masked_invhessian,
                                big_step=self.big_step,
                              )
                    info(" @NEB: after Damped_BFGS() call", verbosity.debug)

                    # Let's try to make big_step brave enough
                    if quality > 0.8:
                        self.big_step *= 1.3
                        info(" @NEB: increasing big_step to %.6f bohr." % self.big_step,
                             verbosity.debug)
                    elif quality > 0.6:
                        self.big_step *= 1.1
                        info(" @NEB: increasing big_step to %.6f bohr." % self.big_step,
                             verbosity.debug)
                    else:
                        self.big_step *= 0.5
                        if self.big_step <= self.tolerances["position"]:
                            self.big_step = 10 * self.tolerances["position"]
                        info(" @NEB: Step direction far from nebgrad, "
                             "reducing big_step to %.6f bohr." % self.big_step,
                             verbosity.debug)

#                elif self.mode == "bfgs":
#                    info(" @NEB: before BFGS call", verbosity.debug)
#                    BFGS(
#                        x0=self.old_x,
#                        d0=self.d,
#                        fdf=self.nebgm,
#                        fdf0=(self.nebpot, self.nebgrad),
#                        invhessian=masked_invhessian,
#                        big_step=self.big_step,
#                        tol=self.ls_options["tolerance"] * self.tolerances["energy"],
#                        itmax=self.ls_options["iter"],
#                    )
#                    info(" @NEB: after BFGS call", verbosity.debug)
#                elif self.mode == "bfgstrm":
#                    softexit.trigger("BFGSTRM doesn't work yet.")
                else:
                    softexit.trigger("NEBMover has a bug in BFGS-family options.")

                # Restore dimensionality of the invhessian
                self.invhessian[
                    np.ix_(np.tile(self.nebgm.fixatoms_mask, self.beads.nbeads),
                           np.tile(self.nebgm.fixatoms_mask, self.beads.nbeads)
                    )
                ] = masked_invhessian

                # Update positions
                self.beads.q[:] = self.nebgm.dbeads.q
                info(" @NEB: bead positions updated.", verbosity.debug)

#                print(" @NEB: Before re-calling NEBGradientMapper")
                self.nebpot, self.nebgrad = self.nebgm(
                    self.beads.q[:, self.nebgm.fixatoms_mask])
#                print(" @NEB: After re-calling NEBGradientMapper")
#                print("Properties after calling a BFGS-family optimizer:")
#                print("old_nebgrad:")
#                print(old_nebgrad)
#                print("self.nebgrad:")
#                print(self.nebgrad)

                # dx = current position - previous position.
                # Use to determine converged minimization
                dx = np.amax(
                        np.abs(self.beads.q[:, self.nebgm.fixatoms_mask] - self.old_x)
                     )

                # Store old positions
                self.old_x[:] = self.beads.q[:, self.nebgm.fixatoms_mask]
                # This transfers forces from the mapper to the "main" beads,
                # so that recalculation won't be triggered after the step.
                self.forces.transfer_forces(self.nebgm.dforces)


            # TODO: Routines for L-BFGS, SD, CG
            else:
                if self.mode in ["sd", "cg", "lbfgs"]:
                    softexit.trigger(
                        "LBFGS, SD and CG are not implemented for NEB. Try damped_bfgs.")

            self.qtime += time.time()


            # Check convergence criteria
            if (
                (np.amax(np.abs(self.nebpot - old_nebpot)) / (self.beads.nbeads * self.beads.natoms)
                    <= self.tolerances["energy"]
                )
                and (np.amax(np.abs(self.nebgrad)) <= self.tolerances["force"])
                and (dx <= self.tolerances["position"])
            ):
                info("NEB path optimization converged.", verbosity.medium)

                # Initialization of climbing
                if self.use_climb:
                    info("Starting climbing.", verbosity.medium)
                    self.stage = "climb"
                    climb_index = np.argmax(self.forces.pots)
                    if climb_index in [0, self.beads.nbeads -1]:
                        softexit.trigger("ERROR: climbing bead is the endpoint.")
                    self.climbgm.q_prev[:] = self.beads.q[climb_index - 1]
                    self.climbgm.q_next[:] = self.beads.q[climb_index + 1]
                    self.invhessian = np.eye(self.beads.q[0].size)

                    info(" @NEB_CLIMB: calling NEBClimbGrMapper first time.",
                         verbosity.debug)
                    self.nebpot, self.nebgrad = self.climbgm(
                        self.beads.q[climb_index, self.climbgm.fixatoms_mask])

                    self.old_x = self.beads.q.copy()[climb_index,
                                                     self.climbgm.fixatoms_mask]
                    self.stage = "climb"

            else:
                info(
                    " @NEB: Not converged, deltaEnergy = %.8f, tol = %.8f per atom"
                    % (np.amax(np.abs(self.nebpot - old_nebpot)) / (self.beads.nbeads * self.beads.natoms),
                       self.tolerances["energy"]),
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

#======================== C L I M B ========================

        # Climbing image optimization
        elif self.stage == "climb":
#            softexit.trigger("Climbing image NEB is not implemented yet.")
            self.ptime = self.ttime = 0
            self.qtime = -time.time()

            if self.mode == "damped_bfgs":
                # BFGS-family algorithms

                # Self instances will be updated in the optimizer, so we store the copies.
                # old_nebpot is used later as a convergence criterion.
                old_nebpot, old_nebgrad = (self.nebpot.copy(), self.nebgrad.copy())

                # Reducing dimension of the invhessian before passing it to optimizer
                masked_invhessian = self.invhessian[
                    np.ix_(self.climbgm.fixatoms_mask, self.climbgm.fixatoms_mask)
                ]

                if self.mode == "damped_bfgs":
                    info(" @NEB_CLIMB: before Damped_BFGS() call", verbosity.debug)
                    quality = Damped_BFGS(
                                x0=self.old_x,
                                fdf=self.climbgm,
                                fdf0=(self.nebpot, self.nebgrad),
                                invhessian=masked_invhessian,
                                big_step=self.big_step,
                              )
                    info(" @NEB_CLIMB: after Damped_BFGS() call", verbosity.debug)

                    # Let's try to make big_step brave enough
                    if quality > 0.8:
                        self.big_step *= 1.3
                        info(" @NEB_CLIMB: increasing big_step to %.6f bohr." % self.big_step,
                             verbosity.debug)
                    elif quality > 0.6:
                        self.big_step *= 1.1
                        info(" @NEB_CLIMB: increasing big_step to %.6f bohr." % self.big_step,
                             verbosity.debug)
                    else:
                        self.big_step *= 0.5
                        if self.big_step <= self.tolerances["position"]:
                            self.big_step = 10 * self.tolerances["position"]
                        info(" @NEB_CLIMB: Step direction far from climbgrad, "
                             "reducing big_step to %.6f bohr." % self.big_step,
                             verbosity.debug)

                # Restore dimensionality of the invhessian
                self.invhessian[
                    np.ix_(self.climbgm.fixatoms_mask, self.climbgm.fixatoms_mask)
                ] = masked_invhessian

                # Update positions
                self.beads.q[climb_index] = self.climbgm.dbeads.q
                info(" @NEB_CLIMB: bead positions updated.", verbosity.debug)

                self.nebpot, self.nebgrad = self.climbgm(
                    self.beads.q[climb_index, self.climbgm.fixatoms_mask])

                # dx = current position - previous position.
                # Use to determine converged minimization
                dx = np.amax(np.abs(
                    self.beads.q[climb_index, self.climbgm.fixatoms_mask] - self.old_x)
                )

                # Store old positions
                self.old_x[:] = self.beads.q[:, self.climbgm.fixatoms_mask]
                # This transfers forces from the mapper to the "main" beads,
                # so that recalculation won't be triggered after the step.
                print("Before transfer_forces_manual in climb.")
                tmp_forces = self.forces.f.copy()
                tmp_forces[climb_index, :] = self.nebgrad[:]
                self.forces.transfer_forces_manual(
                    new_q=self.beads.q,
                    new_v=self.beads.v,
                    new_forces=tmp_forces)
                print("After transfer_forces_manual in climb.")


            # TODO: Routines for L-BFGS, SD, CG, FIRE, ...
            else:
                softexit.trigger(
                    "Try damped_bfgs, other algorithms are not implemented for NEB.")

            self.qtime += time.time()


            # Check convergence criteria
            if (
                (np.amax(np.abs(self.nebpot - old_nebpot)) / self.beads.natoms
                    <= self.tolerances["energy"]
                )
                and (np.amax(np.abs(self.nebgrad)) <= self.tolerances["force"])
                and (dx <= self.tolerances["position"])
            ):
                info("NEB path optimization converged.", verbosity.medium)

            else:
                info(
                    " @NEB_CLIMB: Not converged, deltaEnergy = %.8f, tol = %.8f per atom"
                    % (np.amax(np.abs(self.nebpot - old_nebpot)) / self.beads.natoms,
                       self.tolerances["energy"]),
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

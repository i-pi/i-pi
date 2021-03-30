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
    """Creation of the multi-dimensional function that will be minimized.
       Functional analog of a GradientMapper in geop.py

    Attributes:
#        x0: initial position
#        d: move direction
#        xold: position from previous step
        kappa: spring constants"""

    def __init__(self):
#        self.x0 = None
#        self.d = None
#        self.xold = None
        self.kappa = None

    def bind(self, ens):
        self.dbeads = ens.beads.copy()
        self.dcell = ens.cell.copy()
        self.dforces = ens.forces.copy(self.dbeads, self.dcell)

    def __call__(self, x):
        """ Returns the potential for all beads and the gradient.
        """

        # Bead positions
        self.dbeads.q = x
        bq = self.dbeads.q
#        print("bq [bohr]:")
#        print(bq)

        # Forces
        bf = self.dforces.f.copy()
#        print("bf [atomic]:")
#        print(bf)
        # Zeroing endpoint forces
        bf[0,:] = bf[-1,:] = 0.

        # Bead energies (needed for improved tangents)
        be = self.dforces.pots.copy()

        # Number of images
        nimg = self.dbeads.nbeads

        # Number of atoms
        nat = self.dbeads.natoms

        # Array for spring constants
        kappa = np.zeros(nimg)

        # Get tangents. End images are distinct, fixed, pre-relaxed configurations
        btau = np.zeros((nimg, 3 * nat), float)
        for ii in range(1, nimg - 1):
            d1 = bq[ii] - bq[ii - 1]  # tau minus
            d2 = bq[ii + 1] - bq[ii]  # tau plus

            # Old implementation of NEB tangents
            # btau[ii] = d1 / npnorm(d1) + d2 / npnorm(d2)
            # btau[ii] /= npnorm(btau[ii])

            # Improved tangent estimate
            # J. Chem. Phys. 113, 9978 (2000) https://doi.org/10.1063/1.1323224

#            print("Energies for tangent %i:" % ii)
#            print((be[ii + 1], be[ii], be[ii - 1]))
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

        print("Tangents in NEBGradientMapper.call:")
        print(btau)

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

#        print("NEB forces in NEBGradientMapper.call:")
#        print(bf)

        # Including spring energy into potential (experimental)
        # The reason for having such term is that many optimization algorithms
        # need energy for line search or TRM update.
        for ii in range(1, nimg - 1):
            dl = (npnorm(bq[ii + 1] - bq[ii]) - npnorm(bq[ii] - bq[ii - 1]))
            be[ii] += kappa[ii] * dl * dl

        # Return forces and the total potential energy of the band
        e = be.sum()
        g = -bf
        print("Sum of E_pot of all beads in NEBGradientMapper.call:")
        print(e)
        print("NEB gradient in NEBGradientMapper.call:")
        print(g)
        return e, g


class NEBMover(Motion):
    """Nudged elastic band routine.

    Attributes:
        mode: minimizer to use for NEB
        biggest_step: maximum step size for BFGS/L-BFGS
        old_force: force from previous iteration (not NEB-projected)
        old_direction: direction from previous iteration
        old_potential: bead potentials from previous iteration
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
        spring:
            varsprings: T/F for variable spring constants
            kappa: single spring constant if varsprings is F
            kappamax: max spring constant if varsprings is T *** NOT YET IMPLEMENTED ***
            kappamin: min spring constant if varsprings is T *** NOT YET IMPLEMENTED ***
        climb: flag for climbing image NEB *** NOT YET IMPLEMENTED ***
    """

    def __init__(
        self,
        fixcom=False,
        fixatoms=None,
        # mode="lbfgs",
        mode="bfgs",
        biggest_step=1.0,
        old_coord=np.zeros(0, float),
        old_force=np.zeros(0, float),
        old_potential=np.zeros(0, float),
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
        climb=False,
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
        self.old_x = old_coord
#        self.old_f = old_force
        self.pot = old_potential
        self.nebgrad = old_nebgradient
        self.d = old_direction
        self.invhessian = invhessian_bfgs
        self.corrections = corrections_lbfgs
        self.qlist = qlist_lbfgs
        self.glist = glist_lbfgs
        self.tr_trm = tr_trm
        self.endpoints = endpoints
        self.spring = spring
        self.climb = climb
        self.scale = scale_lbfgs

        self.neblm = NEBLineMapper()
        self.nebgm = NEBGradientMapper()


    def bind(self, ens, beads, nm, cell, bforce, prng, omaker):

        super(NEBMover, self).bind(ens, beads, nm, cell, bforce, prng, omaker)
#        if self.old_f.shape != beads.q.shape:
#            if self.old_f.shape == (0,):
#                self.old_f = np.zeros(beads.q.shape, float)
#            else:
#                raise ValueError("Force size does not match system size.")
        if self.d.shape != beads.q.shape:
            if self.d.shape == (0,):
                self.d = np.zeros(beads.q.shape, float)
            else:
                raise ValueError("Direction size does not match system size.")
        if (self.mode in ["damped_bfgs", "bfgs", "bfgstrm"]
            and self.invhessian.size != (beads.q.size * beads.q.size)
        ):
            if self.invhessian.size == 0:
                self.invhessian = np.eye(
                    beads.q.flatten().size, beads.q.flatten().size, 0, float
                )
            else:
                raise ValueError("Inverse Hessian size does not match system size.")

#        self.neblm.bind(self)
        self.nebgm.bind(self)


    def step(self, step=None):
        """Does one simulation time step."""

        info(" @NEB STEP %d" % step, verbosity.debug)

        # First, optimization of endpoints, if required
        if self.endpoints["optimize"]:
            softexit.trigger("Optimization of endpoints in NEB is not implemented yet.")

        # Endpoints optimized or optimization is not required
        elif not self.climb:
            # Fetch spring constants
            self.nebgm.kappa = self.spring["kappa"]
#            self.neblm.kappa = self.spring["kappa"]

            self.ptime = self.ttime = 0
            self.qtime = -time.time()

            if self.mode in ["damped_bfgs", "bfgs", "bfgstrm"]:
                # BFGS minimization
                if step == 0:
                    # Initialize direction to the steepest descent direction
                    info(" @NEB: calling NEBGradientMapper at step 0.", verbosity.debug)
                    self.pot, self.nebgrad = self.nebgm(self.beads.q)
                    info(" @NEB: NEBGradientMapper returned pot and nebgrad.",
                         verbosity.debug)

                    # Set the initial direction to the direction of NEB forces
                    self.d = -self.nebgrad
                    self.old_x = self.beads.q.copy()
#                    self.old_f = self.forces.f.copy()

                    # Inverse Hessian for BFGS should be
                    # initialized in bind(), but better to check. (remove in cleanup)
                    if self.invhessian.shape != (self.beads.q.size, self.beads.q.size):
                        print((self.invhessian.shape, self.beads.q.size))
                        softexit.trigger("Invhessian not initialized correctly in NEB.")

                # Self instances will be updated in the optimizer, so we store the copies.
                # old_pot is used later as a convergence criterion.
                old_pot, old_nebgrad = (self.pot.copy(), self.nebgrad.copy())

                # Ensure zero displacements of the endpoints
                self.d[0, :] = self.d[-1, :] = 0.

                if self.mode == "damped_bfgs":
                    info(" @NEB: before Damped_BFGS call", verbosity.debug)
                    Damped_BFGS(
                        x0=self.beads.q,
                        fdf=self.nebgm,
                        fdf0=(self.pot, self.nebgrad),
                        invhessian=self.invhessian,
                        big_step=self.big_step,
                    )
                    info(" @NEB: after Damped_BFGS call", verbosity.debug)
                elif self.mode == "bfgs":
                    info(" @NEB: before BFGS call", verbosity.debug)
                    BFGS(
                        x0=self.beads.q,
                        d0=self.d,
                        fdf=self.nebgm,
                        fdf0=(self.pot, self.nebgrad),
                        invhessian=self.invhessian,
                        big_step=self.big_step,
                        tol=self.ls_options["tolerance"] * self.tolerances["energy"],
                        itmax=self.ls_options["iter"],
                    )
                    info(" @NEB: after BFGS call", verbosity.debug)
                elif self.mode == "bfgstrm":
                    softexit.trigger("BFGSTRM doesn't work yet.")
                    info(" @NEB: before BFGSTRM call", verbosity.debug)
                    BFGSTRM(
                        x0=self.beads.q,
                        u0=self.pot,
                        f0=self.nebgrad,
                        h0=self.invhessian,
                        tr=self.tr_trm,
                        mapper=self.nebgm,
                        big_step=self.big_step,
                    )
                    info(" @NEB: after BFGSTRM call", verbosity.debug)
                else:
                    softexit.trigger("NEBMover has a bug in BFGS-family options.")

                # Update positions and forces
                self.beads.q = self.nebgm.dbeads.q
                info(" @NEB: bead positions updated.", verbosity.debug)
                # This transfers forces from the mapper
                # so that recalculation won't be triggered.
                self.forces.transfer_forces(self.nebgm.dforces)

# Temporarily commented, but it is needed for convergence.
# It seems that this block breaks relaxation process...
#                # I call the mapper again and hope that it won't trigger recalculation of forces
#                print(" @NEB: Before re-calling NEBGradientMapper")
#                self.pot, self.nebgrad = self.nebgm(self.beads.q)
#                print(" @NEB: After re-calling NEBGradientMapper")

                print("Properties after (damped) BFGS call:")
                print("old_nebgrad:")
                print(old_nebgrad)
                print("self.nebgrad:")
                print(self.nebgrad)
#                print("self.old_f:")
#                print(self.old_f)
                print("self.forces.f:")
                print(self.forces.f)

                # dx = current position - previous position.
                # Use to determine converged minimization
                dx = np.amax(np.abs(self.beads.q - self.old_x))

                # Store old positions
                self.old_x[:] = self.beads.q


            # TODO: Routines for L-BFGS, SD, CG
            else:
                if self.mode in ["sd", "cg", "lbfgs"]:
                    softexit.trigger(
                        "LBFGS, SD and CG are not implemented for NEB. Try damped_bfgs.")

            self.qtime += time.time()

            # Check convergence criteria
            if (
                (np.amax(np.abs(self.pot - old_pot)) / (self.beads.nbeads * self.beads.natoms)
                    <= self.tolerances["energy"]
                )
                and (np.amax(np.abs(self.nebgrad)) <= self.tolerances["force"])
                and (dx <= self.tolerances["position"])
            ):
                softexit.trigger("Geometry optimization converged. Exiting simulation")
            else:
                info(
                    " @NEB: Not converged, deltaEnergy = %.8f, tol = %.8f per atom"
                    % (np.amax(np.abs(self.pot - old_pot)) / (self.beads.nbeads * self.beads.natoms),
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

        # Climbing image optimization
        elif self.climb:
            softexit.trigger("Climbing image NEB is not implemented yet.")

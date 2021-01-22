"""
Contains classes for instanton  calculations.

Algorithms implemented by Yair Litman and Mariana Rossi, 2017
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np
import time
import sys

from ipi.engine.motion import Motion
from ipi.utils.depend import dstrip, dobject
from ipi.utils.softexit import softexit
from ipi.utils.messages import verbosity, info
from ipi.utils import units
from ipi.utils.mintools import nichols, Powell
from ipi.engine.motion.geop import L_BFGS
from ipi.utils.instools import (
    banded_hessian,
    invmul_banded,
    red2comp,
    get_imvector,
    print_instanton_geo,
)
from ipi.utils.instools import print_instanton_hess, diag_banded, ms_pathway
from ipi.utils.hesstools import get_hessian, clean_hessian, get_dynmat
from ipi.engine.beads import Beads

__all__ = ["InstantonMotion"]


class InstantonMotion(Motion):
    """Instanton motion class.

    Attributes:
        mode: minimization algorithm to use
        biggest_step: max allowed step size
        old_force: force on previous step
        hessian:

        mode= type of instanton calculation
        tolerances:
            energy: change in energy tolerance for ending minimization
            force: force/change in force tolerance foe ending minimization
            position: change in position tolerance for ending minimization}
        biggest_step: The maximum step size during the optimization.
        old_pos: The previous step positions during the optimization.
        old_pot: The previous step potential energy during the optimization
        old_force:  The previous step force during the optimization
        opt: The geometry optimization algorithm to be used
        discretization: Allows for non uniform time discretization
        alt_out: (Alternative output) Prints different formatting of outputs for geometry, hessian and bead potential energies.
        All quantities are also accessible from typical i-pi output infrastructure. Default to 1, which prints
        every step. -1 will suppress the output (except the last one). Any other positive number will set the frequency (in steps) with
        which the quantities are written to file.
        prefix: Prefix of the output files.
        delta: Initial stretch amplitude.
        hessian_init: Boolean which decides whether the initial hessian is going to be computed.
        hessian: Stored  Hessian matrix
        hessian_update: The way to update the hessian after each movement
        hessian_asr: Removes the zero frequency vibrational modes depending on the symmerty of the system.
        glist_lbfgs: List of previous gradients (g_n+1 - g_n) for L-BFGS. Number of entries = corrections_lbfgs
        qlist_lbfgs: List of previous positions (x_n+1 - x_n) for L-BFGS. Number of entries = corrections_lbfgs
        scale_lbfgs: Scale choice for the initial hessian.
        corrections_lbfgs: Number of corrections to be stored for L-BFGS
        ls_options: Options for line search methods.
        hessian_final:  Boolean which decides whether the hessian after the optimization will be computed.
        energy_shift: zero of energy (usually it corresponds to reactant state)
    """

    def __init__(
        self,
        fixcom=False,
        fixatoms=None,
        mode="None",
        tolerances={"energy": 1e-5, "force": 1e-4, "position": 1e-3},
        biggest_step=0.3,
        old_pos=np.zeros(0, float),
        old_pot=np.zeros(0, float),
        old_force=np.zeros(0, float),
        opt="None",
        max_e=0.0,
        max_ms=0.0,
        discretization=np.zeros(0, float),
        alt_out=1,
        prefix="instanton",
        delta=np.zeros(0, float),
        hessian_init=None,
        hessian=np.eye(0, 0, 0, float),
        hessian_update=None,
        hessian_asr=None,
        qlist_lbfgs=np.zeros(0, float),
        glist_lbfgs=np.zeros(0, float),
        scale_lbfgs=1,
        corrections_lbfgs=5,
        ls_options={"tolerance": 1e-1, "iter": 100},
        old_direction=np.zeros(0, float),
        hessian_final="False",
        energy_shift=np.zeros(0, float),
    ):
        """Initialises InstantonMotion."""

        super(InstantonMotion, self).__init__(fixcom=fixcom, fixatoms=fixatoms)

        self.options = {}  # Optimization options
        self.optarrays = {}  # Optimization arrays

        # Optimization mode
        self.options["mode"] = mode

        # Generic optimization
        # self.big_step = biggest_step
        # self.tolerances = tolerances

        self.options["tolerances"] = tolerances
        self.options["save"] = alt_out
        self.options["prefix"] = prefix
        self.options["hessian_final"] = hessian_final

        self.options["max_e"] = max_e
        self.options["max_ms"] = max_ms
        self.options["discretization"] = discretization
        self.optarrays["big_step"] = biggest_step
        self.optarrays["energy_shift"] = energy_shift
        self.optarrays["delta"] = delta
        self.optarrays["old_x"] = old_pos
        self.optarrays["old_u"] = old_pot
        self.optarrays["old_f"] = old_force

        # We set the default optimization algorithm depending on the mode.
        if mode == "rate":
            if opt == "None":
                opt = "nichols"
            self.options["opt"] = opt

        elif mode == "splitting":
            if opt == "None":
                opt = "lbfgs"
            self.options["opt"] = opt

        if (
            self.options["opt"] == "nichols"
            or self.options["opt"] == "NR"
            or self.options["opt"] == "lanczos"
        ):

            self.options["hessian_update"] = hessian_update
            self.options["hessian_asr"] = hessian_asr
            self.options["hessian_init"] = hessian_init
            self.optarrays["hessian"] = hessian

            if self.options["opt"] == "nichols":
                self.optimizer = NicholsOptimizer()
            elif self.options["opt"] == "NR":
                self.optimizer = NROptimizer()
            else:
                self.optimizer = LanczosOptimizer()

        elif self.options["opt"] == "lbfgs":
            self.optimizer = LBFGSOptimizer()
            self.optarrays["hessian"] = hessian  # Only for initial (to spread) or final
            self.options["hessian_asr"] = hessian_asr

            self.options["corrections"] = corrections_lbfgs
            self.options["scale"] = scale_lbfgs
            self.options["ls_options"] = ls_options

            self.optarrays["qlist"] = qlist_lbfgs
            self.optarrays["glist"] = glist_lbfgs
            self.optarrays["d"] = old_direction

        if self.options["opt"] == "NR":
            info(
                "Note that we need scipy to use NR. If storage and diagonalization of the full hessian is not a "
                "problem use nichols even though it may not be as efficient.",
                verbosity.low,
            )

    def bind(self, ens, beads, nm, cell, bforce, prng, omaker):
        """Binds beads, cell, bforce and prng to InstantonMotion

        Args:
        beads: The beads object from which the bead positions are taken.
        nm: A normal modes object used to do the normal modes transformation.
        cell: The cell object from which the system box is taken.
        bforce: The forcefield object from which the force and virial are taken.
        prng: The random number generator object which controls random number generation.
        """

        super(InstantonMotion, self).bind(ens, beads, nm, cell, bforce, prng, omaker)
        # Binds optimizer

        self.optimizer.bind(self)

    def step(self, step=None):
        self.optimizer.step(step)


class Fix(object):
    """Class that applies a fixatoms type constrain"""

    def __init__(self, natoms, fixatoms, nbeads=1):
        self.natoms = natoms
        self.nbeads = nbeads
        self.fixatoms = fixatoms

        self.mask0 = np.delete(np.arange(self.natoms), self.fixatoms)

        mask1 = np.ones(3 * self.natoms, dtype=bool)
        for i in range(3):
            mask1[3 * self.fixatoms + i] = False
        self.mask1 = np.arange(3 * self.natoms)[mask1]

        mask2 = np.tile(mask1, self.nbeads)
        self.mask2 = np.arange(3 * self.natoms * self.nbeads)[mask2]

    def get_mask(self, m):

        if m == 0:
            return self.mask0
        elif m == 1:
            return self.mask1
        elif m == 2:
            return self.mask2
        else:
            raise ValueError("Mask number not valid")

    def get_active_array(self, arrays):
        """Functions that gets the subarray corresponding to the active degrees-of-freedom of the
        full dimensional array"""

        activearrays = {}
        for key in arrays:

            if (
                key == "old_u"
                or key == "big_step"
                or key == "delta"
                or key == "energy_shift"
                or key == "initial_hessian"
            ):
                t = -1
            elif key == "old_x" or key == "old_f" or key == "d":
                t = 1
            elif key == "hessian":
                t = 2
            elif key == "qlist" or key == "glist":
                t = 3
            else:
                raise ValueError(
                    "@get_active_array: There is an array that we can't recognize"
                )

            activearrays[key] = self.get_active_vector(arrays[key], t)

        return activearrays

    def get_full_vector(self, vector, t):
        """Set 0 the degrees of freedom (dof) corresponding to the fix atoms
        IN:
            fixatoms   indexes of the fixed atoms
            vector     vector to be reduced
            t          type of array:
                type=-1 : do nothing
                type=0 : names (natoms )
                type=1 : pos , force or m3 (nbeads,dof)
                type=2 : hessian (dof, nbeads*dof)
                type=3 : qlist or glist (corrections, nbeads*dof)
        OUT:
            clean_vector  reduced vector
        """
        if len(self.fixatoms) == 0 or t == -1:
            return vector

        if t == 1:

            full_vector = np.zeros((self.nbeads, 3 * self.natoms))
            full_vector[:, self.get_mask(1)] = vector

            return full_vector

        elif t == 2:

            full_vector = np.zeros((3 * self.natoms, 3 * self.natoms * self.nbeads))

            ii = 0
            for i in self.get_mask(1):
                full_vector[i, self.get_mask(2)] = vector[ii]
                ii += 1

            return full_vector

        elif t == 3:

            full_vector = np.zeros((vector.shape[0], 3 * self.natoms * self.nbeads))
            full_vector[:, self.fix.get_mask(2)] = vector

            return full_vector

        else:

            raise ValueError("@apply_fix_atoms: type number is not valid")

    def get_active_vector(self, vector, t):
        """Delete the degrees of freedom (dof) corresponding to the fix atoms
        IN:
            fixatoms   indexes of the fixed atoms
            vector     vector to be reduced
            t          type of array:
                type=-1 : do nothing
                type=0 : names (natoms )
                type=1 : pos , force or m3 (nbeads,dof)
                type=2 : hessian (dof, nbeads*dof)
                type=3 : qlist or glist (corrections, nbeads*dof)
        OUT:
            clean_vector  reduced vector
        """
        if len(self.fixatoms) == 0 or t == -1:
            return vector
        if t == 0:
            return vector[self.mask0]
        elif t == 1:
            return vector[:, self.mask1]
        elif t == 2:
            aux = vector[self.mask1]
            return aux[:, self.mask2]
        elif t == 3:
            return vector[:, self.mask2]
        else:
            raise ValueError("@apply_fix_atoms: type number is not valid")


class GradientMapper(object):

    """Creation of the multi-dimensional function to compute the physical potential and forces

    Attributes:
        dbeads:  copy of the bead object
        dcell:   copy of the cell object
        dforces: copy of the forces object
    """

    def __init__(self):
        self.fcount = 0
        pass

    def bind(self, dumop, discretization, max_ms, max_e):

        self.dbeads = dumop.beads.copy()
        self.dcell = dumop.cell.copy()
        self.dforces = dumop.forces.copy(self.dbeads, self.dcell)
        self.fix = Fix(dumop.beads.natoms, dumop.fixatoms, dumop.beads.nbeads)
        self.set_coef(discretization)
        if max_ms > 0 or max_e > 0:
            self.spline = True

            if max_ms > 0:
                self.max_ms = max_ms
            else:
                self.max_ms = 1000000
            if max_e > 0:
                self.max_e = max_e
            else:
                self.max_e = 10000000
        else:
            self.spline = False

    def set_coef(self, coef):
        self.coef = coef.reshape(-1, 1)

    def set_pos(self, x):
        """Set the positions """
        self.dbeads.q = x

    def save(self, e, g):
        self.pot = e
        self.f = -g

    def __call__(self, x, full=False, new_disc=True):
        """computes energy and gradient for optimization step"""
        self.fcount += 1
        full_q = x.copy()
        full_mspath = ms_pathway(full_q, self.dbeads.m3)

        if self.spline:
            try:
                from scipy.interpolate import interp1d
            except ImportError:
                softexit.trigger("Scipy required to use  max_ms >0")

            indexes = list()
            indexes.append(0)
            old_index = 0
            for i in range(1, self.dbeads.nbeads):
                if (full_mspath[i] - full_mspath[old_index] > self.max_ms) or (
                    np.absolute(self.pot[i] - self.pot[old_index]) > self.max_e
                ):
                    indexes.append(i)
                    old_index = i
            if self.dbeads.nbeads - 1 not in indexes:
                indexes.append(self.dbeads.nbeads - 1)
            info(
                "The reduced RP for this step has {} beads.".format(len(indexes)),
                verbosity.low,
            )
            if len(indexes) <= 2:
                softexit.trigger(
                    "Too few beads fulfill criteria. Please reduce max_ms or max_e"
                )
        else:
            indexes = np.arange(self.dbeads.nbeads)

        # Create reduced bead and force object and evaluate forces
        reduced_b = Beads(self.dbeads.natoms, len(indexes))
        reduced_b.q[:] = full_q[indexes]
        reduced_b.m[:] = self.dbeads.m
        reduced_b.names[:] = self.dbeads.names

        reduced_cell = self.dcell.copy()
        reduced_forces = self.dforces.copy(reduced_b, reduced_cell)

        rpots = reduced_forces.pots  # reduced energy
        rforces = reduced_forces.f  # reduced gradient

        # Interpolate if necessary to get full pot and forces
        if self.spline:
            red_mspath = full_mspath[indexes]
            spline = interp1d(red_mspath, rpots.T, kind="cubic")
            full_pot = spline(full_mspath).T
            spline = interp1d(red_mspath, rforces.T, kind="cubic")
            full_forces = spline(full_mspath).T
        else:
            full_pot = rpots
            full_forces = rforces

        # This forces the update of the forces
        self.dbeads.q[:] = x[:]
        self.dforces.transfer_forces_manual([full_q], [full_pot], [full_forces])

        # e = self.dforces.pot   # Energy
        # g = -self.dforces.f    # Gradient
        e = np.sum(full_pot)  # Energy
        g = -full_forces  # Gradient

        self.save(full_pot, g)

        if not full:
            g = self.fix.get_active_vector(g, 1)

        # APPLY OTHERS CONSTRAIN?

        # Discretization
        if new_disc:
            e = e * (self.coef[1:] + self.coef[:-1]) / 2
            g = g * (self.coef[1:] + self.coef[:-1]) / 2

        return e, g


class SpringMapper(object):
    """Creation of the multi-dimensional function to compute full or half ring polymer pot
    and forces.
    """

    def __init__(self):

        self.pot = None
        self.f = None
        self.dbeads = None
        self.fix = None

    def bind(self, dumop, discretization):

        self.temp = dumop.temp
        self.fix = Fix(dumop.beads.natoms, dumop.fixatoms, dumop.beads.nbeads)
        self.dbeads = Beads(
            dumop.beads.natoms - len(dumop.fixatoms), dumop.beads.nbeads
        )
        self.dbeads.q[:] = self.fix.get_active_vector(dumop.beads.copy().q, 1)
        self.dbeads.m[:] = self.fix.get_active_vector(dumop.beads.copy().m, 0)
        self.dbeads.names[:] = self.fix.get_active_vector(dumop.beads.copy().names, 0)
        self.set_coef(discretization)

        if dumop.options["mode"] == "rate":
            self.omega2 = (
                self.temp
                * (2 * self.dbeads.nbeads)
                * units.Constants.kb
                / units.Constants.hbar
            ) ** 2
        elif dumop.options["mode"] == "splitting":
            self.omega2 = (
                self.temp
                * self.dbeads.nbeads
                * units.Constants.kb
                / units.Constants.hbar
            ) ** 2

        if (
            dumop.options["opt"] == "nichols"
            or dumop.options["opt"] == "NR"
            or dumop.options["opt"] == "lanczos"
        ):
            self.h = self.spring_hessian(
                natoms=self.dbeads.natoms,
                nbeads=self.dbeads.nbeads,
                m3=self.dbeads.m3[0],
                omega2=self.omega2,
                coef=self.coef,
            )

    def set_coef(self, coef):
        """ Sets coefficients for non-uniform instanton calculation """
        self.coef = coef.reshape(-1, 1)

    def save(self, e, g):
        """ Stores potential and forces in this class for convenience """
        self.pot = e
        self.f = -g

    @staticmethod
    def spring_hessian(natoms, nbeads, m3, omega2, mode="half", coef=None):
        """Compute the 'spring hessian'

        OUT    h       = hessian with only the spring terms ('spring hessian')
        """
        if coef is None:
            coef = np.ones(nbeads + 1).reshape(-1, 1)

        # Check size of discretization:
        if coef.size != nbeads + 1:
            print("@spring_hessian: discretization size error")
            sys.exit()

        info(" @spring_hessian", verbosity.high)
        ii = natoms * 3
        h = np.zeros([ii * nbeads, ii * nbeads])

        if nbeads == 1:
            return h

        # Diagonal
        h_sp = m3 * omega2
        diag1 = np.diag(h_sp)
        # diag2 = np.diag(2.0 * h_sp)

        if mode == "half":
            i = 0
            h[i * ii : (i + 1) * ii, i * ii : (i + 1) * ii] += diag1 / coef[1]
            i = nbeads - 1
            h[i * ii : (i + 1) * ii, i * ii : (i + 1) * ii] += diag1 / coef[-2]
            for i in range(1, nbeads - 1):
                h[i * ii : (i + 1) * ii, i * ii : (i + 1) * ii] += diag1 * (
                    1.0 / coef[i] + 1.0 / coef[i + 1]
                )
        elif mode == "splitting" or mode == "full":
            for i in range(0, nbeads):
                h[i * ii : (i + 1) * ii, i * ii : (i + 1) * ii] += diag1 * (
                    1.0 / coef[i] + 1.0 / coef[i + 1]
                )
        else:
            raise ValueError("We can't compute the spring hessian.")

        # Non-Diagonal
        ndiag = np.diag(-h_sp)
        # Quasi-band
        for i in range(0, nbeads - 1):
            h[i * ii : (i + 1) * ii, (i + 1) * ii : (i + 2) * ii] += ndiag * (
                1.0 / coef[i + 1]
            )
            h[(i + 1) * ii : (i + 2) * ii, i * ii : (i + 1) * ii] += ndiag * (
                1.0 / coef[i + 1]
            )

        # Corner
        if mode == "full":
            h[0:ii, (nbeads - 1) * ii : (nbeads) * ii] += ndiag / coef[0]
            h[(nbeads - 1) * ii : (nbeads) * ii, 0:ii] += ndiag / coef[0]

        return h

    def __call__(self, x, ret=True, new_disc=True):
        """Computes spring energy and gradient for instanton optimization step"""

        x = self.fix.get_active_vector(x, 1).copy()

        if new_disc:
            coef = self.coef
        elif new_disc == "one":
            coef = np.ones(self.coef.shape)
        else:
            coef = new_disc.reshape(self.coef.shape)

        if x.shape[0] == 1:  # only one bead
            self.dbeads.q = x
            e = 0.0
            g = np.zeros(x.shape[1])
            self.save(e, g)

        else:
            self.dbeads.q = x
            e = 0.00
            g = np.zeros(self.dbeads.q.shape, float)

            # OLD reference
            # for i in range(self.dbeads.nbeads - 1):
            #    dq = self.dbeads.q[i + 1, :] - self.dbeads.q[i, :]
            #    e += self.omega2 * 0.5 * np.dot(self.dbeads.m3[0] * dq, dq)
            # for i in range(0, self.dbeads.nbeads - 1):
            #    g[i, :] += self.dbeads.m3[i, :] * self.omega2 * (self.dbeads.q[i, :] - self.dbeads.q[i + 1, :])
            # for i in range(1, self.dbeads.nbeads):
            #    g[i, :] += self.dbeads.m3[i, :] * self.omega2 * (self.dbeads.q[i, :] - self.dbeads.q[i - 1, :])

            # With new discretization
            for i in range(self.dbeads.nbeads - 1):
                dq = (self.dbeads.q[i + 1, :] - self.dbeads.q[i, :]) / np.sqrt(
                    coef[i + 1]
                )  # coef[0] and coef[-1] do not enter
                e += self.omega2 * 0.5 * np.dot(self.dbeads.m3[0] * dq, dq)
            for i in range(0, self.dbeads.nbeads - 1):
                g[i, :] += (
                    self.dbeads.m3[i, :]
                    * self.omega2
                    * (
                        self.dbeads.q[i, :] / coef[i + 1]
                        - self.dbeads.q[i + 1, :] / coef[i + 1]
                    )
                )
            for i in range(1, self.dbeads.nbeads):
                g[i, :] += (
                    self.dbeads.m3[i, :]
                    * self.omega2
                    * (
                        self.dbeads.q[i, :] / coef[i]
                        - self.dbeads.q[i - 1, :] / coef[i]
                    )
                )

            self.save(e, g)

        if ret:
            return e, g


class FullMapper(object):
    """Creation of the multi-dimensional function to compute the physical and the spring forces."""

    def __init__(self, im, gm, esum=False):

        self.im = im
        self.gm = gm
        self.esum = esum

    def __call__(self, x):

        e1, g1 = self.im(x)
        e2, g2 = self.gm(x)
        e = e1 + e2
        g = np.add(g1, g2)

        if self.esum:
            e = np.sum(e)

        return e, g


class DummyOptimizer(dobject):
    """ Dummy class for all optimization classes """

    def __init__(self):
        """Initialises object for GradientMapper (physical potential, forces and Hessian)
        and SpringMapper (spring potential, forces and Hessian)"""

        self.options = {}  # Optimization options
        self.optarrays = {}  # Optimization arrays

        self.gm = GradientMapper()
        self.im = SpringMapper()
        self.fm = FullMapper(self.im, self.gm)
        self.fix = None

        self.exit = False
        self.init = False

    def bind(self, geop):
        """
        Bind optimization options and call bind function of Mappers (get beads, cell, forces)
        check whether force size, Hessian size from  match system size
        """

        self.beads = geop.beads
        self.cell = geop.cell
        self.forces = geop.forces
        self.fixcom = geop.fixcom
        self.fixatoms = geop.fixatoms
        self.nm = geop.nm
        # self.ensemble = geop.ens
        self.output_maker = geop.output_maker
        # The resize action must be done before the bind

        if geop.optarrays["old_x"].size != self.beads.q.size:
            if geop.optarrays["old_x"].size == 0:
                geop.optarrays["old_x"] = np.zeros(
                    (self.beads.nbeads, 3 * self.beads.natoms), float
                )
            else:
                raise ValueError("Old positions size does not match system size")
        if geop.optarrays["old_u"].size != self.beads.nbeads:
            if geop.optarrays["old_u"].size == 0:
                geop.optarrays["old_u"] = np.zeros(self.beads.nbeads, float)
            else:
                raise ValueError("Old potential energy size does not match system size")
        if geop.optarrays["old_f"].size != self.beads.q.size:
            if geop.optarrays["old_f"].size == 0:
                geop.optarrays["old_f"] = np.zeros(
                    (self.beads.nbeads, 3 * self.beads.natoms), float
                )
            else:
                raise ValueError("Old forces size does not match system size")

        # Temperature
        self.temp = geop.ensemble.temp
        if geop.ensemble.temp == -1.0 or geop.ensemble.temp == 1.0:
            # This is due to a little inconsistency on the default value
            if self.beads.nbeads != 1:
                raise ValueError(
                    "Temperature must be specified for an Instanton calculation "
                )

        # Optimization mode
        self.options["mode"] = geop.options["mode"]

        # Generic optimization
        if geop.options["discretization"].size != self.beads.nbeads + 1:
            if geop.options["discretization"].size == 0:
                geop.options["discretization"] = np.ones(self.beads.nbeads + 1, float)
            else:
                raise ValueError("Discretization coefficients do not match system size")

        self.options["max_ms"] = geop.options["max_ms"]
        self.options["max_e"] = geop.options["max_e"]
        self.options["discretization"] = geop.options["discretization"]
        self.options["tolerances"] = geop.options["tolerances"]
        self.optarrays["big_step"] = geop.optarrays["big_step"]
        self.optarrays["old_x"] = geop.optarrays["old_x"]
        self.optarrays["old_u"] = geop.optarrays["old_u"]
        self.optarrays["old_f"] = geop.optarrays["old_f"]
        self.options["opt"] = geop.options["opt"]  # optimization algorithm

        # Generic instanton
        self.options["save"] = geop.options["save"]
        self.options["prefix"] = geop.options["prefix"]
        self.optarrays["delta"] = geop.optarrays["delta"]
        self.options["hessian_final"] = geop.options["hessian_final"]
        self.optarrays["energy_shift"] = geop.optarrays["energy_shift"]

        self.gm.bind(
            self,
            self.options["discretization"],
            self.options["max_ms"],
            self.options["max_e"],
        )
        self.im.bind(self, self.options["discretization"])
        self.fix = Fix(geop.beads.natoms, geop.fixatoms, geop.beads.nbeads)

    def initial_geo(self):
        # TODO : add linear interpolation

        info(
            " @GEOP: We stretch the initial geometry with an 'amplitude' of {:4.2f}".format(
                self.optarrays["delta"]
            ),
            verbosity.low,
        )

        fix_onebead = Fix(self.beads.natoms, self.fixatoms, 1)
        active_hessian = fix_onebead.get_active_vector(
            self.optarrays["initial_hessian"], 2
        )
        active_imvector = get_imvector(active_hessian, self.im.dbeads.m3[0].flatten())
        imvector = fix_onebead.get_full_vector(active_imvector, 1).flatten()

        for i in range(self.beads.nbeads):
            self.beads.q[i, :] += (
                self.optarrays["delta"]
                * np.cos(i * np.pi / float(self.beads.nbeads - 1))
                * imvector[:]
            )

    def exitstep(self, d_x_max, step):
        """ Exits the simulation step. Computes time, checks for convergence. """
        self.qtime += time.time()

        tolerances = self.options["tolerances"]
        d_u = self.forces.pot - self.optarrays["old_u"].sum()
        active_force = self.fix.get_active_vector(self.forces.f, 1) + self.im.f

        fff = (
            self.fix.get_active_vector(self.forces.f, 1)
            * (self.im.coef[1:] + self.im.coef[:-1])
            / 2
        )
        active_force = fff + self.im.f

        info(
            " @Exit step: Energy difference: {:4.2e}, (condition: {:4.2e})".format(
                np.absolute(d_u / self.im.dbeads.natoms), tolerances["energy"]
            ),
            verbosity.low,
        )
        info(
            " @Exit step: Maximum force component: {:4.2e}, (condition: {:4.2e})".format(
                np.amax(np.absolute(active_force)), tolerances["force"]
            ),
            verbosity.low,
        )
        info(
            " @Exit step: Maximum component step component: {:4.2e}, (condition: {:4.2e})".format(
                d_x_max, tolerances["position"]
            ),
            verbosity.low,
        )

        if (
            (np.absolute(d_u / self.im.dbeads.natoms) <= tolerances["energy"])
            and (
                (np.amax(np.absolute(active_force)) <= tolerances["force"])
                or (
                    np.linalg.norm(
                        self.forces.f.flatten() - self.optarrays["old_f"].flatten()
                    )
                    <= 1e-08
                )
            )
            and (d_x_max <= tolerances["position"])
        ):

            print_instanton_geo(
                self.options["prefix"] + "_FINAL",
                step,
                self.beads.nbeads,
                self.beads.natoms,
                self.beads.names,
                self.beads.q,
                self.forces.f,
                self.forces.pots,
                self.cell,
                self.optarrays["energy_shift"],
                self.output_maker,
            )
            if self.options["hessian_final"] != "true":
                info("We are not going to compute the final hessian.", verbosity.low)
                info(
                    "Warning, The current hessian is not the real hessian is only an approximation .",
                    verbosity.low,
                )

            else:
                info("We are going to compute the final hessian", verbosity.low)
                active_hessian = get_hessian(
                    self.gm,
                    self.beads.q.copy(),
                    self.beads.natoms,
                    self.beads.nbeads,
                    self.fixatoms,
                )
                self.optarrays["hessian"][:] = self.fix.get_full_vector(
                    active_hessian, 2
                )
                print_instanton_hess(
                    self.options["prefix"] + "_FINAL",
                    step,
                    self.optarrays["hessian"],
                    self.output_maker,
                )

            return True
            # If we just exit here, the last step (including the last hessian) will not be in the RESTART file

        return False

    def update_pos_for(self):
        """ Update positions and forces """

        self.beads.q[:] = self.gm.dbeads.q[:]

        # This forces the update of the forces
        self.forces.transfer_forces(self.gm.dforces)

    def update_old_pos_for(self):
        # Update "old" positions and forces
        self.optarrays["old_x"][:] = self.beads.q
        self.optarrays["old_u"][:] = self.forces.pots
        self.optarrays["old_f"][:] = self.forces.f

    def print_geo(self, step):
        # Print current instanton geometry
        if (
            self.options["save"] > 0 and np.mod(step, self.options["save"]) == 0
        ) or self.exit:
            print_instanton_geo(
                self.options["prefix"],
                step,
                self.beads.nbeads,
                self.beads.natoms,
                self.beads.names,
                self.beads.q,
                self.forces.f,
                self.forces.pots,
                self.cell,
                self.optarrays["energy_shift"],
                self.output_maker,
            )

    def pre_step(self, step=None, adaptative=False):
        """ General tasks that have to be performed before actual step"""

        if self.exit:
            softexit.trigger("Geometry optimization converged. Exiting simulation")

        if not self.init:
            self.initialize(step)

        if adaptative:
            softexit.trigger("Adaptative discretization is not fully implemented")
            # new_coef = <implement_here>
            # self.im.set_coef(coef)
            # self.gm.set_coef(coef)
            raise NotImplementedError

        self.qtime = -time.time()
        info("\n Instanton optimization STEP {}".format(step), verbosity.low)

        activearrays = self.fix.get_active_array(self.optarrays)

        return activearrays

    def step(self, step=None):
        """Dummy simulation time step which does nothing."""
        pass

    def opt_coef(self, coef):
        # func = lambda x: 2 * np.sum(x) - x[0] - x[-1]
        def func(x):
            return 2 * np.sum(x) - x[0] - x[-1]

        coef = np.absolute(coef)
        s = func(coef)
        coef *= 2 * self.im.dbeads.nbeads / s
        # c0   = 2*self.im.dbeads.nbeads - 2*np.sum(coef)
        # coef = np.insert(coef,0,c0)

        self.im.set_coef(coef)

        fphys = self.gm.dforces.f * ((coef[1:] + coef[:-1]) / 2).reshape(-1, 1)
        e, gspring = self.im(self.im.dbeads.q)
        return np.amax(np.absolute(-gspring + fphys))


class HessianOptimizer(DummyOptimizer):
    """ Instanton Rate calculation"""

    def bind(self, geop):
        # call bind function from DummyOptimizer
        super(HessianOptimizer, self).bind(geop)

        self.options["hessian_update"] = geop.options["hessian_update"]
        self.options["hessian_asr"] = geop.options["hessian_asr"]

        if len(self.fixatoms) > 0:
            info(" 'fixatoms' is enabled. Setting asr to None", verbosity.low)
            self.options["hessian_asr"] = "none"
        #        self.output_maker = geop.output_maker
        self.options["hessian_init"] = geop.options["hessian_init"]
        self.optarrays["initial_hessian"] = None

        if geop.optarrays["hessian"].size != (
            self.beads.natoms * 3 * self.beads.q.size
        ):
            if geop.optarrays["hessian"].size == (self.beads.natoms * 3) ** 2:
                self.optarrays["initial_hessian"] = geop.optarrays["hessian"].copy()
                geop.optarrays["hessian"] = np.zeros(
                    (self.beads.natoms * 3, self.beads.q.size), float
                )

            elif (
                geop.optarrays["hessian"].size == 0
                and geop.options["hessian_init"] == "true"
            ):
                info(
                    " Initial hessian is not provided. We are going to compute it.",
                    verbosity.low,
                )
                geop.optarrays["hessian"] = np.zeros(
                    (self.beads.natoms * 3, self.beads.q.size)
                )

                if (
                    (self.beads.q - self.beads.q[0]) == 0
                ).all() and self.beads.nbeads > 1:
                    raise ValueError(
                        """We need an initial hessian in order to create our initial
                    instanton geometry. Please provide a (1-bead) hessian or an initial instanton geometry."""
                    )

            else:
                raise ValueError(
                    " 'Hessian_init' is false, an initial hessian (of the proper size) must be provided."
                )

        self.optarrays["hessian"] = geop.optarrays["hessian"]

    def initialize(self, step):

        if step == 0:

            info(" @GEOP: Initializing INSTANTON", verbosity.low)

            if self.beads.nbeads == 1:

                info(" @GEOP: Classical TS search", verbosity.low)

            else:
                # If the coordinates in all the imaginary time slices are the same
                if ((self.beads.q - self.beads.q[0]) == 0).all():

                    self.initial_geo()
                    self.options["hessian_init"] = "true"

                else:

                    info(
                        " @GEOP: Starting from the provided geometry in the extended phase space",
                        verbosity.low,
                    )
                    if not (
                        (self.optarrays["initial_hessian"] is None)
                        or (
                            self.optarrays["initial_hessian"].size
                            == (self.beads.natoms * 3) ** 2
                        )
                    ):
                        raise ValueError(
                            " You have to provide a hessian with size (3 x natoms)^2 but also geometry in"
                            " the extended phase space (nbeads>1). Please check the inputs\n"
                        )

        self.gm.save(self.forces.pots, self.forces.f)

        if self.options["hessian_init"] == "true":
            active_hessian = get_hessian(
                self.gm,
                self.beads.q.copy(),
                self.beads.natoms,
                self.beads.nbeads,
                self.fixatoms,
            )
            self.optarrays["hessian"][:] = self.fix.get_full_vector(active_hessian, 2)

        if self.im.f is None:
            self.im(self.beads.q, ret=False)  # Init instanton mapper

        self.gm.save(self.forces.pots, self.forces.f)
        self.update_old_pos_for()

        self.init = True

    def update_hessian(self, update, active_hessian, new_x, d_x, d_g):
        """ Update hessian """

        if update == "powell":

            i = self.im.dbeads.natoms * 3
            for j in range(self.im.dbeads.nbeads):
                aux = active_hessian[:, j * i : (j + 1) * i]
                dg = d_g[j, :]
                dx = d_x[j, :]
                Powell(dx, dg, aux)

        elif update == "recompute":
            active_hessian = get_hessian(
                self.gm, new_x, self.beads.natoms, self.beads.nbeads, self.fixatoms
            )

        self.optarrays["hessian"][:] = self.fix.get_full_vector(active_hessian, 2)

    def print_hess(self, step):
        if (
            self.options["save"] > 0 and np.mod(step, self.options["save"]) == 0
        ) or self.exit:
            print_instanton_hess(
                self.options["prefix"],
                step,
                self.optarrays["hessian"],
                self.output_maker,
            )

    def post_step(self, step, new_x, d_x, activearrays):
        """ General tasks that have to be performed after the  actual step"""

        d_x_max = np.amax(np.absolute(d_x))
        info("Current step norm = {}".format(d_x_max), verbosity.medium)

        # Get new energy (u)  and forces(f) using mapper
        self.im(new_x, ret=False, new_disc=True)  # Only to update the mapper

        u, g2 = self.gm(new_x, new_disc=False)
        f = -g2
        d_g = np.subtract(activearrays["old_f"], f)

        # Update
        self.update_hessian(
            self.options["hessian_update"], activearrays["hessian"], new_x, d_x, d_g
        )
        self.update_pos_for()

        #  Print
        self.print_geo(step)
        self.print_hess(step)

        # Check Exit and only then update old arrays
        self.exit = self.exitstep(d_x_max, step)
        self.update_old_pos_for()


class NicholsOptimizer(HessianOptimizer):
    """ Class that implements a nichols optimizations. It can find first order saddle points or minimum"""

    def bind(self, geop):
        # call bind function from HessianOptimizer
        super(NicholsOptimizer, self).bind(geop)

    def initialize(self, step):
        # call initialize function from HessianOptimizer
        super(NicholsOptimizer, self).initialize(step)

    def step(self, step=None):
        """ Does one simulation time step."""

        activearrays = self.pre_step(step)

        # First construct complete hessian from reduced
        h0 = red2comp(
            activearrays["hessian"],
            self.im.dbeads.nbeads,
            self.im.dbeads.natoms,
            self.im.coef,
        )

        # Add spring terms to the physical hessian
        h1 = np.add(self.im.h, h0)

        # Get eigenvalues and eigenvector.
        d, w = clean_hessian(
            h1,
            self.im.dbeads.q,
            self.im.dbeads.natoms,
            self.im.dbeads.nbeads,
            self.im.dbeads.m,
            self.im.dbeads.m3,
            self.options["hessian_asr"],
        )

        # d,w =np.linalg.eigh(h1) #Cartesian
        info(
            "\n@Nichols: 1st freq {} cm^-1".format(
                units.unit_to_user(
                    "frequency", "inversecm", np.sign(d[0]) * np.sqrt(np.absolute(d[0]))
                )
            ),
            verbosity.medium,
        )
        info(
            "@Nichols: 2nd freq {} cm^-1".format(
                units.unit_to_user(
                    "frequency", "inversecm", np.sign(d[1]) * np.sqrt(np.absolute(d[1]))
                )
            ),
            verbosity.medium,
        )
        info(
            "@Nichols: 3rd freq {} cm^-1".format(
                units.unit_to_user(
                    "frequency", "inversecm", np.sign(d[2]) * np.sqrt(np.absolute(d[2]))
                )
            ),
            verbosity.medium,
        )
        # info('@Nichols: 4th freq {} cm^-1'.format(units.unit_to_user('frequency','inversecm',np.sign(d[3])*np.sqrt(np.absolute(d[3])))),verbosity.medium)
        # info('@Nichols: 8th freq {} cm^-1\n'.format(units.unit_to_user('frequency','inversecm',np.sign(d[7])*np.sqrt(np.absolute(d[7])))),verbosity.medium)

        # Find new movement direction
        if self.options["mode"] == "rate":
            f = activearrays["old_f"] * (self.im.coef[1:] + self.im.coef[:-1]) / 2
            d_x = nichols(
                f, self.im.f, d, w, self.im.dbeads.m3, activearrays["big_step"]
            )
        elif self.options["mode"] == "splitting":
            d_x = nichols(
                activearrays["old_f"],
                self.im.f,
                d,
                w,
                self.im.dbeads.m3,
                activearrays["big_step"],
                mode=0,
            )

        # Rescale step if necessary
        if np.amax(np.absolute(d_x)) > activearrays["big_step"]:
            info(
                "Step norm, scaled down to {}".format(activearrays["big_step"]),
                verbosity.low,
            )
            d_x *= activearrays["big_step"] / np.amax(np.absolute(d_x))

        # Get the new full-position
        d_x_full = self.fix.get_full_vector(d_x, t=1)
        new_x = self.optarrays["old_x"].copy() + d_x_full

        self.post_step(step, new_x, d_x, activearrays)


class NROptimizer(HessianOptimizer):
    """ Class that implements a Newton-Raphson optimizations. It can find first order saddle points or minima"""

    def bind(self, geop):
        # call bind function from HessianOptimizer
        super(NROptimizer, self).bind(geop)

    def initialize(self, step):
        # call initialize function from HessianOptimizer
        super(NROptimizer, self).initialize(step)

    def step(self, step=None):
        """ Does one simulation time step."""
        activearrays = self.pre_step(step)

        dyn_mat = get_dynmat(
            activearrays["hessian"], self.im.dbeads.m3, self.im.dbeads.nbeads
        )
        h_up_band = banded_hessian(
            dyn_mat, self.im, masses=False, shift=0.0000001
        )  # create upper band matrix

        fff = activearrays["old_f"] * (self.im.coef[1:] + self.im.coef[:-1]) / 2
        f = (fff + self.im.f).reshape(
            self.im.dbeads.natoms * 3 * self.im.dbeads.nbeads, 1
        )
        f = np.multiply(f, self.im.dbeads.m3.reshape(f.shape) ** -0.5)

        d_x = invmul_banded(h_up_band, f).reshape(self.im.dbeads.q.shape)
        d_x = np.multiply(d_x, self.im.dbeads.m3 ** -0.5)

        # Rescale step if necessary
        if np.amax(np.absolute(d_x)) > activearrays["big_step"]:
            info(
                "Step norm, scaled down to {}".format(activearrays["big_step"]),
                verbosity.low,
            )
            d_x *= activearrays["big_step"] / np.amax(np.absolute(d_x))

        # Get the new full-position
        d_x_full = self.fix.get_full_vector(d_x, t=1)
        new_x = self.optarrays["old_x"].copy() + d_x_full

        self.post_step(step, new_x, d_x, activearrays)


class LanczosOptimizer(HessianOptimizer):
    """Class that implements a modified Nichols algorithm based on Lanczos diagonalization to avoid constructing and diagonalizing
    the full (3*natoms*nbeads)^2 matrix"""

    def bind(self, geop):
        # call bind function from HessianOptimizer
        super(LanczosOptimizer, self).bind(geop)

    def initialize(self, step):
        # call initialize function from HessianOptimizer
        super(LanczosOptimizer, self).initialize(step)

    def step(self, step=None):
        """ Does one simulation time step."""

        activearrays = self.pre_step(step)

        fff = activearrays["old_f"] * (self.im.coef[1:] + self.im.coef[:-1]) / 2
        f = (fff + self.im.f).reshape(
            self.im.dbeads.natoms * 3 * self.im.dbeads.nbeads, 1
        )

        banded = False
        banded = True
        if banded:
            # BANDED Version
            # MASS-scaled
            dyn_mat = get_dynmat(
                activearrays["hessian"], self.im.dbeads.m3, self.im.dbeads.nbeads
            )
            h_up_band = banded_hessian(
                dyn_mat, self.im, masses=False, shift=0.000000001
            )  # create upper band matrix
            f = np.multiply(f, self.im.dbeads.m3.reshape(f.shape) ** -0.5)
            # CARTESIAN
            # h_up_band = banded_hessian(activearrays["hessian"], self.im,masses=True)  # create upper band matrix

            d = diag_banded(h_up_band)
        else:
            # FULL dimensions version
            h_0 = red2comp(
                activearrays["hessian"],
                self.im.dbeads.nbeads,
                self.im.dbeads.natoms,
                self.im.coef,
            )
            h_test = np.add(self.im.h, h_0)  # add spring terms to the physical hessian
            d, w = clean_hessian(
                h_test,
                self.im.dbeads.q,
                self.im.dbeads.natoms,
                self.im.dbeads.nbeads,
                self.im.dbeads.m,
                self.im.dbeads.m3,
                None,
            )
            # CARTESIAN
            # d,w =np.linalg.eigh(h_test) #Cartesian
        info(
            "\n@Lanczos: 1st freq {} cm^-1".format(
                units.unit_to_user(
                    "frequency", "inversecm", np.sign(d[0]) * np.sqrt(np.absolute(d[0]))
                )
            ),
            verbosity.medium,
        )
        info(
            "@Lanczos: 2nd freq {} cm^-1".format(
                units.unit_to_user(
                    "frequency", "inversecm", np.sign(d[1]) * np.sqrt(np.absolute(d[1]))
                )
            ),
            verbosity.medium,
        )
        info(
            "@Lanczos: 3rd freq {} cm^-1\n".format(
                units.unit_to_user(
                    "frequency", "inversecm", np.sign(d[2]) * np.sqrt(np.absolute(d[2]))
                )
            ),
            verbosity.medium,
        )

        if d[0] > 0:
            if d[1] / 2 > d[0]:
                alpha = 1
                lamb = (2 * d[0] + d[1]) / 4
            else:
                alpha = (d[1] - d[0]) / d[1]
                lamb = (
                    3 * d[0] + d[1]
                ) / 4  # midpoint between b[0] and b[1]*(1-alpha/2)
        elif d[1] < 0:  # Jeremy Richardson
            if d[1] >= d[0] / 2:
                alpha = 1
                lamb = (d[0] + 2 * d[1]) / 4
            else:
                alpha = (d[0] - d[1]) / d[1]
                lamb = (d[0] + 3 * d[1]) / 4
        # elif d[1] < 0:  #Litman for Second Order Saddle point
        #    alpha = 1
        #    lamb = (d[1] + d[2]) / 4
        #    print 'WARNING: We are not using the standard Nichols'
        #    print 'd_x', d_x[0],d_x[1]

        else:  # Only d[0] <0
            alpha = 1
            lamb = (d[0] + d[1]) / 4

        if banded:
            h_up_band[-1, :] += -np.ones(h_up_band.shape[1]) * lamb
            d_x = invmul_banded(h_up_band, f)
        else:
            h_test = alpha * (h_test - np.eye(h_test.shape[0]) * lamb)
            d_x = np.linalg.solve(h_test, f)

        d_x.shape = self.im.dbeads.q.shape

        # MASS-scaled
        d_x = np.multiply(d_x, self.im.dbeads.m3 ** -0.5)

        # Rescale step if necessary
        if np.amax(np.absolute(d_x)) > activearrays["big_step"]:
            info(
                "Step norm, scaled down to {}".format(activearrays["big_step"]),
                verbosity.low,
            )
            d_x *= activearrays["big_step"] / np.amax(np.absolute(d_x))

        # Get the new full-position
        d_x_full = self.fix.get_full_vector(d_x, t=1)
        new_x = self.optarrays["old_x"].copy() + d_x_full

        self.post_step(step, new_x, d_x, activearrays)


class LBFGSOptimizer(DummyOptimizer):
    def bind(self, geop):
        # call bind function from DummyOptimizer
        super(LBFGSOptimizer, self).bind(geop)

        if geop.optarrays["hessian"].size == (self.beads.natoms * 3) ** 2:
            self.optarrays["initial_hessian"] = geop.optarrays["hessian"].copy()
            geop.optarrays["hessian"] = np.zeros(
                (self.beads.natoms * 3, self.beads.q.size)
            )

        if geop.options["hessian_final"] == "true":
            self.options["hessian_asr"] = geop.options["hessian_asr"]
            if geop.optarrays["hessian"].size == 0:
                geop.optarrays["hessian"] = np.zeros(
                    (self.beads.natoms * 3, self.beads.q.size)
                )
            self.optarrays["hessian"] = geop.optarrays["hessian"]

        self.im.bind(self, self.options["discretization"])

        # Specific for LBFGS
        self.options["corrections"] = geop.options["corrections"]
        self.options["ls_options"] = geop.options["ls_options"]
        if geop.optarrays["qlist"].size != (
            self.options["corrections"] * self.beads.q.size
        ):
            if geop.optarrays["qlist"].size == 0:
                geop.optarrays["qlist"] = np.zeros(
                    (self.options["corrections"], self.beads.q.size), float
                )
            else:
                raise ValueError("qlist size does not match system size")
        if geop.optarrays["glist"].size != (
            self.options["corrections"] * self.beads.q.size
        ):
            if geop.optarrays["glist"].size == 0:
                geop.optarrays["glist"] = np.zeros(
                    (self.options["corrections"], self.beads.q.size), float
                )
            else:
                raise ValueError("qlist size does not match system size")

        self.optarrays["qlist"] = geop.optarrays["qlist"]
        self.optarrays["glist"] = geop.optarrays["glist"]

        if geop.options["scale"] not in [0, 1, 2]:
            raise ValueError("Scale option is not valid")

        self.options["scale"] = geop.options["scale"]

        if geop.optarrays["d"].size != self.beads.q.size:
            if geop.optarrays["d"].size == 0:
                geop.optarrays["d"] = np.zeros(
                    (self.beads.nbeads, 3 * self.beads.natoms), float
                )
            else:
                raise ValueError("Initial direction size does not match system size")

        self.optarrays["d"] = geop.optarrays["d"]

        self.fm.esum = True

    def initialize(self, step):

        if step == 0:
            info(" @GEOP: Initializing instanton", verbosity.low)

            if self.beads.nbeads == 1:
                raise ValueError(
                    "We can not perform an splitting calculation with nbeads =1"
                )

            else:
                if ((self.beads.q - self.beads.q[0]) == 0).all():
                    # If the coordinates in all the imaginary time slices are the same
                    self.initial_geo()
                else:
                    info(
                        " @GEOP: Starting from the provided geometry in the extended phase space",
                        verbosity.low,
                    )

        # This must be done after the stretching and before the self.d.
        if self.im.f is None:
            self.im(self.beads.q, ret=False)  # Init instanton mapper

        if (
            self.optarrays["old_x"]
            == np.zeros((self.beads.nbeads, 3 * self.beads.natoms), float)
        ).all():
            self.optarrays["old_x"][:] = self.beads.q

        # Specific for LBFGS
        if np.linalg.norm(self.optarrays["d"]) == 0.0:
            f = self.forces.f + self.im.f
            self.optarrays["d"] += dstrip(f) / np.sqrt(np.dot(f.flatten(), f.flatten()))

        self.update_old_pos_for()
        self.init = True

    def post_step(self, step, activearrays):

        """ General tasks that have to be performed after the  actual step"""

        # Update
        self.optarrays["qlist"][:] = self.fix.get_full_vector(
            activearrays["qlist"], t=3
        )
        self.optarrays["glist"][:] = self.fix.get_full_vector(
            activearrays["glist"], t=3
        )
        self.optarrays["d"][:] = self.fix.get_full_vector(activearrays["d"], t=1)

        self.update_pos_for()

        self.print_geo(step)

        # Check Exit and only then update old arrays
        d_x_max = np.amax(
            np.absolute(np.subtract(self.beads.q, self.optarrays["old_x"]))
        )
        self.exit = self.exitstep(d_x_max, step)
        self.update_old_pos_for()

    def step(self, step=None):
        """ Does one simulation time step."""

        activearrays = self.pre_step(step)

        e, g = self.fm(self.beads.q)
        fdf0 = (e, g)

        # Do one step. Update the position and force inside the mapper.
        L_BFGS(
            activearrays["old_x"],
            activearrays["d"],
            self.fm,
            activearrays["qlist"],
            activearrays["glist"],
            fdf0,
            activearrays["big_step"],
            self.options["ls_options"]["tolerance"]
            * self.options["tolerances"]["energy"],
            self.options["ls_options"]["iter"],
            self.options["corrections"],
            self.options["scale"],
            step,
        )

        self.post_step(step, activearrays)

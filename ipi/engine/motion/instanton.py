"""
Contains classes for instanton  calculations.

Algorithms implemented by Yair Litman and Mariana Rossi, 2017
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np
import time

from ipi.engine.motion import Motion
from ipi.utils.depend import dstrip, dobject
from ipi.utils.softexit import softexit
from ipi.utils.messages import verbosity, info
from ipi.utils import units
from ipi.utils.mintools import nichols, Powell
from ipi.engine.motion.geop import L_BFGS
from ipi.utils.instools import banded_hessian, invmul_banded, red2comp, get_hessian, clean_hessian, get_imvector, print_instanton_geo, print_instanton_hess

__all__ = ['InstantonMotion']


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
        alt_out: (Alternative output) Prints different formatting of outputs for geometry, hessian and bead potential
         energies.
        All quantities are also accessible from typical i-pi output infrastructure. Default to 1, which prints 
        every step. -1 will suppress the output (except the last one). Any other positive number will set the frequency 
        (in steps) with which the quantities are written to file.
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

    def __init__(self, fixcom=False, fixatoms=None,
                 mode='None',
                 tolerances={"energy": 1e-5, "force": 1e-4, "position": 1e-3},
                 biggest_step=0.3,
                 old_pos=np.zeros(0, float),
                 old_pot=np.zeros(0, float),
                 old_force=np.zeros(0, float),
                 opt='None',
                 alt_out=1,
                 prefix="INSTANTON",
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
                 hessian_final='False',
                 energy_shift=np.zeros(0, float)):
        """Initialises InstantonMotion.
        """

        super(InstantonMotion, self).__init__(fixcom=fixcom, fixatoms=fixatoms)

        self.options = {}  #Optimization options
        self.optarrays  = {} #Optimization arrays

        # Optimization mode
        self.options["mode"] = mode

        # Generic optimization

        self.options["tolerances"] = tolerances
        self.options["save"] = alt_out
        self.options["prefix"] = prefix
        self.options["hessian_final"] = hessian_final

        self.optarrays["big_step"] = biggest_step
        self.optarrays["energy_shift"] = energy_shift
        self.optarrays["delta"] = delta
        self.optarrays["old_x"] = old_pos
        self.optarrays["old_u"] = old_pot
        self.optarrays["old_f"] = old_force


        # We set the default optimization algorithm depending on the mode.
        if mode == 'rate':
            if opt == 'None':
                opt = 'nichols'
            self.options["opt"] = opt

        elif mode == 'splitting':
            if opt == 'None':
                opt = 'lbfgs'
            self.options["opt"] = opt

        if self.options["opt"] == 'nichols':
            self.optimizer = NicholsOptimizer()
            self.options["hessian_update"] = hessian_update
            self.options["hessian_asr"] = hessian_asr
            self.options["hessian_init"] = hessian_init
            self.optarrays["hessian"] = hessian


        elif self.options["opt"] == 'NR':
            self.optimizer = NicholsOptimizer()
            self.options["hessian_update"] = hessian_update
            self.options["hessian_asr"] = hessian_asr
            self.options["hessian_init"] = hessian_init
            self.optarrays["hessian"] = hessian

        elif self.opt == 'lbfgs':
            self.optimizer = LBFGSOptimizer()
            self.optarrays["hessian"] = hessian  # Only for initial (to spread) or final
            self.options["hessian_asr"] = hessian_asr

            self.options["corrections"] = corrections_lbfgs
            self.options["scale"] = scale_lbfgs
            self.options["ls_options"] = ls_options

            self.optarrays["qlist"] = qlist_lbfgs
            self.optarrays["glist"] = glist_lbfgs
            self.optarrays["d"] = old_direction

        if self.options["opt"] == 'NR':
            info("Note that we need scipy to use NR. If storage and diagonalization of the full hessian is not a "
                 "problem use nichols even though it may not be as efficient.", verbosity.low)

    def bind(self, ens, beads, nm, cell, bforce, prng):
        """Binds beads, cell, bforce and prng to InstantonMotion

            Args:
            beads: The beads object from whcih the bead positions are taken.
            nm: A normal modes object used to do the normal modes transformation.
            cell: The cell object from which the system box is taken.
            bforce: The forcefield object from which the force and virial are taken.
            prng: The random number generator object which controls random number generation.
        """

        super(InstantonMotion, self).bind(ens, beads, nm, cell, bforce, prng)
        # Binds optimizer

        self.optimizer.bind(self)

    def step(self, step=None):
        self.optimizer.step(step)


class GradientMapper(object):

    """Creation of the multi-dimensional function to compute the physical potential and forces

    Attributes:
        dbeads:  copy of the bead object
        dcell:   copy of the cell object
        dforces: copy of the forces object
    """

    def __init__(self):
        pass

    def bind(self, dumop):
        self.dbeads = dumop.beads.copy()
        self.dcell = dumop.cell.copy()
        self.dforces = dumop.forces.copy(self.dbeads, self.dcell)

    def set_pos(self, x):
        """Set the positions """
        self.dbeads.q = x

    def __call__(self, x):
        """computes energy and gradient for optimization step"""

        self.dbeads.q = x
        e = self.dforces.pot   # Energy
        g = -self.dforces.f   # Gradient

        return e, g


class SpringMapper(object):
    """Creation of the multi-dimensional function to compute full or half ring polymer pot
       and forces.
    """

    def __init__(self):
        self.pot = None
        self.f = None
        pass

    def bind(self, dumop):
        self.dbeads = dumop.beads.copy()
        self.temp = dumop.temp
        if dumop.options["mode"] == 'rate':
            self.omega2 = (self.temp * (2 * self.dbeads.nbeads) * units.Constants.kb / units.Constants.hbar) ** 2
        elif dumop.options["mode"] == 'splitting':
            self.omega2 = (self.temp * (self.dbeads.nbeads) * units.Constants.kb / units.Constants.hbar) ** 2

        if dumop.options["opt"] == 'nichols' or dumop.options["opt"] == 'NR':
            self.h = self.spring_hessian(self.dbeads.natoms, self.dbeads.nbeads, self.dbeads.m3[0], self.omega2)

    def save(self, e, g):
        self.pot = e
        self.f = -g

    @staticmethod
    def spring_hessian(natoms, nbeads, m3, omega2, mode='half'):
        """Compute the 'spring hessian'

           OUT    h       = hessian with only the spring terms ('spring hessian')
            """

        info(" @spring_hessian", verbosity.high)
        ii = natoms * 3
        h = np.zeros([ii * nbeads, ii * nbeads])

        if nbeads == 1:
            return h

        # Diagonal
        h_sp = m3 * omega2
        diag1 = np.diag(h_sp)
        diag2 = np.diag(2.0 * h_sp)

        if mode == 'half':
            i = 0
            h[i * ii:(i + 1) * ii, i * ii:(i + 1) * ii] += diag1
            i = nbeads - 1
            h[i * ii:(i + 1) * ii, i * ii:(i + 1) * ii] += diag1
            for i in range(1, nbeads - 1):
                h[i * ii:(i + 1) * ii, i * ii:(i + 1) * ii] += diag2
        elif mode == 'splitting' or mode == 'full':
            for i in range(0, nbeads):
                h[i * ii:(i + 1) * ii, i * ii:(i + 1) * ii] += diag2
        else:
            raise ValueError("We can't compute the spring hessian.")

        # Non-Diagonal
        ndiag = np.diag(-h_sp)
        # Quasi-band
        for i in range(0, nbeads - 1):
            h[i * ii:(i + 1) * ii, (i + 1) * ii:(i + 2) * ii] += ndiag
            h[(i + 1) * ii:(i + 2) * ii, i * ii:(i + 1) * ii] += ndiag

        # Corner
        if mode == 'full':
            h[0:ii, (nbeads - 1) * ii:(nbeads) * ii] += ndiag
            h[(nbeads - 1) * ii:(nbeads) * ii, 0:ii] += ndiag

        return h

    def __call__(self, x, ret=True):
        """Computes spring energy and gradient for instanton optimization step"""

        if x.shape[0] == 1:  # only one bead
            self.dbeads.q = x.copy()
            e = 0.0
            g = np.zeros(x.shape[1])
            self.save(e, g)
        else:
            self.dbeads.q = x.copy()
            e = 0.00
            g = np.zeros(self.dbeads.q.shape, float)
            for i in range(self.dbeads.nbeads - 1):
                dq = self.dbeads.q[i + 1, :] - self.dbeads.q[i, :]
                e += self.omega2 * 0.5 * np.dot(self.dbeads.m3[0] * dq, dq)
            for i in range(0, self.dbeads.nbeads - 1):
                g[i, :] += self.dbeads.m3[i, :] * self.omega2 * (self.dbeads.q[i, :] - self.dbeads.q[i + 1, :])
            for i in range(1, self.dbeads.nbeads):
                g[i, :] += self.dbeads.m3[i, :] * self.omega2 * (self.dbeads.q[i, :] - self.dbeads.q[i - 1, :])

            self.save(e, g)

        if ret:
            return e, g


class FullMapper(object):
    """Creation of the multi-dimensional function to compute the physical and the spring forces.
    """

    def __init__(self, im, gm):
        self.im = im
        self.gm = gm

    def __call__(self, x):
        e1, g1 = self.im(x)
        e2, g2 = self.gm(x)
        e = e1 + e2
        g = np.add(g1, g2)
        return (e, g)


class DummyOptimizer(dobject):
    """ Dummy class for all optimization classes """

    def __init__(self):
        """Initialises object for GradientMapper (physical potential, forces and hessian)
        and SpringMapper ( spring potential,forces and hessian) """

        self.options = {}  #Optimization options
        self.optarrays  = {} #Optimization arrays

        self.gm = GradientMapper()
        self.im = SpringMapper()
        self.fm = FullMapper(self.im, self.gm)
        # self.lm           = LineMapper(self.fm)
        self.exit = False
        self.init = False


    def step(self, step=None):
        """Dummy simulation time step which does nothing."""
        pass

    def bind(self, geop):
        """
        Bind optimization options and call bind function of Mappers (get beads, cell,forces)
        check whether force size,  Hessian size from  match system size
        """

        self.beads = geop.beads
        self.cell = geop.cell
        self.forces = geop.forces
        self.fixcom = geop.fixcom
        self.fixatoms = geop.fixatoms

        # The resize action must be done before the bind

        if geop.optarrays["old_x"].size != self.beads.q.size:
            if geop.optarrays["old_x"].size == 0:
                geop.optarrays["old_x"] = np.zeros((self.beads.nbeads, 3 * self.beads.natoms), float)
            else:
                raise ValueError("Old positions size does not match system size")
        if geop.optarrays["old_u"].size != self.beads.nbeads:
            if geop.optarrays["old_u"].size == 0:
                geop.optarrays["old_u"] = np.zeros(self.beads.nbeads, float)
            else:
                raise ValueError("Old potential energy size does not match system size")
        if geop.optarrays["old_f"].size != self.beads.q.size:
            if geop.optarrays["old_f"].size == 0:
                geop.optarrays["old_f"] = np.zeros((self.beads.nbeads, 3 * self.beads.natoms), float)
            else:
                raise ValueError("Old forces size does not match system size")

        # Temperature
        self.temp = geop.ensemble.temp
        if geop.ensemble.temp == -1.0 or geop.ensemble.temp == 1.0:
            # This is due to a little inconsistency on the default value
            if self.beads.nbeads != 1:
                raise ValueError("Temperature must be specified for an Instanton calculation ")

        # Optimization mode
        self.options["mode"] = geop.options["mode"]

        # Generic optimization
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
        self.gm.bind(self)
        self.optarrays["energy_shift"] = geop.optarrays["energy_shift"]

    def exitstep(self, fx, fx0, x, exitt, step):
        """ Exits the simulation step. Computes time, checks for convergence. """
        self.qtime += time.time()

        tolerances = self.options["tolerances"]

        info(' @Exit step: Energy difference: %.1e, (condition: %.1e)' % (np.absolute((fx - fx0) / self.beads.natoms), tolerances["energy"]), verbosity.low)
        info(' @Exit step: Maximum force component: %.1e, (condition: %.1e)' % (np.amax(np.absolute(self.forces.f + self.im.f)), tolerances["force"]), verbosity.low)
        info(' @Exit step: Maximum component step component: %.1e, (condition: %.1e)' % (x, tolerances["position"]), verbosity.low)

        if (np.absolute((fx - fx0) / self.beads.natoms) <= tolerances["energy"]) \
                and ((np.amax(np.absolute(self.forces.f + self.im.f)) <= tolerances["force"]) or
                     (np.linalg.norm(self.forces.f.flatten() - self.optarrays["old_f"].flatten()) <= 1e-08)) \
                and (x <= tolerances["position"]):

            print_instanton_geo(self.options["prefix"] + '_FINAL', step, self.im.dbeads.nbeads, self.im.dbeads.natoms, self.im.dbeads.names,
                                self.im.dbeads.q, self.optarrays["old_u"], self.cell, self.optarrays["energy_shift"])

            if self.options["hessian_final"] != 'true':
                info("We are not going to compute the final hessian.", verbosity.low)
                info("Warning, The current hessian is not the real hessian is only an approximation .", verbosity.low)

            else:
                info("We are going to compute the final hessian", verbosity.low)
                get_hessian(self.optarrays["hessian"], self.gm, self.im.dbeads.q)
                print_instanton_hess(self.options["prefix"] + '_FINAL', step, self.optarrays["hessian"])

            exitt = True
            # If we just exit here, the last step (including the last hessian) will not be in the RESTART file

        return exitt


class HessianOptimizer(DummyOptimizer):
    """ Class that implements an optimization when a hessian is available."""

    def bind(self, geop):
        # call bind function from DummyOptimizer
        super(HessianOptimizer, self).bind(geop)

        self.options["hessian_update"] = geop.options["hessian_update"]
        self.options["hessian_asr"] = geop.options["hessian_asr"]
        self.options["hessian_init"] = geop.options["hessian_init"]

        self.im.bind(self)

        self.optarrays["initial_hessian"] = None

        if geop.optarrays["hessian"].size != (self.beads.natoms * 3 * self.beads.q.size):
            if geop.optarrays["hessian"].size == (self.beads.natoms * 3)**2:
                self.optarrays["initial_hessian"] = geop.optarrays["hessian"].copy()
                geop.optarrays["hessian"] = np.zeros((self.beads.natoms * 3, self.beads.q.size), float)

            elif geop.optarrays["hessian"].size == 0 and geop.options["hessian_init"] == 'true':
                info(" Initial hessian is not provided. We are going to compute it.", verbosity.low)
                geop.optarrays["hessian"] = np.zeros((self.beads.natoms * 3, self.beads.q.size))

                if ((self.beads.q - self.beads.q[0]) == 0).all() and self.beads.nbeads > 1:
                    raise ValueError("""We need a initial hessian in order to create our initial
                    instanton geometry. Please provide a (1-bead) hessian or an initial instanton geometry.""")

            else:
                raise ValueError(" 'Hessian_init' is false, an initial hessian (of the proper size) must be provided.")

        self.optarrays["hessian"] = geop.optarrays["hessian"]

    def initialize(self, step):

        if step == 0:

            info(" @GEOP: Initializing INSTANTON", verbosity.low)

            if self.beads.nbeads == 1:
                info(" @GEOP: Classical TS search", verbosity.low)
                if self.options["hessian_init"] == 'true':
                    get_hessian(self.optarrays["hessian"], self.gm, self.beads.q)
            else:
                # If the coordinates in all the imaginary time slices are the same
                if ((self.beads.q - self.beads.q[0]) == 0).all():
                    info(" @GEOP: We stretch the initial geometry with an 'amplitud' of %4.2f" 
                         % self.optarrays["delta"],verbosity.low)
                    
                    imvector = get_imvector(self.optarrays["initial_hessian"], self.beads.m3[0].flatten())
                    for i in range(self.beads.nbeads):
                        self.beads.q[i, :] += self.optarrays["delta"] * np.cos(
                            i * np.pi / float(self.beads.nbeads - 1)) * imvector[:]
                    if self.options["hessian_init"] != 'true':
                        info(
                            " @GEOP: Hessian_init isn't true but we have stretched the polymer so we are going to "
                            "compute the initial hessian anyway.",
                            verbosity.low)
                        self.options["hessian_init"] = 'true'
                else:
                    info(" @GEOP: Starting from the provided geometry in the extended phase space", verbosity.low)
                    if not (self.optarrays["initial_hessian"] is None):
                        raise ValueError(" You have to provided a hessian with size (3 x natoms)^2 but also geometry in"
                                         " the extended phase space (nbeads>1). Please check the inputs\n")

                if self.options["hessian_init"] == 'true':
                    info(" @GEOP: We are computing the initial hessian", verbosity.low)
                    get_hessian(self.optarrays["hessian"], self.gm, self.beads.q)

            # Update positions and forces
            self.optarrays["old_x"][:] = self.beads.q
            self.optarrays["old_u"][:] = self.forces.pots
            self.optarrays["old_f"][:] = self.forces.f

        if type(self.im.f) == type(None):
            self.im(self.beads.q, ret=False)  # Init instanton mapper

        if (self.optarrays["old_x"] == np.zeros((self.beads.nbeads, 3 * self.beads.natoms), float)).all():
            self.optarrays["old_x"][:] = self.beads.q

        self.init=True

    def step(self, step=None):
        """Dummy simulation time step which does nothing."""
        pass


class NicholsOptimizer(HessianOptimizer):
    """ Class that implements a nichols optimizations. It can find first order saddle points or minima"""

    def bind(self, geop):
        # call bind function from HessianOptimizer
        super(NicholsOptimizer, self).bind(geop)

    def step(self, step=None):
        """ Does one simulation time step."""

        if not self.init:
           self.initialize(step)

        self.qtime = -time.time()
        info("\n Instanton optimization STEP %d" % step, verbosity.low)

        if self.exit:
            softexit.trigger("Geometry optimization converged. Exiting simulation")

        # Construct hessian and get eigenvalues and eigenvector.
        # First construct complete hessian from reduced
        h0 = red2comp(self.optarrays["hessian"], self.beads.nbeads, self.beads.natoms)
        h1 = np.add(self.im.h, h0)  # add spring terms to the physical hessian
        d, w = clean_hessian(h1, self.im.dbeads.q, self.im.dbeads.natoms,
                             self.im.dbeads.nbeads, self.im.dbeads.m, self.im.dbeads.m3, self.options["hessian_asr"])

        # Find new movement direction
        if self.options["mode"] == 'rate':
            d_x = nichols(self.optarrays["old_f"], self.im.f, d, w, self.beads.m3, self.optarrays["big_step"])
        elif self.options["mode"] == 'splitting':
            d_x = nichols(self.optarrays["old_f"], self.im.f, d, w, self.beads.m3, self.optarrays["big_step"], mode=0)

        # Rescale step if necessary
        d_x_max = np.amax(np.absolute(d_x))
        info(" @Instanton: Current step norm = %g" % d_x_max, verbosity.medium)
        if np.amax(np.absolute(d_x)) > self.optarrays["big_step"]:
            info(" @Instanton: Attempted step norm = %g, scaled down to %g" 
                 % (d_x_max, self.optarrays["big_step"]), verbosity.low)
            
            d_x *= self.optarrays["big_step"] / np.amax(np.absolute(d_x_max))

        # Make movement and get new energy (u)  and forces(f) using mapper
        x = self.optarrays["old_x"] + d_x
        self.im(x, ret=False)  # Only to update the mapper
        u, g2 = self.gm(x)
        f = -g2

        # Update hessian
        if self.options["hessian_update"] == 'powell':
            d_g = np.subtract(self.optarrays["old_f"], f)
            i = self.im.dbeads.natoms * 3
            for j in range(self.im.dbeads.nbeads):
                aux = self.optarrays["hessian"][:, j * i:(j + 1) * i]
                dg = d_g[j, :]
                dx = d_x[j, :]
                Powell(dx, dg, aux)
        elif self.options["hessian_update"] == 'recompute':
            get_hessian(self.optarrays["hessian"], self.gm, x)

        # Update positions and forces
        self.beads.q = self.gm.dbeads.q
        self.forces.transfer_forces(self.gm.dforces)  # This forces the update of the forces

        # Print current instanton geometry and hessian
        if (self.options["save"] > 0 and np.mod(step, self.options["save"]) == 0) or self.exit:
            print_instanton_geo(self.options["prefix"], step, self.im.dbeads.nbeads, self.im.dbeads.natoms,
                                self.im.dbeads.names,self.im.dbeads.q, self.optarrays["old_u"], self.cell,
                                self.optarrays["energy_shift"])
            
            print_instanton_hess(self.options["prefix"], step, self.optarrays["hessian"])

        # Check Exit
        d_x_max = np.amax(np.absolute(np.subtract(self.beads.q, self.optarrays["old_x"])))
        self.exit = self.exitstep(self.forces.pot, self.optarrays["old_u"].sum(), d_x_max, self.exit, step)

        # Update positions and forces
        self.optarrays["old_x"][:] = self.beads.q
        self.optarrays["old_u"][:] = self.forces.pots
        self.optarrays["old_f"][:] = self.forces.f


class NROptimizer(HessianOptimizer):
    """ Class that implements a Newton-Raphson optimizations. It can find first order saddle points or minima"""

    def bind(self, geop):
        # call bind function from HessianOptimizer
        super(NROptimizer, self).bind(geop)

    def step(self, step=None):
        """ Does one simulation time step."""

        if not self.init:
            self.initialize(step)

        self.qtime = -time.time()
        info("\n Instanton optimization STEP %d" % step, verbosity.low)

        if self.exit:
            softexit.trigger("Geometry optimization converged. Exiting simulation")

        h_up_band = banded_hessian(self.optarrays["hessian"], self.im)  # create upper band matrix
        f = (self.optarrays["old_f"] + self.im.f).reshape(self.beads.natoms * 3 * self.nbeads, 1)

        d_x = invmul_banded(h_up_band, f)
        d_x.shape = self.beads.q.shape

        # Rescale step if necessary
        d_x_max = np.amax(np.absolute(d_x))
        info(" @Instanton: Current step norm = %g" % d_x_max, verbosity.medium)
        if np.amax(np.absolute(d_x)) > self.optarrays["big_step"]:
            info(" @Instanton: Attempted step norm = %g, scaled down to %g"
                 % (d_x_max, self.optarrays["big_step"]), verbosity.low)

            d_x *= self.optarrays["big_step"] / np.amax(np.absolute(d_x_max))

        # Make movement and get new energy (u)  and forces(f) using mapper
        x = self.optarrays["old_x"] + d_x
        self.im(x, ret=False)  # Only to update the mapper
        u, g2 = self.gm(x)
        f = -g2

        # Update hessian
        if self.options["hessian_update"] == 'powell':
            d_g = np.subtract(self.optarrays["old_f"], f)
            i = self.im.dbeads.natoms * 3
            for j in range(self.im.dbeads.nbeads):
                aux = self.optarrays["hessian"][:, j * i:(j + 1) * i]
                dg = d_g[j, :]
                dx = d_x[j, :]
                Powell(dx, dg, aux)
        elif self.options["hessian_update"] == 'recompute':
            get_hessian(self.optarrays["hessian"], self.gm, x)

        # Update positions and forces
        self.beads.q = self.gm.dbeads.q
        self.forces.transfer_forces(self.gm.dforces)  # This forces the update of the forces

        # Print current instanton geometry and hessian
        if (self.options["save"] > 0 and np.mod(step, self.options["save"]) == 0) or self.exit:
            print_instanton_geo(self.options["prefix"], step, self.im.dbeads.nbeads, self.im.dbeads.natoms,
                                self.im.dbeads.names, self.im.dbeads.q, self.optarrays["old_u"], self.cell,
                                self.optarrays["energy_shift"])

            print_instanton_hess(self.options["prefix"], step, self.optarrays["hessian"])

        # Check Exit
        d_x_max = np.amax(np.absolute(np.subtract(self.beads.q, self.optarrays["old_x"])))
        self.exit = self.exitstep(self.forces.pot, self.optarrays["old_u"].sum(), d_x_max, self.exit, step)

        # Update positions and forces
        self.optarrays["old_x"][:] = self.beads.q
        self.optarrays["old_u"][:] = self.forces.pots
        self.optarrays["old_f"][:] = self.forces.f

class LBFGSOptimizer(DummyOptimizer):

    def bind(self, geop):
        # call bind function from DummyOptimizer
        super(LBFGSOptimizer, self).bind(geop)

        if geop.optarrays["hessian"].size == (self.beads.natoms * 3) ** 2:
            self.optarrays["initial_hessian"] = geop.optarrays["hessian"].copy()
            geop.optarrays["hessian"] = np.zeros((self.beads.natoms * 3, self.beads.q.size))

        if geop.options["hessian_final"] == 'true':
            self.options["hessian_asr"] = geop.options["hessian_asr"]
            if geop.optarrays["hessian.size"] == 0:
                geop.optarrays["hessian"] = np.zeros((self.beads.natoms * 3, self.beads.q.size))
            self.optarrays["hessian"] = geop.optarrays["hessian"]

        self.im.bind(self)

        # Specific for LBFGS
        self.options["corrections"] = geop.options["corrections"]
        self.options["ls_options"] = geop.options["ls_options"]
        if geop.optarrays["qlist"].size != (self.options["corrections"] * self.beads.q.size):
            if geop.optarrays["qlist"].size == 0:
                geop.optarrays["qlist"] = np.zeros((self.options["corrections"], self.beads.q.size), float)
            else:
                raise ValueError("qlist size does not match system size")
        if geop.optarrays["glist"].size != (self.options["corrections"] * self.beads.q.size):
            if geop.optarrays["glist"].size == 0:
                geop.optarrays["glist"] = np.zeros((self.corrections, self.beads.q.size), float)
            else:
                raise ValueError("qlist size does not match system size")

        self.optarrays["qlist"] = geop.optarrays["qlist"]
        self.optarrays["glist"] = geop.optarrays["glist"]

        if geop.options["scale"] not in [0, 1, 2]:
            raise ValueError("Scale option is not valid")

        self.options["scale"] = geop.options["scale"]

        if geop.optarrays["d"].size != self.beads.q.size:
            if geop.optarrays["d"].size == 0:
                geop.optarrays["d"] = np.zeros((self.beads.nbeads, 3 * self.beads.natoms), float)
            else:
                raise ValueError("Initial direction size does not match system size")

        self.optarrays["d"] = geop.optarrays["d"]

    def initialize(self, step):

        if step == 0:
            info(" @GEOP: Initializing INSTANTON", verbosity.low)

            if self.beads.nbeads == 1:
                raise ValueError("We can not perform an splitting calculation with nbeads =1")
                # get_hessian(self.hessian, self.gm, self.beads.q)
            else:
                if ((self.beads.q - self.beads.q[0]) == 0).all():
                    # If the coordinates in all the imaginary time slices are the same
                    info(" @GEOP: We stretch the initial geometry with an 'amplitud' of %4.2f" % self.delta, verbosity.low)
                    imvector = get_imvector(self.initial_hessian, self.beads.m3[0].flatten())
                    for i in range(self.beads.nbeads):

                        self.beads.q[i, :] += self.delta * np.cos(i * np.pi / float(self.beads.nbeads - 1)) * imvector[:]
                else:
                    info(" @GEOP: Starting from the provided geometry in the extended phase space", verbosity.low)

            # Update positions and forces
            self.old_x[:] = self.beads.q
            self.old_u[:] = self.forces.pots
            self.old_f[:] = self.forces.f

        # This must be done after the stretching and before the self.d.
        if type(self.im.f) == type(None):
            self.im(self.beads.q, ret=False)  # Init instanton mapper

        # Specific for LBFGS
        if np.linalg.norm(self.d) == 0.0:
            f = self.forces.f + self.im.f  # ALBERTO1
            self.d += dstrip(f) / np.sqrt(np.dot(f.flatten(), f.flatten()))

        if (self.old_x == np.zeros((self.beads.nbeads, 3 * self.beads.natoms), float)).all():
            self.old_x[:] = self.beads.q

    def step(self, step=None):
        """ Does one simulation time step."""

        self.qtime = -time.time()
        info("\n Instanton optimization STEP %d" % step, verbosity.low)

        if not self.init:
            self.initialize(step)

        if self.exit:
            softexit.trigger("Geometry optimization converged. Exiting simulation")

        e, g = self.fm(self.beads.q)
        fdf0 = (e, g)

        # Do one step. Update hessian for the new position. Update the position and force inside the mapper.
        L_BFGS(self.old_x, self.d, self.fm, self.qlist, self.glist,
               fdf0, self.big_step, self.ls_options["tolerance"] * self.tolerances["energy"],
               self.ls_options["iter"], self.corrections, self.scale, step)

        # Update positions and forces
        self.beads.q = self.gm.dbeads.q
        self.forces.transfer_forces(self.gm.dforces)  # This forces the update of the forces

        # Print current instanton geometry
        if (self.options["save"] > 0 and np.mod(step, self.options["save"]) == 0) or self.exit:
            print_instanton_geo(self.options["prefix"], step, self.im.dbeads.nbeads, self.im.dbeads.natoms,
                                self.im.dbeads.names, self.im.dbeads.q, self.optarrays["old_u"], self.cell,
                                self.optarrays["energy_shift"])

        # Check exit
        d_x_max = np.amax(np.absolute(np.subtract(self.beads.q, self.optarrays["old_x"])))
        self.exit = self.exitstep(self.forces.pot, self.optarrays["old_u"].sum(), d_x_max, self.exit, step)

        # Update positions and forces
        self.old_x[:] = self.beads.q
        self.old_u[:] = self.forces.pots
        self.old_f[:] = self.forces.f

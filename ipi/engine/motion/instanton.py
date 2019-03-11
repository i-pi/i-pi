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
        alt_out: (Alternative outpu) Prints different formatting of outputs for geometry, hessian and bead potential energies.
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

    def __init__(self, fixcom=False, fixatoms=None,
                 mode='None',
                 tolerances={"energy": 1e-5, "force": 1e-4, "position": 1e-3},
                 biggest_step=0.3,
                 old_pos=np.zeros(0, float),
                 old_pot=np.zeros(0, float),
                 old_force=np.zeros(0, float),
                 opt='None',
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
                 hessian_final='False',
                 energy_shift=np.zeros(0, float)):
        """Initialises InstantonMotion.
        """

        super(InstantonMotion, self).__init__(fixcom=fixcom, fixatoms=fixatoms)

        # Optimization mode
        self.mode = mode

        # Generic optimization
        self.big_step = biggest_step
        self.tolerances = tolerances

        self.old_x = old_pos
        self.old_u = old_pot
        self.old_f = old_force

        # Generic instanton
        self.save = alt_out
        self.prefix = prefix
        self.delta = delta
        self.hessian_final = hessian_final
        self.energy_shift = energy_shift

        # We set the default optimization algorithm depending on the mode.
        if mode == 'rate':
            if opt == 'None':
                opt = 'nichols'
            if opt == 'lbfgs':
                raise ValueError("lbfgs option is not compatible with a rate calculation")
            self.opt = opt

        elif mode == 'splitting':
            if opt == 'None':
                opt = 'lbfgs'
            self.opt = opt

        if opt == 'nichols' or opt == 'NR':
            self.optimizer = HessianOptimizer()
            self.hessian_init = hessian_init
            self.hessian = hessian
            self.hessian_update = hessian_update
            self.hessian_asr = hessian_asr

        if self.opt == 'lbfgs':
            self.optimizer = LBFGSOptimizer()
            self.hessian = hessian  # Only for initial (to spread) or final
            self.hessian_asr = hessian_asr
            self.corrections = corrections_lbfgs
            self.scale = scale_lbfgs
            self.qlist = qlist_lbfgs
            self.glist = glist_lbfgs
            self.ls_options = ls_options
            self.d = old_direction

        # Do we put a warning to say that NR use Scipy? ALBERTO
        if self.opt == 'NR':
            info("Note that we need scipy to use NR. If storage and diagonalization of the full hessian is not a problem use nichols even though it may not be as efficient.", verbosity.low)

    def bind(self, ens, beads, nm, cell, bforce, prng, omaker):
        """Binds beads, cell, bforce and prng to InstantonMotion

            Args:
            beads: The beads object from whcih the bead positions are taken.
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
        if dumop.mode == 'rate':
            self.omega2 = (self.temp * (2 * self.dbeads.nbeads) * units.Constants.kb / units.Constants.hbar) ** 2
        elif dumop.mode == 'splitting':
            self.omega2 = (self.temp * (self.dbeads.nbeads) * units.Constants.kb / units.Constants.hbar) ** 2

        if dumop.opt == 'nichols' or dumop.opt == 'NR':
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

        self.gm = GradientMapper()
        self.im = SpringMapper()
        self.fm = FullMapper(self.im, self.gm)
        # self.lm           = LineMapper(self.fm)
        self.exit = False

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
        self.output_maker = geop.output_maker

        # The resize action must be done before the bind
        if geop.old_x.size != self.beads.q.size:
            if geop.old_x.size == 0:
                geop.old_x = np.zeros((self.beads.nbeads, 3 * self.beads.natoms), float)
            else:
                raise ValueError("Old positions size does not match system size")
        if geop.old_u.size != self.beads.nbeads:
            if geop.old_u.size == 0:
                geop.old_u = np.zeros(self.beads.nbeads, float)
            else:
                raise ValueError("Old potential energy size does not match system size")
        if geop.old_f.size != self.beads.q.size:
            if geop.old_f.size == 0:
                geop.old_f = np.zeros((self.beads.nbeads, 3 * self.beads.natoms), float)
            else:
                raise ValueError("Old forces size does not match system size")

        # Temperature
        self.temp = geop.ensemble.temp
        if geop.ensemble.temp == -1.0 or geop.ensemble.temp == 1.0:  # This is due to a little inconsistency on the default value
            if self.beads.nbeads != 1:
                raise ValueError("Temperature must be specified for an Instanton calculation ")

        # Optimization mode
        self.mode = geop.mode

        # Generic optimization
        self.tolerances = geop.tolerances
        self.big_step = geop.big_step
        self.old_x = geop.old_x
        self.old_u = geop.old_u
        self.old_f = geop.old_f
        self.opt = geop.opt  # optimization algorithm

        # Generic instanton
        self.save = geop.save
        self.prefix = geop.prefix
        self.delta = geop.delta
        self.hessian_final = geop.hessian_final
        self.gm.bind(self)
        self.energy_shift = geop.energy_shift

    def exitstep(self, fx, fx0, x, exitt, step):
        """ Exits the simulation step. Computes time, checks for convergence. """
        self.qtime += time.time()

        #f = open('STEP', 'a+')
        #print >>f, 'STEP %i' % step
        #print >>f, 'Energy difference: %.1e, (condition: %.1e)' % (np.absolute((fx - fx0) / self.beads.natoms), self.tolerances["energy"] )
        #print >>f, 'Maximum force component: %.1e, (condition: %.1e)' % (np.amax(np.absolute(self.forces.f+self.im.f)), self.tolerances["force"])
        #print >>f, 'Maximum component step component: %.1e, (condition: %.1e)' % (x, self.tolerances["position"])
        #print >>f, ' '
        # f.close()

        info(' @Exit step: Energy difference: %.1e, (condition: %.1e)' % (np.absolute((fx - fx0) / self.beads.natoms), self.tolerances["energy"]), verbosity.low)
        info(' @Exit step: Maximum force component: %.1e, (condition: %.1e)' % (np.amax(np.absolute(self.forces.f + self.im.f)), self.tolerances["force"]), verbosity.low)
        info(' @Exit step: Maximum component step component: %.1e, (condition: %.1e)' % (x, self.tolerances["position"]), verbosity.low)

        if (np.absolute((fx - fx0) / self.beads.natoms) <= self.tolerances["energy"]) \
                and ((np.amax(np.absolute(self.forces.f + self.im.f)) <= self.tolerances["force"]) or
                     (np.linalg.norm(self.forces.f.flatten() - self.old_f.flatten()) <= 1e-08)) \
                and (x <= self.tolerances["position"]):

            print_instanton_geo(self.prefix + '_FINAL', step, self.im.dbeads.nbeads, self.im.dbeads.natoms, self.im.dbeads.names,
                                self.im.dbeads.q, self.old_u, self.cell, self.energy_shift, self.output_maker)

            if self.hessian_final != 'true':
                info("We are not going to compute the final hessian.", verbosity.low)
                info("Warning, The current hessian is not the real hessian is only an approximation .", verbosity.low)

            else:
                info("We are going to compute the final hessian", verbosity.low)
                get_hessian(self.hessian, self.gm, self.im.dbeads.q)
                print_instanton_hess(self.prefix + '_FINAL', step, self.hessian, self.output_maker)

            exitt = True  # If we just exit here, the last step (including the last hessian) will not be in the RESTART file

        return exitt


class HessianOptimizer(DummyOptimizer):
    """ INSTANTON Rate calculation"""

    def bind(self, geop):
        # call bind function from DummyOptimizer
        super(HessianOptimizer, self).bind(geop)

        # Specific for RateOptimizer

        self.hessian_update = geop.hessian_update
        self.hessian_asr = geop.hessian_asr
        self.hessian_init = geop.hessian_init
#        self.output_maker = geop.output_maker

        self.im.bind(self)

        # Hessian
        self.initial_hessian = None

        if geop.hessian.size != (self.beads.natoms * 3 * self.beads.q.size):
            if geop.hessian.size == (self.beads.natoms * 3)**2:
                self.initial_hessian = geop.hessian.copy()
                geop.hessian = np.zeros((self.beads.natoms * 3, self.beads.q.size), float)
            elif geop.hessian.size == 0 and geop.hessian_init == 'true':
                info(" Initial hessian is not provided. We are going to compute it.", verbosity.low)
                geop.hessian = np.zeros((self.beads.natoms * 3, self.beads.q.size))
                if ((self.beads.q - self.beads.q[0]) == 0).all() and self.beads.nbeads > 1:
                    raise ValueError("""We need a initial hessian in order to create our initial
                    instanton geometry. Please provide a (1-bead) hessian or an initial instanton geometry.""")
            else:
                raise ValueError(" 'Hessian_init' is false, an initial hessian (of the proper size) must be provided.")

        self.hessian = geop.hessian

    def step(self, step=None):
        """ Does one simulation time step."""

        self.qtime = -time.time()
        info("\n Instanton optimization STEP %d" % step, verbosity.low)

        if step == 0:
            info(" @GEOP: Initializing INSTANTON", verbosity.low)

            if self.beads.nbeads == 1:
                info(" @GEOP: Classical TS search", verbosity.low)
                if self.hessian_init == 'true':
                    get_hessian(self.hessian, self.gm, self.beads.q)
            else:
                if ((self.beads.q - self.beads.q[0]) == 0).all():  # If the coordinates in all the imaginary time slices are the same
                    info(" @GEOP: We stretch the initial geometry with an 'amplitud' of %4.2f" % self.delta, verbosity.low)
                    imvector = get_imvector(self.initial_hessian, self.beads.m3[0].flatten())
                    for i in range(self.beads.nbeads):
                        self.beads.q[i, :] += self.delta * np.cos(i * np.pi / float(self.beads.nbeads - 1)) * imvector[:]
                    if self.hessian_init != 'true':
                        info(" @GEOP: Hessian_init isn't true but we have stretched the polymer so we are going to compute the initial hessian anyway.", verbosity.low)
                        self.hessian_init = 'true'
                else:
                    info(" @GEOP: Starting from the provided geometry in the extended phase space", verbosity.low)
                    if not (self.initial_hessian is None):
                        raise ValueError(" You have to provided a hessian with size (3xnatoms)^2 but also geometry in the extended phase space (nbeads>1). Please check the inputs\n")

                if self.hessian_init == 'true':
                    info(" @GEOP: We are computing the initial hessian", verbosity.low)
                    get_hessian(self.hessian, self.gm, self.beads.q)

            # Update positions and forces
            self.old_x[:] = self.beads.q
            self.old_u[:] = self.forces.pots
            self.old_f[:] = self.forces.f

        if type(self.im.f) == type(None):
            self.im(self.beads.q, ret=False)  # Init instanton mapper

        if (self.old_x == np.zeros((self.beads.nbeads, 3 * self.beads.natoms), float)).all():
            self.old_x[:] = self.beads.q
        if self.exit:
            softexit.trigger("Geometry optimization converged. Exiting simulation")

        if len(self.fixatoms) > 0:
            for dqb in self.old_f:
                dqb[self.fixatoms * 3] = 0.0
                dqb[self.fixatoms * 3 + 1] = 0.0
                dqb[self.fixatoms * 3 + 2] = 0.0

        # Do one step. Update hessian for the new position. Update the position and force inside the mapper.
        Instanton(self.old_x, self.old_f, self.im.f, self.hessian, self.hessian_update, self.hessian_asr, self.im, self.gm, self.big_step, self.opt, self.mode)

        # Update positions and forces
        self.beads.q = self.gm.dbeads.q
        self.forces.transfer_forces(self.gm.dforces)  # This forces the update of the forces

        # Print current instanton geometry and hessian
        if (self.save > 0 and np.mod(step, self.save) == 0) or self.exit:
            print_instanton_geo(self.prefix, step, self.im.dbeads.nbeads, self.im.dbeads.natoms, self.im.dbeads.names,
                                self.im.dbeads.q, self.old_u, self.cell, self.energy_shift, self.output_maker)
            print_instanton_hess(self.prefix, step, self.hessian, self.output_maker)

        # Exit simulation step
        d_x_max = np.amax(np.absolute(np.subtract(self.beads.q, self.old_x)))
        self.exit = self.exitstep(self.forces.pot, self.old_u.sum(), d_x_max, self.exit, step)

        # Update positions and forces
        self.old_x[:] = self.beads.q
        self.old_u[:] = self.forces.pots
        self.old_f[:] = self.forces.f


def Instanton(x0, f0, f1, h, update, asr, im, gm, big_step, opt, m):
    """Do one step. Update hessian for the new position. Update the position and force inside the mapper.

       Input:  x0 = last positions
               f0 = last physical forces
               f1 = last spring forces
                h = physical hessian
           update = how to update the hessian
              asr = how to clean the hessian
               im = spring mapper
               gm = gradient  mapper
         big_step = limit on step length
              opt = optimization algorithm to use
              m   = type of calculation: rate or splitting"""

    info(" @Instanton_step", verbosity.high)

    if opt == 'nichols':
        # Construct hessian and get eigenvalues and eigenvector
        h0 = red2comp(h, im.dbeads.nbeads, im.dbeads.natoms)    # construct complete hessian from reduced
        h1 = np.add(im.h, h0)                                 # add spring terms to the physical hessian
        d, w = clean_hessian(h1, im.dbeads.q, im.dbeads.natoms, im.dbeads.nbeads, im.dbeads.m, im.dbeads.m3, asr)

        # Find new movement direction
        if m == 'rate':
            d_x = nichols(f0, f1, d, w, im.dbeads.m3, big_step)
        elif m == 'splitting':
            d_x = nichols(f0, f1, d, w, im.dbeads.m3, big_step, mode=0)

    elif opt == 'NR':
        h_up_band = banded_hessian(h, im)  # create upper band matrix
        f = (f0 + f1).reshape(im.dbeads.natoms * 3 * im.dbeads.nbeads, 1)

        d_x = invmul_banded(h_up_band, f)
        d_x.shape = im.dbeads.q.shape

    # Rescale step
    d_x_max = np.amax(np.absolute(d_x))
    info(" @Instanton: Current step norm = %g" % d_x_max, verbosity.medium)
    if np.amax(np.absolute(d_x)) > big_step:
        info(" @Instanton: Attempted step norm = %g, scaled down to %g" % (d_x_max, big_step), verbosity.low)
        d_x *= big_step / np.amax(np.absolute(d_x_max))

    # Make movement and get new energy (u)  and forces(f) using mapper
    x = x0 + d_x
    im(x, ret=False)  # Only to update the mapper
    u, g2 = gm(x)
    f = -g2

    # Update hessian
    if update == 'powell':
        d_g = np.subtract(f0, f)
        i = im.dbeads.natoms * 3
        for j in range(im.dbeads.nbeads):
            aux = h[:, j * i:(j + 1) * i]
            dg = d_g[j, :]
            dx = d_x[j, :]
            Powell(dx, dg, aux)
    elif update == 'recompute':
        get_hessian(h, gm, x)


class LBFGSOptimizer(DummyOptimizer):

    def bind(self, geop):
        # call bind function from DummyOptimizer
        super(LBFGSOptimizer, self).bind(geop)

        if geop.hessian.size == (self.beads.natoms * 3) ** 2:
            self.initial_hessian = geop.hessian.copy()
            geop.hessian = np.zeros((self.beads.natoms * 3, self.beads.q.size))
        if geop.hessian_final == 'true':
            self.hessian_asr = geop.hessian_asr
            if geop.hessian.size == 0:
                geop.hessian = np.zeros((self.beads.natoms * 3, self.beads.q.size))
            self.hessian = geop.hessian

        self.im.bind(self)

        # Specific for LBFGS
        self.corrections = geop.corrections
        self.ls_options = geop.ls_options

        if geop.qlist.size != (self.corrections * self.beads.q.size):
            if geop.qlist.size == 0:
                geop.qlist = np.zeros((self.corrections, self.beads.q.size), float)
            else:
                raise ValueError("qlist size does not match system size")
        if geop.glist.size != (self.corrections * self.beads.q.size):
            if geop.glist.size == 0:
                geop.glist = np.zeros((self.corrections, self.beads.q.size), float)
            else:
                raise ValueError("qlist size does not match system size")

        self.qlist = geop.qlist
        self.glist = geop.glist

        if geop.scale not in [0, 1, 2]:
            raise ValueError("Scale option is not valid")

        self.scale = geop.scale

        if geop.d.size != self.beads.q.size:
            if geop.d.size == 0:
                geop.d = np.zeros((self.beads.nbeads, 3 * self.beads.natoms), float)
            else:
                raise ValueError("Initial direction size does not match system size")

        self.d = geop.d

    def step(self, step=None):
        """ Does one simulation time step."""

        self.qtime = -time.time()
        info("\n Instanton optimization STEP %d" % step, verbosity.low)

        if step == 0:
            info(" @GEOP: Initializing INSTANTON", verbosity.low)

            if self.beads.nbeads == 1:
                raise ValueError("We can not perform an splitting calculation with nbeads =1")
                # get_hessian(self.hessian, self.gm, self.beads.q)
            else:
                if ((self.beads.q - self.beads.q[0]) == 0).all():  # If the coordinates in all the imaginary time slices are the same
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

        if self.exit:
            softexit.trigger("Geometry optimization converged. Exiting simulation")

        if len(self.fixatoms) > 0:
            for dqb in self.old_f:
                dqb[self.fixatoms * 3] = 0.0
                dqb[self.fixatoms * 3 + 1] = 0.0
                dqb[self.fixatoms * 3 + 2] = 0.0

        e, g = self.fm(self.beads.q)
        fdf0 = (e, g)

        # Do one step. Update hessian for the new position. Update the position and force inside the mapper.
        L_BFGS(self.old_x, self.d, self.fm, self.qlist, self.glist,
               fdf0, self.big_step, self.ls_options["tolerance"] * self.tolerances["energy"],
               self.ls_options["iter"], self.corrections, self.scale, step)
        # ALBERTO2

        # Update positions and forces
        self.beads.q = self.gm.dbeads.q
        self.forces.transfer_forces(self.gm.dforces)  # This forces the update of the forces

        # Exit simulation step
        d_x_max = np.amax(np.absolute(np.subtract(self.beads.q, self.old_x)))
        self.exit = self.exitstep(self.forces.pot, self.old_u.sum(), d_x_max, self.exit, step)

        # Update positions and forces
        self.old_x[:] = self.beads.q
        self.old_u[:] = self.forces.pots
        self.old_f[:] = self.forces.f

        # Print current instanton geometry and hessian
        if (self.save > 0 and np.mod(step, self.save) == 0) or self.exit:
            print_instanton_geo(self.prefix, step, self.im.dbeads.nbeads, self.im.dbeads.natoms, self.im.dbeads.names,
                                self.im.dbeads.q, self.old_u, self.cell, self.energy_shift, self.output_maker)

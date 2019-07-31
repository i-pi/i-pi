"""Contains the classes that deal with the different dynamics required in
different types of ensembles.

Holds the algorithms required for normal mode propagators, and the objects to
do the constant temperature and pressure algorithms. Also calculates the
appropriate conserved energy quantity for the ensemble of choice.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import time

import numpy as np

from ipi.engine.motion import Motion, Dynamics
from ipi.engine.motion.dynamics import DummyIntegrator

from ipi.utils.depend import *
from ipi.engine.thermostats import Thermostat
from ipi.engine.barostats import Barostat

class ConstrainedDynamics(Dynamics):

    """self (path integral) molecular dynamics class.

    Gives the standard methods and attributes needed in all the
    dynamics classes.

    Attributes:
        beads: A beads object giving the atoms positions.
        cell: A cell object giving the system box.
        forces: A forces object giving the virial and the forces acting on
            each bead.
        prng: A random number generator object.
        nm: An object which does the normal modes transformation.

    Depend objects:
        econs: The conserved energy quantity appropriate to the given
            ensemble. Depends on the various energy terms which make it up,
            which are different depending on the ensemble.he
        temp: The system temperature.
        dt: The timestep for the algorithms.
        ntemp: The simulation temperature. Will be nbeads times higher than
            the system temperature as PIMD calculations are done at this
            effective classical temperature.
    """

    def __init__(self, timestep, mode="nve", splitting="obabo",
                thermostat=None, barostat=None, fixcom=False, fixatoms=None,
                nmts=None, nsteps_geo=1, nsteps_o = 1,
                constraint_list=[], csolver=None):

        """Initialises a "dynamics" motion object.

        Args:
            dt: The timestep of the simulation algorithms.
            fixcom: An optional boolean which decides whether the centre of mass
                motion will be constrained or not. Defaults to False.
        """

        super(Dynamics, self).__init__(fixcom=fixcom, fixatoms=fixatoms)
        dd(self).dt = depend_value(name='dt', value=timestep)

        if thermostat is None:
            self.thermostat = Thermostat()
        else:
            self.thermostat = thermostat

        if nmts is None or len(nmts) == 0:
            dd(self).nmts = depend_array(name="nmts", value=np.asarray([1], int))
        else:
            dd(self).nmts = depend_array(name="nmts", value=np.asarray(nmts, int))

        if barostat is None:
            self.barostat = Barostat()
        else:
            self.barostat = barostat
        self.enstype = mode


        if nmts is None or len(nmts) == 0:
            if self.enstype == "nve":
                self.integrator = NVEIntegrator_constraint()
            elif self.enstype == "nvt":
                self.integrator = NVTIntegrator_constraint()
            else:
                ValueError("No valid ensemble provided")
        else:
            if self.enstype == "nve":
                self.integrator = NVEIntegrator_constraintMTS()
            elif self.enstype == "nvt":
                self.integrator = NVTIntegrator_constraintMTS()
            else:
                ValueError("No valid ensemble provided")


        # splitting mode for the integrators
        dd(self).splitting = depend_value(name='splitting', value=splitting)

        # constraints
        self.fixcom = fixcom
        if fixatoms is None:
            self.fixatoms = np.zeros(0, int)
        else:
            self.fixatoms = fixatoms

        # The list of constraints coming from the input is an actual list of independent
        # constraints, and should be treated as such
        if constraint_list is not None:
            self.constraint_list = constraint_list


        if csolver is None:
            self.csolver = ConstraintSolver(self.constraint_list, tolerance=0.0001,maxit=10000,norm_order=2)
        else:
            if csolver["norm_order"] == -1:
                norm_order = float('inf')
            else:
                norm_order = csolver["norm_order"]
            if csolver["mode"] == "full":
                self.csolver = ConstraintSolver(self.constraint_list, tolerance=csolver["tolerance"],
                                                maxit=csolver["maxit"],
                                                norm_order=norm_order)
            else:
                self.csolver = SparseConstraintSolver(self.constraint_list, tolerance=csolver["tolerance"],
                                                      maxit=csolver["maxit"],
                                                      norm_order=norm_order)



        # MS: should all these be initialized as dependent arrays/values.?
        self.nsteps_geo = nsteps_geo
        self.nsteps_o = nsteps_o


    def bind(self, ens, beads, nm, cell, bforce, prng, omaker):
        """Binds ensemble beads, cell, bforce, and prng to the dynamics.

        This takes a beads object, a cell object, a forcefield object and a
        random number generator object and makes them members of the ensemble.
        It also then creates the objects that will hold the data needed in the
        ensemble algorithms and the dependency network. Note that the conserved
        quantity is defined in the init, but as each ensemble has a different
        conserved quantity the dependencies are defined in bind.

        Args:
            beads: The beads object from whcih the bead positions are taken.
            nm: A normal modes object used to do the normal modes transformation.
            cell: The cell object from which the system box is taken.
            bforce: The forcefield object from which the force and virial are
                taken.
            prng: The random number generator object which controls random number
                generation.
        """
        super(ConstrainedDynamics, self).bind(ens, beads, nm, cell, bforce, prng, omaker)

        # now binds the constraints
        for c in self.constraint_list:
            c.bind(beads)
        self.csolver.bind(beads)

        """
        if len(self.nmts) > 1 or (len(self.nmts) == 1 and self.nmts[0] != 1):
            raise ValueError("MTS with constrains has not been implemented.")
        """

class ConstraintBase(dobject):
    """ Constraint class for MD. Base class."""

    def __init__(self, constrained_indices, constraint_values, ncons=0):
        self.constrained_indices = constrained_indices

        dd(self).constraint_values = depend_array(name="constraint_values",
                value=np.asarray(constraint_values).copy())
        if ncons == 0:
            self.ncons = len(constraint_values)

        # determines the list of unique indices of atoms that are affected by this constraint
        self.i_unique = np.unique(self.constrained_indices.flatten())
        self.mk_idmaps()

    def mk_idmaps(self):

        # makes lookup dictionary and lists to quickly access the portions of arrays that are affected by this constraint
        self.i_reverse = { value : i for i,value in enumerate(self.i_unique)}
        self.n_unique = len(self.i_unique)
        self.i3_unique = np.zeros(self.n_unique*3, int)
        for i in range(self.n_unique):
            self.i3_unique[3*i:3*(i+1)] = [3*self.i_unique[i], 3*self.i_unique[i]+1, 3*self.i_unique[i]+2]
        # this can be used to access the position array based on the constrained_indices list
        self.i3_indirect = np.zeros((len(self.constrained_indices.flatten()),3), int)
        for ri, i in enumerate(self.constrained_indices.flatten()):
            rri = self.i_reverse[i]
            self.i3_indirect[ri] = [3*rri, 3*rri+1, 3*rri+2]

    def bind(self, beads):

        self.beads = beads
        dself = dd(self)
        dself.m3 = depend_array(name="g", value=np.zeros(self.n_unique*3),
                    func=(lambda: self.beads.m3[0,self.i3_unique]),
                    dependencies=[dd(self.beads).m3])
        # The constraint function is computed at iteratively updated 
        # coordinates during geodesic integration; this array holds such
        # updates
        dself.q = depend_array(name="q", value=np.zeros(self.n_unique*3))
        dself.g = depend_array(name="g", value=np.zeros(self.ncons),
                               func=self.gfunc, 
                               dependencies=[dself.q, dself.constraint_values])
        # The gradient of the constraint function is computed at the configuration
        # obtained from the previous converged SHAKE step; this array holds
        # a local copy of that configuration
        dself.qprev = depend_array(name="qprev", 
                                   value=np.zeros(self.n_unique*3))
        dself.Dg = depend_array(name="Dg", 
                                value=np.zeros((self.ncons, self.n_unique*3)),
                                func=self.Dgfunc, 
                                dependencies=[dself.qprev,
                                              dself.constraint_values])
        dself.GramChol = depend_array(name="GramChol",
                                      value=np.zeros((self.ncons,self.ncons)),
                                      func=self.GCfunc, 
                                      dependencies=[dself.Dg] )

    def gfunc(self):
        raise NotImplementedError()

    def Dgfunc(self):
        """
        Calculates the Jacobian of the constraint.
        """
        raise NotImplementedError()

    def GCfunc(self):
        dg = dstrip(self.Dg).copy()
        dgm = dg/self.m3
        return np.linalg.cholesky(np.dot(dg, dgm.T))


class RigidBondConstraint(ConstraintBase):
    """ Constraint class for MD.
        Specialized for rigid bonds. This can actually hold a *list*
        of rigid bonds, i.e. there will be a list of pairs of atoms and
        a list of bond lengths. """

    def __init__(self,constrained_indices,constraint_values):

        super(RigidBondConstraint,self).__init__(constrained_indices, constraint_values)
        self.constrained_indices.shape = (self.ncons, 2)
        self.i3_indirect.shape = (self.ncons, 2, 3)

    def gfunc(self):
        """
        Calculates the constraint.
        """

        q = dstrip(self.q)
        r = np.zeros(self.ncons)
        constraint_distances = dstrip(self.constraint_values)
        for i in range(self.ncons):
            c_atoms = self.i3_indirect[i]
            c_dist = constraint_distances[i]
            #print q[c_atoms[0]], q[c_atoms[1]], c_dist
            r[i] = np.sum((q[c_atoms[0]] - q[c_atoms[1]])**2) - c_dist**2
        if q[0] == float('inf'):
            ValueError("fgfgf")
            print("autsch")
            exit()
        #print("gfunc", r)
        return r

    def Dgfunc(self, reduced=False):
        """
        Calculates the Jacobian of the constraint.
        """

        q = dstrip(self.qprev)
        #constrained_indices = self.constrained_indices
        r = np.zeros((self.ncons, self.n_unique*3))
        for i in range(self.ncons):
            c_atoms = self.i3_indirect[i]
            inst_position_vector = q[c_atoms[0]] - q[c_atoms[1]]
            r[i][c_atoms[0]] =   2.0 * inst_position_vector
            r[i][c_atoms[1]] = - 2.0 * inst_position_vector
        return r

    """
    def bind(self, integrator):
        '''
        Reference all the variables for simpler access
        '''
        super(RigidBondConstraint, self).bind(integrator)

        dself = dd(self)


        dself.constrained_indices = depend_value(name="constrained_indices", func=lambda: integrator.constrained_indices.reshape((len(integrator.constrained_indices) / 2, 2)))
        dself.constrained_distances = depend_value(name="constrained_distances", func=lambda: integrator.constrained_distances)
        dself.G = depend_array(name="G", func=self.gfunc, value=np.zeros(len(self.constrained_distances)), dependencies=[dself.beads.q, dself.constrained_indices, dself.constrained_distances])
        dself.DG = depend_array(name="DG", func=self.Dgfunc, value=np.zeros((len(self.constrained_distances), 3 * self.beads.natoms)), dependencies=[dself.G, dself.beads.q, dself.constrained_indices, dself.constrained_distances])
        #MS: why not dself.G
    """


class ConstraintList(ConstraintBase):
    """ Constraint class for MD"""

    def __init__(self, constraint_list):
        self.constraint_list = constraint_list
        self.ncons = sum([constr.ncons for constr in constraint_list])

        # determines the list of unique indices of atoms that are affected 
        # by this list of constraint
        self.i_unique = np.zeros(0)
        for c in self.constraint_list:
            self.i_unique = np.union1d(
                    c.constrained_indices.flatten(),self.i_unique
                    )
        self.constrained_indices = self.i_unique
        self.mk_idmaps()

        # must now find the mapping from the unique indices in each constraint 
        # ic_map[i] gives the position where the atoms involved 
        # in the i-th constraint are stored in the compact list
        self.ic_map = []
        self.ic3_map = []
        for c in self.constraint_list:
            i_map = np.array([self.i_reverse[i] for i in c.i_unique])
            self.ic_map.append(i_map)
            self.ic3_map.append(np.array([[i*3, i*3+1, i*3+2] for i in i_map]).flatten())

    def bind(self, beads):

        # this is special because it doesn't hold constraint_values so we have to really specialize
        self.beads = beads
        dself = dd(self)
        # this holds all of the atoms in this list of constraints
        dself.q = depend_array(name="q", value=np.zeros(self.n_unique*3))
        # this holds the configurations of the listed atom obtained
        # at the end of the previous step
        dself.qprev = depend_array(name="qprev", value=np.zeros(self.n_unique*3))
        dself.g = depend_array(name="g", value=np.zeros(self.ncons),
                               func=self.gfunc)
        dself.Dg = depend_array(name="Dg", 
                                value=np.zeros((self.ncons, self.n_unique*3)), 
                                func=self.Dgfunc)
        dself.GramChol = depend_array(name="GramChol", 
                                      value=np.zeros((self.ncons,self.ncons)),
                                      func=self.GCfunc, dependencies=[dself.Dg])
        # we link all of the sub-constraints to the lists of unique q and qprev,
        # so that each constraint gets automatically updated
        def make_qgetter(k):
            return lambda: self.q[self.ic3_map[k]]
        def make_qprevgetter(k):
            return lambda: self.qprev[self.ic3_map[k]]
        for ic, c in enumerate(self.constraint_list):
            c.bind(beads)
            # deal with constraint values
            dq = dd(c).q
            dq.add_dependency(dself.q)
            dq._func = make_qgetter(ic)
            dself.g.add_dependency(dd(c).g)
            # ...and gradients
            dqprev = dd(c).qprev
            dqprev.add_dependency(dself.qprev)
            dqprev._func = make_qprevgetter(ic)
            dself.Dg.add_dependency(dd(c).Dg)

    def gfunc(self):
        """
        Compute the constraint function.
        """
        r = np.zeros(self.ncons)
        si = 0
        for constr in self.constraint_list:
            r[si:si+constr.ncons] = constr.g
            si += constr.ncons
        return r

    def Dgfunc(self):
        """
        Compute the Jacobian of the constraint function.
        """
        r = np.zeros((self.ncons, np.size(q)))
        si = 0
        for ic, constr in enumerate(self.constraint_list):
            r[si:si+constr.ncons, self.ic_map[ic]] = constr.Dg
            si += constr.ncons
        return r

    def get_iai(self):
        iai = []
        for constr in self.constraint_list:
            iai += list(constr.get_iai())
        return np.unique(iai)

class ConstraintSolverBase(dobject):

    def __init__(self, constraint_list):
        self.constraint_list = constraint_list

    def update_constraints(self, beads):
        raise NotImplementedError()

    def proj_cotangent(self, beads):
        raise NotImplementedError()

    def proj_manifold(self, beads, stepsize=None, proj_p=True):
        raise NotImplementedError()


class ConstraintSolver(ConstraintSolverBase):

    def __init__(self, constraints, tolerance=0.0001,maxit=1000,norm_order=2):
        super(ConstraintSolver,self).__init__(constraints)

        self.tolerance = tolerance
        self.maxit = maxit
        self.norm_order = norm_order

    def update_constraints(self, beads):

        self.Dg = self.constraints.Dgfunc(dstrip(beads.q[0]))
        self.Gram = np.dot(self.Dg, np.dot(np.diagflat(1.0/beads.m3[0]), np.transpose(self.Dg)))
        self.GramChol = np.linalg.cholesky(self.Gram)
        self.ciu = True

    def proj_cotangent(self, beads):


        m3 = dstrip(beads.m3[0])
        p = dstrip(beads.p[0]).copy()



        b = np.dot(self.Dg, p/m3[0])
        x = np.linalg.solve(np.transpose(self.GramChol),np.linalg.solve(self.GramChol, b))
        p += - np.dot(np.transpose(self.Dg),x)

        beads.p[0] = p

    def proj_manifold(self, beads, stepsize=None, proj_p=True):
        '''
        projects onto Manifold using the Gram matrix defined by self.Dg and self.Gram
        '''
        m3 = dstrip(beads.m3[0])
        p = dstrip(beads.p[0]).copy()
        q = dstrip(beads.q[0]).copy()

        if proj_p and stepsize is None:
            stepsize = self.dt

        i = 0
        if self.constraints.ncons > 0:
            self.g = self.constraints.gfunc(q)
            while (i < self.maxit and self.tolerance <= np.linalg.norm(self.g, ord=self.norm_order)):
                dlambda = np.linalg.solve(np.transpose(self.GramChol),np.linalg.solve(self.GramChol, self.g))
                delta = np.dot(np.transpose(self.Dg),dlambda)
                q += - delta /m3
                self.g = self.constraints.gfunc(q)
                if proj_p:
                    update_diff = - delta / stepsize
                    p += update_diff
                i += 1

                if (i == self.maxit):
                    print('No convergence in Newton iteration for positional component');

        beads.p[0] = p
        beads.q[0] = q


class SparseConstraintSolver(ConstraintSolverBase):

    def __init__(self, constraint_list, tolerance=0.001,maxit=1000,norm_order=2):
        super(SparseConstraintSolver,self).__init__(constraint_list)

        self.tolerance = tolerance
        self.maxit = maxit
        self.norm_order = norm_order
        #self.ic_list = [constr.get_ic() for constr in self.constraint_list]

    def bind(self, beads):

        self.beads = beads

        # sets the initial value of the constraint positions
        q = dstrip(beads.q[0])
        for constr in self.constraint_list:
            constr.qprev = q[constr.i3_unique.flatten()]

    def update_constraints(self, beads):

        #m3 = dstrip(beads.m3[0])
        #self.Dg_list = [constr.Dg[:,constr] for constr in self.constraint_list ]
        #self.Gram_list = [np.dot(dg,(dg/m3[ic]).T) for dg, ic in zip(self.Dg_list,self.ic_list)]
        #self.GramChol_list = [ np.linalg.cholesky(G) for G in self.Gram_list]
        self.ciu = True


    def proj_cotangent(self, beads):

        m3 = dstrip(beads.m3[0])
        p = dstrip(beads.p[0]).copy()

        if len(self.constraint_list) > 0:
            #print("number of constraints: ", len(self.constraint_list))
            for constr in self.constraint_list : #zip(self.Dg_list, self.ic_list, self.GramChol_list):
                dg = dstrip(constr.Dg)
                ic = constr.i3_unique
                gramchol = dstrip(constr.GramChol)

                b = np.dot(dg, p[ic]/constr.m3)
                x = np.linalg.solve(np.transpose(gramchol),np.linalg.solve(gramchol, b))
                p[ic] += - np.dot(np.transpose(dg),x)
        beads.p[0] = p

    def proj_manifold(self, beads, stepsize=None, proj_p=True):
        '''
        projects onto Manifold using the Gram matrix defined by self.Dg and self.Gram
        '''
        m3 = dstrip(beads.m3[0])
        p = dstrip(beads.p[0]).copy()
        q = dstrip(beads.q[0]).copy()

        #for constr in self.constraint_list:
                #print "before", constr.qprev
                #print "vs", q[constr.i3_unique.flatten()]

        i = 0
        if len(self.constraint_list) > 0:
            for constr in self.constraint_list: # zip(self.Dg_list, self.ic_list, self.GramChol_list, self.constraint_list):
                # these must only be computed on the manifold so we store them and don't update them
                dg = dstrip(constr.Dg)#.copy()
                gramchol = dstrip(constr.GramChol)#.copy()
                ic = constr.i3_unique
                constr.q = q[ic]
                g = dstrip(constr.g)
                #print "g vector", g
                while (i < self.maxit and self.tolerance <= np.linalg.norm(g, ord=self.norm_order)):
                    dlambda = np.linalg.solve(np.transpose(gramchol),np.linalg.solve(gramchol, g))
                    delta = np.dot(np.transpose(dg),dlambda)
                    q[ic] += - delta / m3[ic]
                    constr.q = q[ic] # updates the constraint to recompute g
                    g = dstrip(constr.g)
                    if proj_p:
                        update_diff = - delta / stepsize
                        p[ic] += update_diff
                        i += 1
                    if (i == self.maxit):
                        print('No convergence in Newton iteration for positional component');
            # in theory each constraint update should be independent but we cannot be sure...
            # after all constraints have been applied, q is on the manifold and we can update the constraint
            # positions
            for constr in self.constraint_list:
                constr.qprev = q[constr.i3_unique.flatten()]
                #print "after", constr.qprev

        beads.p[0] = p
        beads.q[0] = q


class ConstrainedIntegrator(DummyIntegrator):
    """ No-op integrator for (PI)MD """

    def __init__(self):
        # Fix attribute input
        #print "**************** ConstrainedIntegrator init **************"
        super(ConstrainedIntegrator,self).__init__()


    def get_qdt(self):
        return self.dt * 0.5

    def get_pdt(self):
        return np.array([self.dt * 0.5])

    def get_tdt(self):
        tdt = super(ConstrainedIntegrator,self).get_tdt()/self.nsteps_o
        return super(ConstrainedIntegrator,self).get_tdt()/self.nsteps_o

    def bind(self, motion):
        """ Reference all the variables for simpler access."""

        super(ConstrainedIntegrator,self).bind(motion)

        self.constraint_list = motion.constraint_list
        self.ciu = False

        if motion.nsteps_geo is None or len(motion.nsteps_geo) == 0:
            dd(self).nsteps_geo = depend_array(name="nsteps_geo", value=np.asarray([1], int))
        else:
            dd(self).nsteps_geo = depend_array(name="nsteps_geo", value=np.asarray(motion.nsteps_geo, int))

        if motion.nsteps_o is None or len(motion.nsteps_o) == 0:
            dd(self).nsteps_o = depend_array(name="nsteps_o", value=np.asarray([1], int))
        else:
            dd(self).nsteps_o = depend_array(name="nsteps_o", value=np.asarray(motion.nsteps_o, int))

        dd(self).csolver = motion.csolver


        #tolerance = depend_value(name="tolerance", func=lambda: motion.tolerance)
        #dd(self).maxit = depend_value(name="maxit", func=lambda: motion.maxit)
        #dd(self).norm_order = depend_value(name="norm_order", func=lambda: motion.norm_order)

    #print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Number of constraints: ", self.constraints.ncons)


    def step_A(self, stepsize=None):
        if stepsize is None:
            stepsize = self.qdt
        """Unconstrained A-step"""
        self.beads.q[0] += self.beads.p[0] / dstrip(self.beads.m3)[0]  * stepsize
        self.ciu = False

    def step_B(self, stepsize=None):
        if stepsize is None:
            stepsize = self.pdt[0]
        """Unconstrained B-step"""
        self.beads.p[0] += self.forces.forces_mts(0)[0] * stepsize

    def step_Bc(self, stepsize=None):
        """Unconstrained B-step followed by a projection into the cotangent space"""
        self.step_B(stepsize)
        self.proj_cotangent(self.beads)

    def step_BAc(self, stepsize=None):
        """
        half unconstrained B step and full A step followed by a projection onto the manifold
        """
        if stepsize is None:
            stepsize = self.dt

        if not self.ciu:
            self.update_constraints(self.beads)

        self.step_B(.5 * stepsize)
        self.step_A(stepsize)
        self.proj_manifold(self.beads, stepsize)

    def update_constraints(self, beads):
        pass

    def proj_cotangent(self, beads):
        self.csolver.proj_cotangent(beads)


    def proj_manifold(self, beads, stepsize=None, proj_p=True):
        if proj_p and stepsize is None:
            stepsize = self.dt
        self.csolver.proj_manifold(beads, stepsize, proj_p)


    def step_Ag(self, stepsize=None, Nsteps=None):
        '''
        Geodesic flow

        When called:
        -self.d_params is  assumed to satisfy contraint
        '''
        if stepsize is None:
            stepsize = self.qdt
        if Nsteps is None:
            Nsteps = self.n_geoflow

        substepsize = (1.0 * stepsize)/Nsteps

        '''
        Resolve momentum constraint and update Gram matrix if neccesary
        '''


        self.proj_cotangent(self.beads)

        for i in range(Nsteps):

            self.step_A(substepsize)
            self.proj_manifold(self.beads, substepsize)
            self.update_constraints(self.beads)
            self.proj_cotangent(self.beads)


    def step(self, step=None):
        """Dummy simulation time step which raises an error if called."""
        raise NotImplementedError()

    def pconstraints(self):
        """Dummy centroid momentum step which does nothing."""
        pass



class NVEIntegrator_constraint(ConstrainedIntegrator):



    """ Integrator object for constant energy simulations.

    Has the relevant conserved quantity and normal mode propagator for the
    constant energy ensemble. Note that a temperature of some kind must be
    defined so that the spring potential can be calculated.

    Attributes:
        ptime: The time taken in updating the velocities.
        qtime: The time taken in updating the positions.
        ttime: The time taken in applying the thermostat steps.

    Depend objects:
        econs: Conserved energy quantity. Depends on the bead kinetic and
            potential energy, and the spring potential energy.
    """

    def __init__(self):
        super(NVEIntegrator_constraint,self).__init__()


    def bind(self, motion):
        """ Reference all the variables for simpler access."""

        if len(motion.nmts) > 1 or motion.nmts[0] != 1 :
            raise ValueError("NVEIntegrator_constraint does not support multiple time stepping. Use NVEIntegrator_constraintMTS instead")

        ConstrainedIntegrator.bind(self, motion)


    def step(self, step=None):
        """Does one simulation time step."""

        #print("NVE step called")
        if self.splitting == "rattle":
             #self.constraint_step()

            self.step_BAc(self.dt)
            self.step_Bc(.5 * self.dt)
            #print("rattle")

        elif self.splitting == "geodesic":
            self.step_Bc(.5 * self.dt)
            self.step_Ag(stepsize=self.dt, Nsteps=self.nsteps_geo)
            self.step_Bc(.5 * self.dt)
            #print("geodesic")

class NVTIntegrator_constraint(NVEIntegrator_constraint):

    """Integrator object for constant temperature simulations.

    Has the relevant conserved quantity and normal mode propagator for the
    constant temperature ensemble. Contains a thermostat object containing the
    algorithms to keep the temperature constant.

    Attributes:
        thermostat: A thermostat object to keep the temperature constant.
    """
    def __init__(self):
        #print "**************** NVEIntegrator_constraint init **************"
        super(NVTIntegrator_constraint,self).__init__()
    #print("~~~~~~~~~~~~~~~~~~~~~~ tau = ", self.thermostat.tau)

    def bind(self, motion):
        """ Reference all the variables for simpler access."""

        if len(motion.nmts) > 1 or motion.nmts[0] != 1 :
            raise ValueError("NVTIntegrator_constraint does not support multiple time stepping. Use NVTIntegrator_constraintMTS instead")

        ConstrainedIntegrator.bind(self, motion)


    def step_Oc(self):

        for i in xrange(self.nsteps_o):
            self.thermostat.step()

            p = dstrip(self.beads.p).copy()
            sm = dstrip(self.thermostat.sm)
            p /= sm
            self.ensemble.eens += np.dot(p.flatten(), p.flatten()) * 0.5

            self.proj_cotangent(self.beads)

            p = dstrip(self.beads.p).copy()
            p /= sm
            self.ensemble.eens -= np.dot(p.flatten(), p.flatten()) * 0.5


    def step(self, step=None):
        """Does one simulation time step."""
        #print("~~~~~~~~~~~~~~~~~~~~~~ tau = ", self.thermostat.tau)
        #print("~~~~~~~~~~~~~~~~~~~~~~ dt = ", self.thermostat.dt)
        if self.splitting == "obabo":
            self.step_Oc()
            self.step_Bc(.5 * self.dt)
            self.step_Ag(self.dt,Nsteps=self.nsteps_geo)
            self.step_Bc(.5 * self.dt)
            self.step_Oc()


        elif self.splitting == "baoab":
            self.step_Bc(.5 * self.dt)
            self.step_Ag(.5 * self.dt,Nsteps=self.nsteps_geo)
            self.step_Oc()
            self.step_Ag(.5 * self.dt,Nsteps=self.nsteps_geo)
            self.step_Bc(.5 * self.dt)


class ConstrainedIntegratorMTS(NVTIntegrator_constraint):

    def bind(self, motion):
        """ Reference all the variables for simpler access."""
        ConstrainedIntegrator.bind(self, motion)


    def step_A(self, stepsize=None, level=0):
        if stepsize is None:
            stepsize = self.qdt/np.prod(self.nmts[:(level+1)])

        super(ConstrainedIntegratorMTS,self).step_A(stepsize)

    def step_B(self, stepsize=None, level=0):
        if stepsize is None:
            stepsize = self.pdt[0]/np.prod(self.nmts[:(level+1)])
        """Unconstrained B-step"""
        self.beads.p[0] += self.forces.forces_mts(level)[0] * stepsize

    def step_Bc(self, stepsize=None, level=0):
        """Unconstrained B-step followed by a projection into the cotangent space"""
        self.step_B(stepsize=stepsize, level=level)
        self.proj_cotangent(self.beads)

    def step_BAc(self, stepsize=None, level=0):
        """
            half unconstrained B step and full A step followed by a projection onto the manifold
        """
        if stepsize is None:
            stepsize = self.dt/np.prod(self.nmts[:(level+1)])

        if not self.ciu:
            self.update_constraints(self.beads)

        self.step_B(.5 * stepsize,level=level)
        self.step_A(stepsize)
        self.proj_manifold(self.beads, stepsize)


    def step_Ag(self, stepsize=None, Nsteps=None, level=0):
        '''
        Geodesic flow

        When called: -self.d_params is  assumed to satisfy contraint
        '''
        if stepsize is None:
            stepsize = self.qdt/np.prod(self.nmts[:(level+1)])
        if Nsteps is None:
            Nsteps = self.n_geoflow

        super(ConstrainedIntegratorMTS,self).step_Ag(stepsize=stepsize, Nsteps=Nsteps)

    def step_rattle(self, stepsize=None, level=0):
        if stepsize is None:
            stepsize = self.dt/np.prod(self.nmts[:(level+1)])

        self.step_BAc(stepsize=stepsize, level=level)
        self.step_Bc(stepsize = .5 * stepsize, level=level)

    def step_respa(self, level=0):

        #print("level: ", level)
        stepsize = self.dt/np.prod(self.nmts[:(level+1)])

        #print("stepsize: ", stepsize)
        self.step_Bc(stepsize = .5 * stepsize, level = level)
        #self.step_B(stepsize = .5 * stepsize, level = level)
        if level+1 < np.size(self.nmts):
            for i in range(self.nmts[level+1]):
                self.step_respa(level=level+1)
        else:
            self.step_A(stepsize)
            self.proj_manifold(self.beads, stepsize=stepsize, proj_p=True)
        #self.step_B(stepsize = .5 * stepsize, level = level)
        self.step_Bc(stepsize = .5 * stepsize, level = level)

    def step_grespa(self, level=0):

        stepsize = self.dt/np.prod(self.nmts[:(level+1)])

        self.step_Bc(stepsize = .5 * stepsize, level = level)
        if level+1 < np.size(self.nmts):
            for i in range(self.nmts[level+1]):
                self.step_grespa(level=level+1)
        else:
            self.step_Ag(stepsize, Nsteps=self.nsteps_geo[level])
        self.step_Bc(stepsize = .5 * stepsize, level = level)

    def step(self, step=None):
        raise NotImplementedError()

class NVEIntegrator_constraintMTS(ConstrainedIntegratorMTS):

    def step(self, step=None):
        """Does one simulation time step."""
        #print("~~~~~~~~~~~~~~~~~~~~~~ tau = ", self.thermostat.tau)
        #print("~~~~~~~~~~~~~~~~~~~~~~ dt = ", self.thermostat.dt)
        if self.splitting == "respa":
            self.step_respa()


        elif self.splitting == "grespa":
            self.step_grespa()

        else:
            raise ValueError("No valid splitting method spefified for NVE integration. Choose \"respa\" or \"grespa\".")


class NVTIntegrator_constraintMTS(ConstrainedIntegratorMTS):


    def get_tdt(self):
        if self.splitting in ["o-respa-o", "o-grespa-o"]:
            return self.dt * 0.5
        else:
            raise ValueError("Invalid splitting requested. Only \"o-respa-o\" and \"o-grespa-o\" are supported.")

    def step(self, step=None):
        """Does one simulation time step."""
        #print("~~~~~~~~~~~~~~~~~~~~~~ tau = ", self.thermostat.tau)
        #print("~~~~~~~~~~~~~~~~~~~~~~ dt = ", self.thermostat.dt)
        if self.splitting == "o-respa-o":
            self.step_Oc()
            self.step_respa()
            self.step_Oc()

        elif self.splitting == "o-grespa-o":
            self.step_grespa()

        else:
            raise ValueError("No valid splitting method spefified for NVE integration. Choose \"o-respa-o\" or \"o-grespa-o\".")

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

        # Create constraint object
        if constraint_list is not None:
            self.constraints = ConstraintList(constraint_list)

            
        if csolver is None:
            self.csolver = ConstraintSolver(self.constraints, tolerance=0.0001,maxit=10000,norm_order=2)
        else:
            if csolver["norm_order"] == -1:
                norm_order = float('inf')
            else:
                norm_order = csolver["norm_order"]
            if csolver["mode"] == "full":
                self.csolver = ConstraintSolver(self.constraints, tolerance=csolver["tolerance"],
                                                maxit=csolver["maxit"],
                                                norm_order=norm_order)
            else:
                self.csolver = SparseConstraintSolver(self.constraints, tolerance=csolver["tolerance"],
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

        """
        if len(self.nmts) > 1 or (len(self.nmts) == 1 and self.nmts[0] != 1):
            raise ValueError("MTS with constrains has not been implemented.")
        """

class ConstraintBase(dobject):
    """ Constraint class for MD. Base class."""

    def __init__(self, ncons):
        self.ncons = ncons

    def gfunc(self,q):
        raise NotImplementedError()

    def Dgfunc(self,q):
        """
        Calculates the Jacobian of the constraint.
        """
        raise NotImplementedError()
    
    def get_ic(self):
        """
        returns the indices constrained coordinate indices
        """
        return [ i for j in self.get_iai() for i in range(3*j,3*j+3) ]

    def get_iai(self):
        """
        get indices of atoms to which the constraint applies. Implementation of this function is required if constraints are handled by a subclass of SparseConstraintSolver
            
        """
        raise NotImplementedError()


class RigidBondConstraint(ConstraintBase):
    """ Constraint class for MD.
        Specialized for rigid bonds. This can actually hold a *list*
        of rigid bonds, i.e. there will be a list of pairs of atoms and
        a list of bond lengths. """

    def __init__(self,constrained_indices,constrained_distances):

        super(RigidBondConstraint,self).__init__(ncons=len(constrained_distances))
        self.constrained_indices = constrained_indices.reshape([self.ncons, 2])
        self.constrained_distances = constrained_distances
        
        self.iai = self.get_iai()
        self.iai_inv = { value : i for i,value in enumerate(self.iai)}
    
    def get_iai(self):
        return np.unique(self.constrained_indices.flatten())
    

    def gfunc(self, q):
        """
        Calculates the constraint.
        """
        r = np.zeros(self.ncons)
        constrained_indices = dstrip(self.constrained_indices)
        constrained_distances = dstrip(self.constrained_distances)
        for i in range(self.ncons):
            c_atoms = constrained_indices[i]
            c_dist = constrained_distances[i]
            #print c_atoms, self.constrained_indices
            r[i] = np.sum((q[c_atoms[0] * 3:c_atoms[0] * 3 + 3] - q[c_atoms[1] * 3: c_atoms[1] * 3 + 3])**2) - c_dist**2
                #print("q", beads.q[0])
        if q[0] == float('inf'):
            ValueError("fgfgf")
            print("autsch")
            exit()
        #print("gfunc", r)
        return r

    def Dgfunc(self, q, reduced=False):
        """
        Calculates the Jacobian of the constraint.
        """
        constrained_indices = dstrip(self.constrained_indices)
        r = np.zeros([self.ncons, np.size(q)])
        for i in range(self.ncons):
            c_atoms = constrained_indices[i]
            inst_position_vector = q[c_atoms[0] * 3:c_atoms[0] * 3 + 3] - q[c_atoms[1] * 3 : c_atoms[1] * 3 + 3]
            r[i,c_atoms[0] * 3 : c_atoms[0] * 3 + 3] =   2.0 * inst_position_vector
            r[i,c_atoms[1] * 3 : c_atoms[1] * 3 + 3] = - 2.0 * inst_position_vector
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

    def gfunc(self, q):
        """
        Compute the constraint function.
        """
        r = np.zeros(self.ncons)
        si = 0
        for constr in self.constraint_list:
            r[si:si+constr.ncons] = constr.gfunc(q)
            si += constr.ncons
        return r

    def Dgfunc(self, q):
        """
        Compute the Jacobian of the constraint function.
        """
        r = np.zeros((self.ncons, np.size(q)))
        si = 0
        for constr in self.constraint_list:
            r[si:si+constr.ncons,:] = constr.Dgfunc(q)
            si += constr.ncons
        return r

    def get_iai(self):
        iai = []
        for constr in self.constraint_list:
            iai += list(constr.get_iai())
        return np.unique(iai)

class ConstraintSolverBase(dobject):
    
    def __init__(self, constraints):
        self.constraints = constraints

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
    
    def __init__(self, constraints, tolerance=0.001,maxit=1000,norm_order=2):
        super(SparseConstraintSolver,self).__init__(constraints)
        
        self.tolerance = tolerance
        self.maxit = maxit
        self.norm_order = norm_order
        
        self.constraint_list = constraints.constraint_list
        self.ic_list = [constr.get_ic() for constr in self.constraint_list]

    
    def update_constraints(self, beads):
        
        m3 = dstrip(beads.m3[0])
        self.Dg_list = [constr.Dgfunc(dstrip(beads.q[0]))[:,ic] for constr, ic in zip(self.constraint_list,self.ic_list) ]
        self.Gram_list = [np.dot(dg,(dg/m3[ic]).T) for dg, ic in zip(self.Dg_list,self.ic_list)]
        self.GramChol_list = [ np.linalg.cholesky(G) for G in self.Gram_list]
        self.ciu = True
                             
                                                  
    def proj_cotangent(self, beads):
        
        m3 = dstrip(beads.m3[0])
        p = dstrip(beads.p[0]).copy()
        
        if len(self.constraint_list) > 0:
            #print("number of constraints: ", len(self.constraint_list))
            for dg, ic, gramchol in zip(self.Dg_list, self.ic_list, self.GramChol_list):
                b = np.dot(dg, p[ic]/m3[ic])
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
        
        i = 0
        if len(self.constraint_list) > 0:
            for dg, ic, gramchol, constr in zip(self.Dg_list, self.ic_list, self.GramChol_list, self.constraint_list):
                g = constr.gfunc(q)
                while (i < self.maxit and self.tolerance <= np.linalg.norm(g, ord=self.norm_order)):
                    dlambda = np.linalg.solve(np.transpose(gramchol),np.linalg.solve(gramchol, g))
                    delta = np.dot(np.transpose(dg),dlambda)
                    q[ic] += - delta /m3[ic]
                    g = constr.gfunc(q)
                    if proj_p:
                        update_diff = - delta / stepsize
                        p[ic] += update_diff
                        i += 1
                    if (i == self.maxit):
                        print('No convergence in Newton iteration for positional component');
    
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

        self.constraints = motion.constraints
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
        if self.constraints.ncons > 0:
            self.csolver.update_constraints(beads)
            self.ciu = True
    
    def proj_cotangent(self, beads):
        if self.constraints.ncons > 0:
            if not self.ciu:
                self.update_constraints(beads)
            self.csolver.proj_cotangent(beads)


    def proj_manifold(self, beads, stepsize=None, proj_p=True):
        if self.constraints.ncons > 0:
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
            self.proj_cotangent(self.beads)

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

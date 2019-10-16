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
from ipi.utils.constrtools import *
from ipi.engine.thermostats import Thermostat
from ipi.engine.barostats import Barostat

class ConstrainedDynamics(Dynamics):
    """Constrained molecular dynamics class.

    Provides the basic infrastructure to use constraints in molecular dynamics simulations.

    Attributes (on top of those inherited):
        nsteps_geo: number of geodesic integrator iterations
        nsteps_o: number of steps for the stochastic thermostat integrator 
                (needed because O does not commute with the constraint)
        constraint_list: list of constraints that must be applied
        csolver: solver to be used for the constraining

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
                nmts=None, nsteps_geo=1, nsteps_o=1,
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

        if self.enstype == "nve":
            self.integrator = NVEConstrainedIntegrator()
        elif self.enstype == "nvt":
            self.integrator = NVTConstrainedIntegrator()
        else:
            ValueError("{:s} is not a valid ensemble for constrained dynamics".format(self.enstype))

        # splitting mode for the integrators
        dd(self).splitting = depend_value(name='splitting', value=splitting)

        # GT: these are currently not used in constrained propagation
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
            self.csolver = ConstraintSolver(
                    self.constraint_list, 
                    tolerance=csolver["tolerance"],
                    maxit=csolver["maxit"],
                    norm_order=norm_order)

        # parameters of the geodesic integrator. will probably never need    
        # to be changed dynamically, so for the moment we don't consider them depend objects
        self.nsteps_geo = nsteps_geo
        self.nsteps_o = nsteps_o


    def bind(self, ens, beads, nm, cell, bforce, prng, omaker):
        """Binds the constrained dynamics object.
        Just passes control to the parent class, and then binds the constraints
        and the solver. 
        """
        super(ConstrainedDynamics, self).bind(ens, beads, nm, cell, bforce, prng, omaker)

        # now binds the constraints
        for c in self.constraint_list:
            c.bind(beads)
        self.csolver.bind(beads)


class ConstraintSolverBase(dobject):
    """ Empty base class for the constraint solver. Provides the interface
    that must be used to offer constraint functionalities to an integrator.    
    """

    def __init__(self, constraint_list, dt=1.0):
        self.constraint_list = constraint_list
        
        # time step - will have to be linked to the dynamics time step
        dd(self).dt = depend_value(name="dt", value=dt)

    def bind(self, beads):

        self.beads = beads
        # Sets the initial value of the constraint positions
        q = dstrip(beads.q[0])
        for constr in self.constraint_list:
            constr.qprev = q[constr.i3_unique.flatten()]
            
    def proj_cotangent(self, beads):
        """Set the momenta conjugate to the constrained degrees of freedom
        to zero by projecting onto the cotangent space of the constraint
        manifold.
        """
        
        raise NotImplementedError()

    def proj_manifold(self, beads, stepsize=None, proj_p=True):
        """Enforce the constraints on the positions (project onto the 
        constraint manifold using the Gramian matrix).
        """
                    
        raise NotImplementedError()


class ConstraintSolver(ConstraintSolverBase):
    """ An implementation of a constraint solver that uses M-RATTLE to
    impose constraints onto the momenta and M-SHAKE to impose constraints
    onto the positions.
    """

    def __init__(self, constraint_list, dt=1.0, 
                 tolerance=0.001, maxit=1000, norm_order=2):
        super(ConstraintSolver,self).__init__(constraint_list, dt)
        self.tolerance = tolerance
        self.maxit = maxit
        self.norm_order = norm_order

    def proj_cotangent(self):
        
        m3 = dstrip(self.beads.m3[0])
        p = dstrip(self.beads.p[0]).copy()
        self.beads.p.hold()
        
        for constr in self.constraint_list:
            dg = dstrip(constr.Dg)
            ic = constr.i3_unique
            gramchol = dstrip(constr.GramChol)
            b = np.dot(dg, p[ic]/constr.m3)
            x = np.linalg.solve( 
                    np.transpose(gramchol),np.linalg.solve(gramchol, b))
            p[ic] -= np.dot(np.transpose(dg),x)
        self.beads.p[0] = p
        self.beads.p.resume()

    def proj_manifold(self):
        
        m3 = dstrip(self.beads.m3[0])
        p = dstrip(self.beads.p[0]).copy()
        q = dstrip(self.beads.q[0]).copy()
        for constr in self.constraint_list:
            dg = dstrip(constr.Dg)
            gramchol = dstrip(constr.GramChol)
            ic = constr.i3_unique
            constr.q = q[ic]
            g = dstrip(constr.g)
            i = 0
            while (i < self.maxit and self.tolerance <= np.linalg.norm(g, ord=self.norm_order)):
                dlambda = np.linalg.solve( 
                        np.transpose(gramchol), 
                        np.linalg.solve(gramchol, g))
                delta = np.dot(np.transpose(dg), dlambda)
                q[ic] -= delta / m3[ic]
                constr.q = q[ic] # updates the constraint to recompute g
                g = dstrip(constr.g)
                p[ic] -= delta / self.dt
                i += 1
                if (i == self.maxit):
                    print('No convergence in Newton iteration for positional component');
        # after all constraints have been applied, q is on the manifold and we can update the constraint positions
        for constr in self.constraint_list:
            constr.qprev = q[constr.i3_unique.flatten()]
        self.beads.p[0] = p
        self.beads.q[0] = q


class ConstrainedIntegrator(DummyIntegrator):
    """ No-op integrator for classical constrained propagation """

    def __init__(self):
        super(ConstrainedIntegrator,self).__init__()
        
    def pconstraints(self):
        """Dummy centroid momentum step which does nothing."""
        
        if len(self.fixatoms) > 0:
            raise ValueError("Cannot fix atoms together with constrained MD")
        super(ConstrainedIntegrator,self).pconstraints()

    def get_qdt(self):
        # get the base dt for doing 
        return super(ConstrainedIntegrator,self).get_qdt()/self.nsteps_geo

    def get_tdt(self):
        return super(ConstrainedIntegrator,self).get_tdt()/self.nsteps_o

    def bind(self, motion):
        """ Reference all the variables for simpler access."""

        super(ConstrainedIntegrator,self).bind(motion)

        self.constraint_list = motion.constraint_list
        dself = dd(self)
        if motion.nsteps_geo is None:
            dself.nsteps_geo = depend_value(name="nsteps_geo", value=1)
        else:
            dself.nsteps_geo = depend_value(name="nsteps_geo", value=motion.nsteps_geo)
        print self.nsteps_geo
        dself.qdt.add_dependency(dself.nsteps_geo)
        if motion.nsteps_o is None:
            dself.nsteps_o = depend_value(name="nsteps_o", value=1)
        else:
            dself.nsteps_o = depend_value(name="nsteps_o", value=motion.nsteps_o)
        print self.nsteps_o
        dself.tdt.add_dependency(dself.nsteps_o)
        dd(self).csolver = motion.csolver
        dpipe(dself.qdt, dd(self.csolver).dt)
    
    def proj_cotangent(self):
        self.csolver.proj_cotangent()

    def proj_manifold(self):
        self.csolver.proj_manifold()
        

class NVEConstrainedIntegrator(ConstrainedIntegrator):
    """An propagator for constant-energy integration under a set of
    constraints.
    """
    
    def step_A(self):
        """Unconstrained A-step"""
        self.beads.q[0] += self.beads.p[0] / dstrip(self.beads.m3)[0] * self.qdt

    def step_B(self, level=0):
        """Unconstrained B-step"""
        self.beads.p[0] += self.forces.forces_mts(level)[0] * self.pdt[level]

    def step_Bc(self, level=0):
        """Unconstrained B-step followed by a projection into the cotangent space"""
        self.step_B(level)
        self.proj_cotangent()

    def step_Ag(self):
        """
        Geodesic flow 
        """
        # Resolve momentum constraint and update Gram matrix if neccesary
        # GT: is the following line necessary?
        self.proj_cotangent()

        for i in range(self.nsteps_geo):
            self.step_A()
            self.proj_manifold()
            self.proj_cotangent()
            
    def mtsprop_ba(self, index):
        """ Recursive MTS step -- this is adapted directly from the 
        NVEIntegrator class"""

        mk = int(self.nmts[index] / 2)

        for i in range(mk):  # do nmts/2 full sub-steps
            self.step_Bc(index)
            self.pconstraints() # currently does nothing
            if index == self.nmtslevels - 1:
                # call Q propagation for dt/alpha at the inner step
                self.step_Ag()
                self.step_Ag()
            else:
                self.mtsprop(index + 1)

            self.step_Bc(index)
            self.pconstraints()

        if self.nmts[index] % 2 == 1:
            # propagate p for dt/2alpha with force at level index
            self.step_Bc(index)
            self.pconstraints()
            if index == self.nmtslevels - 1:
                # call Q propagation for dt/alpha at the inner step
                self.step_Ag()
            else:
                self.mtsprop_ba(index + 1)

    def mtsprop_ab(self, index):
        """ Recursive MTS step -- this is adapted directly from the 
        NVEIntegrator class"""

     
        if self.nmts[index] % 2 == 1:
            if index == self.nmtslevels - 1:
                # call Q propagation for dt/alpha at the inner step
                self.step_Ag()
            else:
                self.mtsprop_ab(index + 1)
                
            # propagate p for dt/2alpha with force at level index
            self.step_Bc(index)
            self.pconstraints()
            
        mk = int(self.nmts[index] / 2)
        for i in range(mk):  # do nmts/2 full sub-steps
            self.step_Bc(index)
            self.pconstraints()
            if index == self.nmtslevels - 1:
                # call Q propagation for dt/alpha at the inner step
                self.step_Ag()
                self.step_Ag()
            else:
                self.mtsprop(index + 1)

            self.step_Bc(index)
            self.pconstraints()
            
    def mtsprop(self, index):
        # just calls the two pieces together
        self.mtsprop_ba(index)
        self.mtsprop_ab(index)

    def step(self, step=None):
        """Does one simulation time step."""

        self.mtsprop(0)
    

class NVTConstrainedIntegrator(NVEConstrainedIntegrator):

    """Constrained integrator object for constant temperature simulations.

    Has the relevant conserved quantity for the constant temperature 
    ensemble. Contains a thermostat object that keeps the temperature
    constant.

    Attributes:
        thermostat: A thermostat object to keep the temperature constant.
    """

    def step_Oc(self):
        for i in xrange(self.nsteps_o):
            self.thermostat.step()
            p = dstrip(self.beads.p).copy()
            sm = dstrip(self.thermostat.sm)
            p /= sm
            self.ensemble.eens += np.dot(p.flatten(), p.flatten()) * 0.5
            self.proj_cotangent()
            p = dstrip(self.beads.p).copy()
            p /= sm
            self.ensemble.eens -= np.dot(p.flatten(), p.flatten()) * 0.5

    def step(self, step=None):
        """Does one simulation time step."""
        
        if self.splitting == "obabo":
            self.step_Oc()
            self.pconstraints()
            self.mtsprop(0)
            self.step_Oc()
            self.pconstraints()

        elif self.splitting == "baoab":
            self.mtsprop_ba(0)
            self.step_Oc()
            self.pconstraints()
            self.mtsprop_ab(0)

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

    def __init__(self, timestep, mode="nve", splitting="obabo", thermostat=None, barostat=None, fixcom=False, fixatoms=None, nmts=None, nsteps_geo=1, nsteps_o = 1, constrained_indices=None, constrained_distances=None):

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
            self.integrator = NVEIntegrator_constraint()
        elif self.enstype == "nvt":
            self.integrator = NVTIntegrator_constraint()
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
        if constrained_indices is not None:
            self.constraint = RigidBondConstraint(constrained_indices,constrained_distances)

        self.nsteps_geo = nsteps_geo
        self.nsteps_o = nsteps_o
        #print "constrained_dynamics constructor"

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

        if len(self.nmts) > 1 or (len(self.nmts) == 1 and self.nmts[0] != 1):
            raise ValueError("MTS with constrains has not been implemented.")

class Constraint(dobject):
    """ Constraint class for MD"""

    def __init__(self, ncons):
        self.ncons = ncons

    def gfunc(self):
        raise NotImplementedError()

    def Dgfunc(self):
        """
        Calculates the Jacobian of the constraint.
        """
        raise NotImplementedError()

class ConstraintList(Constraint):
    """ Constraint class for MD"""
    
    def __init__(self, constraint_list):
        self.constraint_list = constraint_list
        self.ncons = sum([constr.ncons for constr in constraint_list])
    
    
    
    def gfunc(self, beads)
        """
        Compute the constraint function.
        """
        r = np.zeros(self.ncons)
        si = 0
        for constr in constraint_list:
            r[si:si+constr.ncons] = constr.gfunc(beads)
            si += constr.ncons
        return r
    
    def Dgfunc(self, beads):
        """
        Compute the Jacobian of the constraint function.
        """
        r = np.zeros((self.ncons, 3 * beads.natoms))
        for constr in constraint_list:
            r[si:si+constr.ncons,:] = constr.Gfunc(beads)
            si += constr.ncons
        return r



class RigidBondConstraint(Constraint):
    """ Constraint class for MD"""

    def __init__(self,constrained_indices,constrained_distances):

        super(RigidBondConstraint,self).__init__(ncons=len(constrained_distances))
        self.constrained_indices = constrained_indices.reshape([self.ncons, 2])
        self.constrained_distances = constrained_distances

        #print "RigidBondConstraint bla bla"

    def gfunc(self, beads):
        """
        Calculates the constraint.
        """
        r = np.zeros(self.ncons)
        for i in xrange(self.ncons):
            c_atoms = self.constrained_indices[i]
            c_dist = self.constrained_distances[i]
            #print c_atoms, self.constrained_indices
            r[i] = np.sum((beads.q[0][c_atoms[0] * 3:c_atoms[0] * 3 + 3] - beads.q[0][c_atoms[1] * 3: c_atoms[1] * 3 + 3])**2) - c_dist**2
                #print("q", beads.q[0])
        if dd(beads).q[0][0] == float('inf'):
            ValueError("fgfgf")
            print("autsch")
            exit()
        #print("gfunc", r)
        return r

    def Dgfunc(self, beads):
        """
        Calculates the Jacobian of the constraint.
        """

        r = np.zeros((self.ncons, 3 * beads.natoms))
        for i in xrange(self.ncons):
            c_atoms = self.constrained_indices[i]
            #print c_atoms, self.constrained_indices
            #print c_atoms[0]
            #print c_atoms[1]
            #print beads.q[0]
            inst_position_vector = beads.q[0][c_atoms[0] * 3:c_atoms[0] * 3 + 3] - beads.q[0][c_atoms[1] * 3 : c_atoms[1] * 3 + 3]
            r[i][c_atoms[0] * 3 : c_atoms[0] * 3 + 3] =   2.0 * inst_position_vector
            r[i][c_atoms[1] * 3 : c_atoms[1] * 3 + 3] = - 2.0 * inst_position_vector
        return r
        #MS: why c_atoms[0] * 3, c_atoms[0] * 3 + 3 ? (the commma...)

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

    def step(self,beads):
        """
        Dummy step.
        """

        for i in xrange(len(self.constrained_distances)):
            c_atoms = self.constrained_indices[i]
            c_dist = self.constrained_distances[i]
            #print "Constrained atoms : ", c_atoms
            #print "Current distance :",  np.linalg.norm(beads.q[0][c_atoms[0] * 3:c_atoms[0] * 3 + 3] - beads.q[0][c_atoms[1] * 3: c_atoms[1] * 3 + 3])
            #print "Target distance :", c_dist


class ConstrainedIntegrator(DummyIntegrator):
    """ No-op integrator for (PI)MD """

    def __init__(self, tol=0.00001, maxit=float('inf')):
        # Fix attribute input
        #print "**************** ConstrainedIntegrator init **************"
        super(ConstrainedIntegrator,self).__init__()
        self.tol = tol
        self.maxit = maxit


    def get_qdt(self):
        return self.dt * 0.5

    def get_pdt(self):
        return np.array([self.dt * 0.5])

    def get_tdt(self):
        tdt = super(ConstrainedIntegrator,self).get_tdt()/self.nsteps_o
        print("nsteps_o", self.nsteps_o)
        print("tdt", tdt)
        return super(ConstrainedIntegrator,self).get_tdt()/self.nsteps_o

    def bind(self, motion):
        """ Reference all the variables for simpler access."""
        #print "************************** ConstrainedIntegrator bind ****************************"
        if len(motion.nmts) > 1 or motion.nmts[0] != 1 :
             raise ValueError("Constrained integrator does not support multiple time stepping")

        super(ConstrainedIntegrator,self).bind(motion)

        self.constraint = motion.constraint # constraint is not a dependent obejct, thus not updated during simulation
        self.ciu = False

        dd(self).nsteps_geo = depend_value(name="nsteps_geo", func=lambda: motion.nsteps_geo)
        dd(self).nsteps_o = depend_value(name="nsteps_o", func=lambda: motion.nsteps_o)
    #print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Number of constraints: ", self.constraint.ncons)

    def update_constraints(self, beads):
        self.Dg = self.constraint.Dgfunc(beads)
        self.Gram = np.dot(self.Dg, np.dot(np.diagflat(1.0/beads.m3[0]), np.transpose(self.Dg)))
        self.GramChol = np.linalg.cholesky(self.Gram)
        self.ciu = True

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
    """
    def proj_cotangent(self, beads):
        pass

    def proj_manifold(self, beads, stepsize=None):
        pass
    """
    def proj_cotangent(self, beads):
        if not self.ciu:
            self.update_constraints(beads)

        if self.constraint.ncons > 0:
            b = np.dot(self.Dg, self.beads.p[0]/beads.m3[0])
            x = np.linalg.solve(np.transpose(self.GramChol),np.linalg.solve(self.GramChol, b))
            beads.p[0] += - np.dot(np.transpose(self.Dg),x)


    def proj_manifold(self, beads, stepsize=None):
        '''
        projects onto Manifold using the Gram matrix defined by self.Dg and self.Gram
        '''
        if stepsize is None:
            stepsize = self.dt

        self.g = self.constraint.gfunc(beads)
        #print("beads.q[0] before",  beads.q[0])
        #print("g :", self.g)
        #print("Dg :",self.Dg)
        #print("Gram :",self.Gram)
        #print("GramChol :",self.GramChol)

        i = 0
        if self.constraint.ncons > 0:
            while (i < self.maxit and self.tol <= np.linalg.norm(self.g)):
                #print("proj_manifold index", i)
                #print("proj_manifold error", np.linalg.norm(self.g))
                #print(self.GramChol)
                #print("g, i :", self.g, i)

                self.constraint.step(beads)
                self.g = self.constraint.gfunc(beads)
                dlambda = np.linalg.solve(np.transpose(self.GramChol),np.linalg.solve(self.GramChol, self.g))
                #print("dlambda", dlambda)
                delta = np.dot(np.transpose(self.Dg),dlambda)
                beads.q[0] += - delta /beads.m3[0]
                #self.d_momentum += - delta/self.p_stepsize
                self.g = self.constraint.gfunc(beads)
                #print("g, i after ud:", self.g, i)
                #print("beads.q[0] after",  beads.q[0])
                #print("delta",  delta)
                #print("stepsize", stepsize)
                #print("beads.m3[0]", beads.m3[0])
                update_diff = - delta / stepsize
                #print("update_diff", update_diff)
                beads.p[0] += update_diff
                i += 1
        if (i == self.maxit):
            print('No convergence in Newton iteration for positional component');

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


        if not self.ciu:
            self.update_constraints(self.beads)
        self.proj_cotangent(self.beads)

        for i in range(Nsteps):

            self.step_A(substepsize)
            self.update_constraints(self.beads)
            self.proj_manifold(self.beads, substepsize)
            self.proj_cotangent(self.beads)


    def step(self, step=None):
        """Dummy simulation time step which raises an error if called."""
        raise NotImplementedError()

    def pconstraints(self):
        """Dummy centroid momentum step which does nothing."""
        pass

    def constraint_step(self):
        """
            Dummy constraint step.
            """

        if len(self.constrained_distances) > 0:
            self.constraint.step()


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
        #print "**************** NVEIntegrator_constraint init **************"
        super(NVEIntegrator_constraint,self).__init__()


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



    def pconstraints(self):
        """This removes the centre of mass contribution to the kinetic energy.

        Calculates the centre of mass momenta, then removes the mass weighted
        contribution from each atom. If the ensemble defines a thermostat, then
        the contribution to the conserved quantity due to this subtraction is
        added to the thermostat heat energy, as it is assumed that the centre of
        mass motion is due to the thermostat.

        If there is a choice of thermostats, the thermostat
        connected to the centroid is chosen.
        """

        if (self.fixcom):
            na3 = self.beads.natoms * 3
            nb = self.beads.nbeads
            p = dstrip(self.beads.p)
            m = dstrip(self.beads.m3)[:, 0:na3:3]
            M = self.beads[0].M
            Mnb = M*nb

            dens = 0
            for i in range(3):
                pcom = p[:, i:na3:3].sum()
                dens += pcom**2
                pcom /= Mnb
                self.beads.p[:, i:na3:3] -= m * pcom

            self.ensemble.eens += dens * 0.5 / Mnb

        if len(self.fixatoms) > 0:
            for bp in self.beads.p:
                m = dstrip(self.beads.m)
                self.ensemble.eens += 0.5 * np.dot(bp[self.fixatoms * 3], bp[self.fixatoms * 3] / m[self.fixatoms])
                self.ensemble.eens += 0.5 * np.dot(bp[self.fixatoms * 3 + 1], bp[self.fixatoms * 3 + 1] / m[self.fixatoms])
                self.ensemble.eens += 0.5 * np.dot(bp[self.fixatoms * 3 + 2], bp[self.fixatoms * 3 + 2] / m[self.fixatoms])
                bp[self.fixatoms * 3] = 0.0
                bp[self.fixatoms * 3 + 1] = 0.0
                bp[self.fixatoms * 3 + 2] = 0.0



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

    def step_Oc(self):

        for i in xrange(self.nsteps_o):
            self.thermostat.step()
            self.proj_cotangent(self.beads)

    def step(self, step=None):
        """Does one simulation time step."""
        #print("~~~~~~~~~~~~~~~~~~~~~~ tau = ", self.thermostat.tau)
        #print("~~~~~~~~~~~~~~~~~~~~~~ dt = ", self.thermostat.dt)
        if self.splitting == "gobabo":
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

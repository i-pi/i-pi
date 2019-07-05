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
                constraint_list=[], tolerance=.0001, maxit=1000, norm_order=2):

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

        # MS: should all these be initialized as dependent arrays/values.?
        self.nsteps_geo = nsteps_geo
        self.nsteps_o = nsteps_o
        self.tolerance = tolerance
        self.maxit = maxit
        if norm_order == -1:
            self.norm_order = float('inf')
        else:
            self.norm_order = norm_order

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
        print("***************************", self.nmts)
        super(ConstrainedDynamics, self).bind(ens, beads, nm, cell, bforce, prng, omaker)

        """
        if len(self.nmts) > 1 or (len(self.nmts) == 1 and self.nmts[0] != 1):
            raise ValueError("MTS with constrains has not been implemented.")
        """

class ConstraintBase(dobject):
    """ Constraint class for MD. Base class."""

    def __init__(self, ncons):
        self.ncons = ncons

    def gfunc(self):
        raise NotImplementedError()

    def Dgfunc(self):
        """
        Calculates the Jacobian of the constraint.
        """
        raise NotImplementedError()

class RigidBondConstraint(ConstraintBase):
    """ Constraint class for MD.
        Specialized for rigid bonds. This can actually hold a *list*
        of rigid bonds, i.e. there will be a list of pairs of atoms and
        a list of bond lengths. """

    def __init__(self,constrained_indices,constrained_distances):

        super(RigidBondConstraint,self).__init__(ncons=len(constrained_distances))
        #print("@@@@@@@: ", self.ncons)
        #print("!!!!!!1: ", constrained_indices.shape)
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


class ConstraintList(dobject):
    """ Constraint class for MD"""

    def __init__(self, constraint_list):
        self.constraint_list = constraint_list
        self.ncons = sum([constr.ncons for constr in constraint_list])

    def gfunc(self, beads):
        """
        Compute the constraint function.
        """
        r = np.zeros(self.ncons)
        si = 0
        for constr in self.constraint_list:
            r[si:si+constr.ncons] = constr.gfunc(beads)
            si += constr.ncons
        return r

    def Dgfunc(self, beads):
        """
        Compute the Jacobian of the constraint function.
        """
        r = np.zeros((self.ncons, 3 * beads.natoms))
        si = 0
        for constr in self.constraint_list:
            r[si:si+constr.ncons,:] = constr.Dgfunc(beads)
            si += constr.ncons
        return r

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

        dd(self).tolerance = depend_value(name="tolerance", func=lambda: motion.tolerance)
        dd(self).maxit = depend_value(name="maxit", func=lambda: motion.maxit)
        dd(self).norm_order = depend_value(name="norm_order", func=lambda: motion.norm_order)

    #print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Number of constraints: ", self.constraints.ncons)

    def update_constraints(self, beads):
        self.Dg = self.constraints.Dgfunc(beads)
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

        if self.constraints.ncons > 0:

            #print(self.GramChol)
            b = np.dot(self.Dg, self.beads.p[0]/beads.m3[0])
            #print(b)
            x = np.linalg.solve(np.transpose(self.GramChol),np.linalg.solve(self.GramChol, b))
            beads.p[0] += - np.dot(np.transpose(self.Dg),x)


    def proj_manifold(self, beads, stepsize=None, proj_p=True):
        '''
        projects onto Manifold using the Gram matrix defined by self.Dg and self.Gram
        '''
        if proj_p and stepsize is None:
            stepsize = self.dt

        self.g = self.constraints.gfunc(beads)
        #print("beads.q[0] before",  beads.q[0])
        #print("g :", self.g)
        #print("Dg :",self.Dg)
        #print("Gram :",self.Gram)
        #print("GramChol :",self.GramChol)

        i = 0
        if self.constraints.ncons > 0:
            while (i < self.maxit and self.tolerance <= np.linalg.norm(self.g, ord=self.norm_order)):
                #print("proj_manifold index", i)
                #print("proj_manifold error", np.linalg.norm(self.g))
                #print(self.GramChol)
                #print("g, i :", self.g, i)

                self.g = self.constraints.gfunc(beads)
                dlambda = np.linalg.solve(np.transpose(self.GramChol),np.linalg.solve(self.GramChol, self.g))
                #print("dlambda", dlambda)
                delta = np.dot(np.transpose(self.Dg),dlambda)
                beads.q[0] += - delta /beads.m3[0]
                #self.d_momentum += - delta/self.p_stepsize
                self.g = self.constraints.gfunc(beads)
                #print("g, i after ud:", self.g, i)
                #print("beads.q[0] after",  beads.q[0])
                #print("delta",  delta)
                #print("stepsize", stepsize)
                #print("beads.m3[0]", beads.m3[0])
                if proj_p:
                    update_diff = - delta / stepsize
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
        #print "**************** NVEIntegrator_constraint init **************"
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

"""Contains the classes that deal with the different constrained dynamics
required in different types of ensembles.

Holds the algorithms for solving projecting phase-space coordinates onto
manifolds defined by the constraints, and the integrators to perform
constant energy and constant temperature dynamics.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np

from ipi.engine.motion import Dynamics
from ipi.engine.motion.dynamics import DummyIntegrator

from ipi.utils.depend import *
from ipi.utils.constrtools import *
from ipi.engine.thermostats import Thermostat
from ipi.engine.barostats import Barostat

from ipi.utils.messages import verbosity, warning


# tries to import scipy to do the cholesky decomposition solver,
# but falls back on numpy if it's not there
try:
    import scipy.linalg as spla
except ImportError:
    spla = None


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

    def __init__(
        self,
        timestep,
        mode="nve",
        splitting="obabo",
        thermostat=None,
        barostat=None,
        fixcom=False,
        fixatoms=None,
        nmts=None,
        nsteps_geo=1,
        nsteps_o=1,
        constraint_list=[],
        csolver=None,
    ):
        """Initialises a "dynamics" motion object.

        Args:
            dt: The timestep of the simulation algorithms.
            fixcom: An optional boolean which decides whether the centre of mass
                motion will be constrained or not. Defaults to False.
        """

        super(Dynamics, self).__init__(fixcom=fixcom, fixatoms=fixatoms)
        self._dt = depend_value(name="dt", value=timestep)

        if thermostat is None:
            self.thermostat = Thermostat()
        else:
            self.thermostat = thermostat

        if nmts is None or len(nmts) == 0:
            self._nmts = depend_array(name="nmts", value=np.asarray([1], int))
        else:
            self._nmts = depend_array(name="nmts", value=np.asarray(nmts, int))

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
            ValueError(
                "{:s} is not a valid ensemble for constrained dynamics".format(
                    self.enstype
                )
            )

        # splitting mode for the integrators
        self._splitting = depend_value(name="splitting", value=splitting)

        # The list of constraints coming from the input is an actual list of independent
        # constraints, and should be treated as such
        if constraint_list is not None:
            self.constraint_list = constraint_list

        has_eckart = False
        for constr in self.constraint_list:
            if isinstance(constr, ConstraintList):
                for c in constr.constraint_list:
                    if isinstance(c, ConstraintList):
                        raise ValueError("Cannot have nested constraint lists!")
                    elif isinstance(c, EckartConstraint):
                        has_eckart = True
            else:
                if isinstance(constr, EckartConstraint):
                    has_eckart = True

        self.fixcom = fixcom
        # If only some of the molecules have Eckart enforced, this will clash
        # with fixing the overall CoM; if *all* the molecules have Eckart
        # enforced, then cell CoM already fixed; in either case, raise
        # an error if fixcom is True and any of the constraints are Eckart.
        if self.fixcom and has_eckart:
            raise ValueError(
                "Cannot simultaneously fix cell CoM and enforce Eckart constraints!"
            )

        if fixatoms is None:
            self.fixatoms = np.zeros(0, int)
        else:
            self.fixatoms = fixatoms
        if len(self.fixatoms) > 0:
            raise ValueError("Cannot fix atoms together with constrained MD")

        if csolver is None:
            self.csolver = ConstraintSolver(
                self.constraint_list, tolerance=0.0001, maxit=10000, norm_order=2
            )
        else:
            if csolver["norm_order"] == -1:
                norm_order = float("inf")
            else:
                norm_order = csolver["norm_order"]
            self.csolver = ConstraintSolver(
                self.constraint_list,
                tolerance=csolver["tolerance"],
                maxit=csolver["maxit"],
                norm_order=norm_order,
            )

        # parameters of the geodesic integrator. will probably never need
        # to be changed dynamically, so for the moment we don't consider them depend objects
        self.nsteps_geo = nsteps_geo
        self.nsteps_o = nsteps_o

    def get_fixdof(self):
        """Calculate the number of fixed degrees of freedom, required for
        temperature and pressure calculations.
        """
        fixdof = super(ConstrainedDynamics, self).get_fixdof()
        for c in self.constraint_list:
            fixdof += c.ncons
        return fixdof

    def bind(self, ens, beads, nm, cell, bforce, prng, omaker):
        """Binds the constrained dynamics object.
        Just passes control to the parent class, and then binds the constraints
        and the solver.
        """
        super(ConstrainedDynamics, self).bind(
            ens, beads, nm, cell, bforce, prng, omaker
        )

        # now binds the constraints
        for c in self.constraint_list:
            c.bind(beads)
        self.csolver.bind(beads)


dproperties(ConstrainedDynamics, ["dt", "nmts", "splitting"])


class ConstraintSolverBase:
    """Empty base class for the constraint solver. Provides the interface
    that must be used to offer constraint functionalities to an integrator.
    """

    def __init__(self, constraint_list, dt=1.0):
        self.constraint_list = constraint_list

        # time step - will have to be linked to the dynamics time step
        self._dt = depend_value(name="dt", value=dt)

    def bind(self, beads):
        if beads.nbeads > 1:
            raise ValueError(
                "Constrained dynamics is only implemented for the case of classical MD (nbeads=1)"
            )

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


dproperties(ConstraintSolverBase, ["dt"])


class ConstraintSolver(ConstraintSolverBase):
    """An implementation of a constraint solver that uses M-RATTLE to
    impose constraints onto the momenta and a quasi-Newton method
    to impose constraints onto the positions. The constraint is applied
    sparsely, i.e. on each block of constraints separately.

    For implementation details of M-RATTLE see
        H. C. Andersen, J. Comput. Phys. 52, 24 (1983),
    and
        M. Tuckerman, "Statistical Mechanics: Theory and Molecular Simulation" pp 106-108 (OUP, Oxford, 2010).

    For the quasi-Newton method in proj_cotangent see
        ???
    """

    def __init__(
        self, constraint_list, dt=1.0, tolerance=0.001, maxit=1000, norm_order=2
    ):
        """Solver options include a tolerance for the projection on the
        manifold, maximum number of iterations for the projection, and
        the order of the norm to estimate the convergence"""

        super(ConstraintSolver, self).__init__(constraint_list, dt)
        self.tolerance = tolerance
        self.maxit = maxit
        self.norm_order = norm_order

    def proj_cotangent(self):
        """Set the momenta conjugate to the constrained degrees of freedom
        to zero non-iteratively by inverting the Gram matrix. For further
        details see H. C. Andersen, J. Comput. Phys. 52, 24 (1983). Note
        that independent groups of constraints are treated separately
        (sparsely).
        """

        p = dstrip(self.beads.p[0]).copy()

        for constr in self.constraint_list:
            dg = dstrip(constr.Dg)
            ic = constr.i3_unique
            b = np.dot(dg, p[ic] / constr.m3)

            if spla is None:
                x = np.linalg.solve(dstrip(constr.Gram), b)
            else:
                x = spla.cho_solve(dstrip(constr.GramChol), b)

            p[ic] -= np.dot(np.transpose(dg), x)
        self.beads.p[0] = p

    def proj_manifold(self):
        """Iteratively enforce the constraints onto the positions by finding
        the Lagrange multipliers using a quasi-Newton solver. Note
        that independent groups of constraints are treated separately
        (sparsely).
        """

        m3 = dstrip(self.beads.m3[0])
        p = dstrip(self.beads.p[0]).copy()
        q = dstrip(self.beads.q[0]).copy()

        for constr in self.constraint_list:
            dg = dstrip(constr.Dg)

            if spla is None:
                gram = dstrip(constr.Gram)
            else:
                chol_gram = dstrip(constr.GramChol)

            ic = constr.i3_unique
            constr.q = q[ic]

            # iterative projection on the manifold
            for i in range(self.maxit):
                g = dstrip(constr.g)
                # bailout condition
                if self.tolerance > np.linalg.norm(g, ord=self.norm_order):
                    break

                if spla is None:
                    dlambda = np.linalg.solve(gram, g)
                else:
                    dlambda = spla.cho_solve(chol_gram, g)
                delta = np.dot(np.transpose(dg), dlambda)
                q[ic] -= delta / m3[ic]
                constr.q = q[ic]  # updates the constraint to recompute g

                p[ic] -= delta / self.dt

            if i == self.maxit:
                warning(
                    "No convergence in Newton iteration for positional component",
                    verbosity.low,
                )

        # after all constraints have been applied, q is on the manifold and we can update the constraint positions
        for constr in self.constraint_list:
            constr.qprev = q[constr.i3_unique.flatten()]
        self.beads.p[0] = p
        self.beads.q[0] = q


class ConstrainedIntegrator(DummyIntegrator):
    """No-op integrator for classical constrained propagation.
    It also incorporates a constraint solver, so as to make the
    integration modular in case one wanted to implement multiple
    solvers.

    Note that this and other constrained integrators inherit
    from the non-constrained dynamics class, so the logic of
    the integration and multiple time stepping machinery is
    better documented in dynamics.py
    """

    def __init__(self):
        super(ConstrainedIntegrator, self).__init__()

    def get_qdt(self):
        # get the base dt for doing q propagation (center of the integrator)
        return super(ConstrainedIntegrator, self).get_qdt() / self.nsteps_geo

    def get_tdt(self):
        return super(ConstrainedIntegrator, self).get_tdt() / self.nsteps_o

    def bind(self, motion):
        """Creates local references to all the variables for simpler access."""

        super(ConstrainedIntegrator, self).bind(motion)

        self.constraint_list = motion.constraint_list
        if motion.nsteps_geo is None:
            self._nsteps_geo = depend_value(name="nsteps_geo", value=1)
        else:
            self._nsteps_geo = depend_value(name="nsteps_geo", value=motion.nsteps_geo)

        self._qdt.add_dependency(self._nsteps_geo)
        if motion.nsteps_o is None:
            self._nsteps_o = depend_value(name="nsteps_o", value=1)
        else:
            self._nsteps_o = depend_value(name="nsteps_o", value=motion.nsteps_o)

        self._tdt.add_dependency(self._nsteps_o)
        self.csolver = motion.csolver
        dpipe(self._qdt, self.csolver._dt)

    def proj_cotangent(self):
        self.csolver.proj_cotangent()

    def proj_manifold(self):
        self.csolver.proj_manifold()


dproperties(ConstrainedIntegrator, ["nsteps_geo", "nsteps_o"])


class NVEConstrainedIntegrator(ConstrainedIntegrator):
    """A propagator for constant-energy integration under a set of constraints.

    Implementation details:
        B. Leimkuhler, C. Matthews Proc. R. Soc. A 472, 20160138, (2016)
    """

    def step_A(self):
        """Unconstrained A-step (coordinate integration)"""
        self.beads.q[0] += self.beads.p[0] / dstrip(self.beads.m3)[0] * self.qdt

    def step_B(self, level=0):
        """Unconstrained B-step (momentum integration)"""
        self.beads.p[0] += self.forces.forces_mts(level)[0] * self.pdt[level]

    def step_Bc(self, level=0):
        """Unconstrained B-step (momentum integration)
        followed by a projection into the cotangent space"""
        self.step_B(level)
        self.proj_cotangent()

    def step_Ag(self):
        """
        Geodesic flow integrator. Makes one A step including manifold
        projection of the position and update of the momenta following
        the orthogonalization. Can be broken down into a MTS-like fashion
        to increase the accuracy without having to recompute forces.
        """
        # Resolve momentum constraint and update Gram matrix if neccesary
        # GT: is the following line necessary?
        self.proj_cotangent()

        for i in range(self.nsteps_geo):
            self.step_A()
            self.proj_manifold()
            self.proj_cotangent()

    def mtsprop_ba(self, index):
        """Recursive MTS step -- this is adapted directly from the
        NVEIntegrator class"""

        mk = int(self.nmts[index] / 2)

        for i in range(mk):  # do nmts/2 full sub-steps
            self.step_Bc(index)
            self.pconstraints()  # currently does nothing
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
        """Recursive MTS step -- this is adapted directly from the
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

    Implementation details:
        B. Leimkuhler, C. Matthews Proc. R. Soc. A 472, 20160138, (2016)

    Attributes:
        thermostat: A thermostat object to keep the temperature constant.
    """

    def step_Oc(self):
        """Constrained stochastic propagation. We solve the problem that
        the thermostat and the projective step do not necessarily commute
        (e.g. in GLE) by doing a MTS splitting scheme"""

        m3 = dstrip(self.beads.m3)
        for i in range(self.nsteps_o):
            self.thermostat.step()

            # accumulates conserved quantity
            p = dstrip(self.beads.p)
            self.ensemble.eens += np.dot(p.flatten(), (p / m3).flatten()) * 0.5

            self.proj_cotangent()

            p = dstrip(self.beads.p).copy()
            self.ensemble.eens -= np.dot(p.flatten(), (p / m3).flatten()) * 0.5

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

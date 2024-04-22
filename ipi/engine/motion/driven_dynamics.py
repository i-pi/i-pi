"""Contains the classes that deal with the different dynamics required in
different types of ensembles.

Holds the algorithms required for normal mode propagators, and the objects to
do the constant temperature and pressure algorithms. Also calculates the
appropriate conserved energy quantity for the ensemble of choice.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

import numpy as np

from ipi.utils.depend import *
from ipi.engine.motion.eda import EDA
from ipi.utils.units import Constants
from ipi.engine.motion.dynamics import NVEIntegrator, DummyIntegrator, Dynamics


class DrivenDynamics(Dynamics):
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

    def __init__(
        self,
        efield=None,
        bec=None,
        *argc,
        **argv,
    ):
        """Initialises a "dynamics" motion object.

        Args:
            dt: The timestep of the simulation algorithms.
            fixcom: An optional boolean which decides whether the centre of mass
                motion will be constrained or not. Defaults to False.
        """

        super().__init__(*argc, **argv)

        if self.enstype == "eda-nve":
            # NVE integrator with an external time-dependent driving (electric field)
            self.integrator = EDANVEIntegrator()
        else:
            self.integrator = DummyIntegrator()

        # if the dynamics is driven, allocate necessary objects
        self.efield = efield
        self.bec = bec
        self.eda = EDA(self.efield, self.bec)

    def bind(self, ens, beads, nm, cell, bforce, prng, omaker):

        super().bind(ens, beads, nm, cell, bforce, prng, omaker)

        self.eda.bind(self.ensemble, self.enstype)

        # now that the timesteps are decided, we proceed to bind the integrator.
        self.integrator.bind(self)

        # applies constraints immediately after initialization.
        self.integrator.pconstraints()

    def step(self, step=None):
        """Advances the dynamics by one time step"""

        super().step(step)
        # self.integrator.step(step)
        # self.ensemble.time += self.dt  # increments internal time

        # Check that these variable are the same.
        # If they are not the same, then there is a bug in the code
        dt = abs(self.ensemble.time - self.integrator.efield_time)
        if dt > 1e-12:
            raise ValueError(
                "The time at which the Electric Field is evaluated is not properly updated!"
            )


dproperties(DrivenDynamics, ["dt", "nmts", "splitting", "ntemp"])


class EDAIntegrator(DummyIntegrator):
    """Integrator object for simulations using the Electric Dipole Approximation (EDA)
    when an external electric field is applied.
    """

    def bind(self, motion):
        """bind variables"""
        super().bind(motion)

        self._time = self.eda._time
        self._mts_time = self.eda._mts_time

        dep = [
            self._time,
            self._mts_time,
            self.eda.Born_Charges._bec,
            self.eda.Electric_Field._Efield,
        ]
        self._EDAforces = depend_array(
            name="EDAforces",
            func=self._eda_forces,
            value=np.zeros((self.beads.nbeads, self.beads.natoms * 3)),
            dependencies=dep,
        )
        pass

    def pstep(self, level=0):
        """Velocity Verlet momentum propagator."""
        self.beads.p += self.EDAforces * self.pdt[level]
        if dstrip(self.efield_time) == dstrip(self.time):
            # it's the first time that 'pstep' is called
            # then we need to update 'efield_time'
            self.efield_time += dstrip(self.dt)
            # the next time this condition will be 'False'
            # so we will avoid to re-compute the EDAforces
        pass

    def _eda_forces(self):
        """Compute the EDA contribution to the forces, i.e. `q_e Z^* @ E(t)`"""
        Z = dstrip(self.eda.Born_Charges.bec)  # tensor of shape (nbeads,3xNatoms,3)
        E = dstrip(self.eda.Electric_Field.Efield)  # vector of shape (3)
        forces = Constants.e * Z @ E  # array of shape (nbeads,3xNatoms)
        return forces

    def step(self, step=None):
        if len(self.nmts) > 1:
            raise ValueError(
                "EDAIntegrator is not implemented with the Multiple Time Step algorithm (yet)."
            )
        # This should call 'NVEIntegrator.step' since 'self' should be an instance of 'EDANVEIntegrator' and not of 'EDAIntegrator'
        super().step(step)


dproperties(EDAIntegrator, ["EDAforces", "efield_time", "time"])


class EDANVEIntegrator(EDAIntegrator, NVEIntegrator):
    """Integrator object for simulations with constant Number of particles, Volume, and Energy (NVE)
    using the Electric Dipole Approximation (EDA) when an external electric field is applied.
    """

    def pstep(self, level):
        # NVEIntegrator does not use 'super()' within 'pstep'
        # then we can not use 'super()' here.
        # We need to call the 'pstep' methods explicitly.
        NVEIntegrator.pstep(
            self, level
        )  # the driver is called here: add nuclear and electronic forces (DFT)
        EDAIntegrator.pstep(self, level)  # add the driving forces, i.e. q_e Z @ E(t)
        pass

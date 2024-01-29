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
from ipi.utils.units import Constants
from ipi.engine.eda import EDA
from ipi.engine.motion.dynamics import NVEIntegrator, DummyIntegrator

__all__ = ["EDANVEIntegrator"]


class EDAIntegrator(DummyIntegrator):
    """Integrator object for simulations using the Electric Dipole Approximation (EDA)
    when an external electric field is applied.
    """

    # author: Elia Stocco
    # motivation: deal with time-dependent external potential

    def __init__(self):
        super().__init__()

    def bind(self, motion):
        """bind variables"""
        super().bind(motion)

        dself = dd(self)

        dself.time = dd(self.eda).time
        dself.mts_time = dd(self.eda).mts_time

        dep = [
            dself.time,
            dself.mts_time,
            dd(self.eda).bec,
            dd(self.eda).Efield,
        ]
        dself.EDAforces = depend_array(
            name="forces",
            func=self._eda_forces,
            value=np.zeros((self.beads.nbeads, self.beads.natoms * 3)),
            dependencies=dep,
        )
        pass

    def pstep(self, level=0):
        """Velocity Verlet momentum propagator."""
        self.beads.p += self.EDAforces * self.pdt[level]
        if dstrip(self.mts_time) == dstrip(self.time):
            # it's the first time that 'pstep' is called
            # then we need to update 'mts_time'
            self.mts_time += dstrip(self.dt)
            # the next time this condition will be 'False'
            # so we will avoid to re-compute the EDAforces
        pass

    def _eda_forces(self):
        """Compute the EDA contribution to the forces due to the polarization"""
        Z = dstrip(self.eda.bec)
        E = dstrip(self.eda.Efield)
        forces = Constants.e * Z @ E
        return forces

    def step(self, step=None):
        if len(self.nmts) > 1:
            raise ValueError(
                "EDAIntegrator is not implemented with the Multiple Time Step algorithm (yet)."
            )
        super().step(step)


class EDANVEIntegrator(EDAIntegrator, NVEIntegrator):
    """Integrator object for simulations with constant Number of particles, Volume, and Energy (NVE)
    using the Electric Dipole Approximation (EDA) when an external electric field is applied.
    """

    # author: Elia Stocco

    def pstep(self, level):
        # NVEIntegrator does not use 'super()' within 'pstep'
        # then we can not use 'super()' here.
        # We need to call the 'pstep' methods explicitly.
        NVEIntegrator.pstep(
            self, level
        )  # the driver is called here: add nuclear and electronic forces (DFT)
        EDAIntegrator.pstep(self, level)  # add the driving forces, i.e. q_e Z @ E(t)
        pass

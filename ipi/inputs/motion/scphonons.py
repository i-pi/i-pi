"""Deals with creating the ensembles class.

Copyright (C) 2013, Joshua More and Michele Ceriotti

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http.//www.gnu.org/licenses/>.

Inputs created by Michele Ceriotti and Benjamin Helfrecht, 2015

Classes:
   InputGeop: Deals with creating the Geop object from a file, and
      writing the checkpoints.
"""

import numpy as np

# import ipi.engine.initializer
from ipi.engine.motion import *
from ipi.utils.inputvalue import *
from ipi.inputs.thermostats import *
from ipi.inputs.initializer import *
from ipi.utils.units import *

__all__ = ["InputSCPhonons"]


class InputSCPhonons(InputDictionary):
    """Calculation options for self-consistent phonons algorithm.

    Contains options related to self consistent phonons method.

    """

    attribs = {
        "mode": (
            InputAttribute,
            {
                "dtype": str,
                "default": "qn",
                "help": "The statistics to be used in the calculation of the free energy. Quantum (qn) or classical (cl) Boltzmann statistics.",
                "options": ["qn", "cl"],
            },
        )
    }
    fields = {
        "prefix": (
            InputValue,
            {"dtype": str, "default": "", "help": "Prefix of the output files."},
        ),
        "asr": (
            InputValue,
            {
                "dtype": str,
                "default": "none",
                "options": ["none", "crystal", "poly"],
                "help": "The method used to project out zero modes coming from continuous symmetries: crystal removes the three translational modes; molecule removes the three rotational modes in addition to the translational ones. none keeps all the modes.",
            },
        ),
        "random_type": (
            InputValue,
            {
                "dtype": str,
                "default": "pseudo",
                "options": ["sobol", "pseudo", "file"],
                "help": "Chooses the type of random numbers.",
            },
        ),
        "displace_mode": (
            InputValue,
            {
                "dtype": str,
                "default": "nmik",
                "options": ["ik", "sd", "nmik", "rnmik"],
                "help": "The type of optimisation strategy for obtaining the mean position. sd stands for a steepest descent algorithm. ik stands for a Newton-Raphson scheme that requires the inverse of the force constant matrix iK. nmik stands for a Newton-Raphson scheme that only displaces along normal modes directions with statistically significant forces. rnmik same as nmik but performs several optimization steps using a reweighted sampling.",
            },
        ),
        "dynmat": (
            InputArray,
            {
                "dtype": float,
                "default": np.zeros(0, float),
                "help": "The dynamical matrix of the trial Hamiltonian.",
            },
        ),
        "max_steps": (
            InputValue,
            {
                "dtype": int,
                "default": None,
                "help": "Maximum number of Monte carlo steps per SCP iteration.",
            },
        ),
        "max_iter": (
            InputValue,
            {"dtype": int, "default": 1, "help": "Maximum number of SCP iterations."},
        ),
        "tau": (
            InputValue,
            {
                "dtype": float,
                "default": 1.0,
                "help": "Step size along the gradient for the sd displace_mode",
            },
        ),
        "wthreshold": (
            InputValue,
            {
                "dtype": float,
                "default": 0.90,
                "help": "Threshold on minimum Boltzmann weights before more statistics must be accumulated.",
            },
        ),
        "precheck": (
            InputValue,
            {
                "dtype": bool,
                "default": True,
                "help": "Flag for checking statistical significance of forces before optimisation of mean position.",
            },
        ),
        "checkweights": (
            InputValue,
            {
                "dtype": bool,
                "default": True,
                "help": "Flag for checking Boltzmann weights for whether more statistics are required.",
            },
        ),
        "chop": (
            InputValue,
            {
                "dtype": float,
                "default": 1e-09,
                "help": "Threshold below which frequencies are set to zero.",
            },
        ),
        "nparallel": (
            InputValue,
            {
                "dtype": int,
                "default": 1,
                "help": "The number of Monte Carlo forces to be evaluated (in parallel) per i-PI step.",
            },
        ),
        "batch_weight_exponent": (
            InputValue,
            {
                "dtype": int,
                "default": 1,
                "help": "The exponent used to suppress low batch weights.",
            },
        ),
    }

    dynamic = {}

    default_help = "Self-consistent phonons class. It variationally optimizes the free energy to calculate the best harmonic approximation to a system."
    default_label = "SCPHONONS"

    def store(self, phonons):
        if phonons == {}:
            return
        self.mode.store(phonons.mode)
        self.prefix.store(phonons.prefix)
        self.asr.store(phonons.asr)
        self.dynmat.store(phonons.dynmatrix)
        self.max_steps.store(phonons.max_steps)
        self.max_iter.store(phonons.max_iter)
        self.tau.store(phonons.tau)
        self.chop.store(phonons.chop)
        self.random_type.store(phonons.random_type)
        self.displace_mode.store(phonons.displace_mode)
        self.nparallel.store(phonons.nparallel)
        self.batch_weight_exponent.store(phonons.batch_weight_exponent)

    def fetch(self):
        rv = super(InputSCPhonons, self).fetch()
        rv["mode"] = self.mode.fetch()
        return rv

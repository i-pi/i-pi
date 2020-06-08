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


Classes:
   InputEnsemble: Deals with creating the Ensemble object from a file, and
      writing the checkpoints.
"""

import numpy as np
from ipi.engine.motion import *
from ipi.utils.inputvalue import *
from ipi.inputs.thermostats import *
from ipi.inputs.initializer import *
from ipi.utils.units import *

__all__ = ["InputNEB"]


class InputNEB(InputDictionary):

    """Geometry optimization options for nudged elastic band (NEB) calculations.

    Contains options related with geometry optimization, such as method,
    thresholds, linear search strategy, etc.

    Also contains options related specifically to NEB, such as spring constants
    and climbing image.

    """

    attribs = {
        "mode": (
            InputAttribute,
            {
                "dtype": str,
                "default": "lbfgs",
                "help": "The geometry optimization algorithm to be used",
                "options": ["sd", "cg", "bfgs", "lbfgs"],
            },
        )
    }

    fields = {
        "ls_options": (
            InputDictionary,
            {
                "dtype": [float, int, float, float],
                "help": """Options for line search methods. Includes:
                              tolerance: stopping tolerance for the search,
                              grad_tolerance: stopping tolerance on gradient for
                              BFGS line search,
                              iter: the maximum number of iterations,
                              step: initial step for bracketing,
                              adaptive: whether to update initial step.
                              """,
                "options": ["tolerance", "iter", "step", "adaptive"],
                "default": [1e-6, 100, 1e-3, 1.0],
                "dimension": ["energy", "undefined", "length", "undefined"],
            },
        ),
        "tolerances": (
            InputDictionary,
            {
                "dtype": float,
                "options": ["energy", "force", "position"],
                "default": [1e-8, 1e-8, 1e-8],
                "dimension": ["energy", "force", "length"],
            },
        ),
        "old_force": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "The previous force in an optimization step.",
                "dimension": "force",
            },
        ),
        "old_direction": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "The previous direction.",
            },
        ),
        "biggest_step": (
            InputValue,
            {
                "dtype": float,
                "default": 100.0,
                "help": "The maximum step size for (L)-BFGS line minimizations.",
            },
        ),
        "scale_lbfgs": (
            InputValue,
            {
                "dtype": int,
                "default": 2,
                "help": """Scale choice for the initial hessian.
                                            0 identity.
                                            1 Use first member of position/gradient list.
                                            2 Use last  member of position/gradient list.""",
            },
        ),
        "invhessian_bfgs": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.eye, args=(0,)),
                "help": "Approximate inverse Hessian for BFGS, if known.",
            },
        ),
        "qlist_lbfgs": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "List of previous position differences for L-BFGS, if known.",
            },
        ),
        "glist_lbfgs": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "List of previous gradient differences for L-BFGS, if known.",
            },
        ),
        "corrections_lbfgs": (
            InputValue,
            {
                "dtype": int,
                "default": 5,
                "help": "The number of past vectors to store for L-BFGS.",
            },
        ),
        "endpoints": (
            InputDictionary,
            {
                "dtype": [bool, str],
                "options": ["optimize", "algorithm"],
                "default": [True, "bfgs"],
                "help": "Geometry optimization of endpoints",
            },
        ),
        "spring": (
            InputDictionary,
            {
                "dtype": [bool, float, float, float],
                "options": ["varsprings", "kappa", "kappamax", "kappamin"],
                "default": [False, 1.0, 1.5, 0.5],
                "help": "Uniform or variable spring constants along the elastic band",
            },
        ),
        "climb": (
            InputValue,
            {"dtype": bool, "default": False, "help": "Use climbing image NEB"},
        ),
    }

    dynamic = {}

    default_help = "Contains the required parameters for performing nudged elastic band (NEB) calculations"
    default_label = "NEB"

    def store(self, neb):
        if neb == {}:
            return
        self.ls_options.store(neb.ls_options)
        self.tolerances.store(neb.tolerances)
        self.mode.store(neb.mode)
        self.old_force.store(neb.old_f)
        self.old_direction.store(neb.old_d)
        self.biggest_step.store(neb.big_step)
        self.invhessian_bfgs.store(neb.invhessian)
        self.qlist_lbfgs.store(neb.qlist)
        self.glist_lbfgs.store(neb.glist)
        self.endpoints.store(neb.endpoints)
        self.spring.store(neb.spring)
        self.climb.store(neb.climb)
        self.scale_lbfgs.store(neb.scale)

    def fetch(self):
        rv = super(InputNEB, self).fetch()
        rv["mode"] = self.mode.fetch()
        return rv

"""Deals with creating the string class.

Copyright (C) 2022, Karen Fidanyan

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

__all__ = ["InputStringMEP"]


class InputStringMEP(InputDictionary):
    """Options for string minimal energy path calculations.

    Contains options related to geometry optimization, such as method,
    thresholds, etc.
    Also contains options related specifically to String method.
    """

    attribs = {
        "mode": (
            InputAttribute,
            {
                "dtype": str,
                "default": "bfgstrm",
                "help": "The geometry optimization algorithm to optimize MEP string",
                "options": [
                    "sd",
                    "cg",
                    "bfgs",
                    "bfgstrm",
                    "damped_bfgs",
                    "lbfgs",
                    "fire",
                    "euler",
                ],
            },
        )
    }

    fields = {
        "tolerances": (
            InputDictionary,
            {
                "dtype": float,
                "options": ["energy", "force", "position"],
                "default": [1e-8, 1e-8, 1e-8],
                "dimension": ["energy", "force", "length"],
                "help": """Tolerance criteria to stop String optimization.
                           If you work with DFT, do not use these defaults.
                        """,
            },
        ),
        "old_coord": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "The previous position in an optimization step.",
                "dimension": "length",
            },
        ),
        "full_force": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "The previous full-dimensional force in an optimization step.",
                "dimension": "force",
            },
        ),
        "full_pots": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "Previous physical potentials of all beads.",
                "dimension": "energy",
            },
        ),
        "old_stringpotential": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "Previous string potential energy.",
            },
        ),
        "old_stringgradient": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "The previous gradient of the string.",
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
                "default": 0.5,
                "help": """The maximum atomic displacement in a single step
                           of optimizations within String MEP procedure.
                           If requested step is larger, it will be downscaled so
                           that maximal atomic displacement won't exceed biggest_step.
                        """,
                "dimension": "length",
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
        "hessian_bfgs": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.eye, args=(0,)),
                "help": "Approximate Hessian for damped_BFGS, if known.",
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
        "tr_trm": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.array, args=(1.0,)),
                "help": "Starting value for the trust radius for BFGSTRM.",
                "dimension": "length",
            },
        ),
        "dtmax_fire": (
            InputValue,
            {
                "dtype": float,
                "default": 1.0,
                "help": "Maximum time interval per step for FIRE.",
            },
        ),
        "v_fire": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "Current velocity for FIRE",
            },
        ),
        "alpha_fire": (
            InputValue,
            {
                "dtype": float,
                "default": 0.1,
                "help": "velocity mixing factor for FIRE",
            },
        ),
        "N_down_fire": (
            InputValue,
            {
                "dtype": int,
                "default": 0,
                "help": "consecutive steps in downhill dierction for FIRE",
            },
        ),
        "N_up_fire": (
            InputValue,
            {
                "dtype": int,
                "default": 0,
                "help": "consecutive steps in uphill direction",
            },
        ),
        "dt_fire": (
            InputValue,
            {
                "dtype": float,
                "default": 0.1,
                "help": "time per step",
            },
        ),
        "endpoints": (
            InputDictionary,
            {
                "dtype": [bool, str],
                "options": ["optimize", "algorithm"],
                "default": [False, "bfgstrm"],
                "help": "Geometry optimization of endpoints (not implemented yet)",
            },
        ),
        "stage": (
            InputValue,
            {
                "dtype": str,
                "options": ["endpoints", "string", "climb"],
                "default": "string",
                "help": "Stage of the String pipeline: optimization of the endpoints, "
                "string opt., climbing image opt.",
            },
        ),
        "use_climb": (
            InputValue,
            {
                "dtype": bool,
                "default": False,
                "help": "Use climbing image String MEP or not",
            },
        ),
        "climb_bead": (
            InputValue,
            {
                "dtype": int,
                "default": -1,
                "help": "The index of the climbing bead.",
            },
        ),
    }

    dynamic = {}

    default_help = (
        "Contains the required parameters "
        "for performing string minimal energy path optimization."
    )
    default_label = "StringMEP"

    def store(self, stringmep):
        if stringmep == {}:
            return
        self.tolerances.store(stringmep.tolerances)
        self.mode.store(stringmep.mode)
        self.old_coord.store(stringmep.old_x)
        self.full_force.store(stringmep.full_f)
        self.full_pots.store(stringmep.full_v)
        self.old_stringpotential.store(stringmep.stringpot)
        self.old_stringgradient.store(stringmep.stringgrad)
        self.old_direction.store(stringmep.d)
        self.biggest_step.store(stringmep.big_step)
        self.hessian_bfgs.store(stringmep.hessian)
        self.qlist_lbfgs.store(stringmep.qlist)
        self.glist_lbfgs.store(stringmep.glist)
        self.tr_trm.store(stringmep.tr_trm)
        self.v_fire.store(stringmep.v)
        self.alpha_fire.store(stringmep.a)
        self.N_down_fire.store(stringmep.N_dn)
        self.N_up_fire.store(stringmep.N_up)
        self.dt_fire.store(stringmep.dt_fire)
        self.dtmax_fire.store(stringmep.dtmax)
        self.endpoints.store(stringmep.endpoints)
        self.stage.store(stringmep.stage)
        self.use_climb.store(stringmep.use_climb)
        self.climb_bead.store(stringmep.cl_indx)
        self.scale_lbfgs.store(stringmep.scale)

    def fetch(self):
        rv = super(InputStringMEP, self).fetch()
        rv["mode"] = self.mode.fetch()
        return rv

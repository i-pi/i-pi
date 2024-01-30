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
                "default": "fire",
                "help": "The geometry optimization algorithm to optimize NEB path",
                "options": [
                    # "sd",
                    # "cg",
                    # "bfgs",
                    "bfgstrm",
                    "damped_bfgs",
                    # "lbfgs",
                    "fire",
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
                "help": """Tolerance criteria to stop NEB optimization.
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
        "old_nebpotential": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "Previous NEB potential energy, which includes spring energy.",
            },
        ),
        "old_nebgradient": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "The previous gradient including NEB spring forces.",
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
                           of optimizations within NEB procedure.
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
                "default": [False, "bfgs"],
                "help": "Geometry optimization of endpoints (not implemented yet)",
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
        "tangent": (
            InputValue,
            {
                "dtype": str,
                "options": ["plain", "improved"],
                "default": "improved",
                "help": "How to calculate tangents: simple averaging from the original 1998 paper, "
                "or the improved tangent estimate from J. Chem. Phys. 113, 9978 (2000)",
            },
        ),
        "stage": (
            InputValue,
            {
                "dtype": str,
                "options": ["endpoints", "neb", "climb"],
                "default": "neb",
                "help": "Stage of the NEB pipeline: optimization of endpoints, "
                "NEB itself, climbing image",
            },
        ),
        "use_climb": (
            InputValue,
            {"dtype": bool, "default": False, "help": "Use climbing image NEB or not"},
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
        "for performing nudged elastic band (NEB) calculations"
    )
    default_label = "NEB"

    def store(self, neb):
        if neb == {}:
            return
        self.tolerances.store(neb.tolerances)
        self.mode.store(neb.mode)
        self.old_coord.store(neb.old_x)
        self.full_force.store(neb.full_f)
        self.full_pots.store(neb.full_v)
        self.old_nebpotential.store(neb.nebpot)
        self.old_nebgradient.store(neb.nebgrad)
        self.old_direction.store(neb.d)
        self.biggest_step.store(neb.big_step)
        self.hessian_bfgs.store(neb.hessian)
        self.qlist_lbfgs.store(neb.qlist)
        self.glist_lbfgs.store(neb.glist)
        self.v_fire.store(neb.v)
        self.alpha_fire.store(neb.a)
        self.N_down_fire.store(neb.N_dn)
        self.N_up_fire.store(neb.N_up)
        self.dt_fire.store(neb.dt_fire)
        self.dtmax_fire.store(neb.dtmax)
        self.endpoints.store(neb.endpoints)
        self.spring.store(neb.spring)
        self.tangent.store(neb.tangent)
        self.stage.store(neb.stage)
        self.use_climb.store(neb.use_climb)
        self.climb_bead.store(neb.cl_indx)
        self.scale_lbfgs.store(neb.scale)

    def fetch(self):
        rv = super(InputNEB, self).fetch()
        rv["mode"] = self.mode.fetch()
        return rv

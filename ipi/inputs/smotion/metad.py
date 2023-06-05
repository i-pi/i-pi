"""Deals with creating the Metadynamics class

Copyright (C) i-PI developers team.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.


Classes:
   InputMetaDyn: Deals with creating the Ensemble object from a file, and
      writing the checkpoints.

"""

import numpy as np
from ipi.engine.motion import *
from ipi.utils.inputvalue import *
from ipi.utils.units import *

__all__ = ["InputMetaDyn"]


class InputMetaDyn(InputDictionary):
    """Metadynamics Options

    Contains options related with MetaDynamics

    """

    fields = {
        "metaff": (
            InputArray,
            {
                "dtype": str,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "List of names of forcefields that should do metadynamics.",
            },
        ),
        "use_energy": (
            InputValue,
            {
                "dtype": bool,
                "default": False,
                "help": """ 
Transfer the potential energy value to PLUMED to use as a collective variable. 
Can only be used with classical simulations because it requires a rather hacky 
mechanism to transfer the energy of the system to the forcefield.
""",
            },
        ),
    }

    default_help = "MetaDynamics"
    default_label = "META"

    def store(self, meta):
        if meta == {}:
            return
        self.metaff.store(meta.metaff)
        self.use_energy.store(meta.use_energy)

    def fetch(self):
        rv = super(InputMetaDyn, self).fetch()
        return rv

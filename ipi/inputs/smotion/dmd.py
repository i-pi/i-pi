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
from ipi.utils.units import *

__all__ = ["InputDMD"]


class InputDMD(InputDictionary):
    """Driven MD Options

    Contains options related with a driven dynamics through a (possibly time-dependent) potential.


    """

    fields = {
        "dmdff": (
            InputArray,
            {
                "dtype": str,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "List of names of forcefields that should do driven MD. Accepts ffdmd forcefield types. Currently implemented (2021) for ffdmd only the driving potential similar to the one described in Bowman, .., Brown JCP 119, 646 (2003).",
            },
        )
    }

    default_help = "DrivenMD with external time-dependent driving potential"
    default_label = "DMD"

    def store(self, drivenmd):
        if drivenmd == {}:
            return
        self.dmdff.store(drivenmd.dmdff)

    def fetch(self):
        rv = super(InputDMD, self).fetch()
        return rv

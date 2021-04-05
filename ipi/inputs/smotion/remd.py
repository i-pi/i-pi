"""Deals with creating the Replica Exchange class

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
   InputReplicaExchange: Deals with creating the Ensemble object from a file, and
      writing the checkpoints.
"""

import numpy as np
from ipi.engine.motion import *
from ipi.utils.inputvalue import *
from ipi.utils.units import *

__all__ = ["InputReplicaExchange"]


class InputReplicaExchange(InputDictionary):
    """Replica Exchange options.

    Contains options related with replica exchange, such as method,
    steps on which REMD should be performed, etc.

    """

    fields = {
        "stride": (
            InputValue,
            {
                "dtype": float,
                "default": 1.0,
                "help": "Every how often to try exchanges (on average).",
            },
        ),
        "krescale": (
            InputValue,
            {
                "dtype": bool,
                "default": True,
                "help": "Rescale kinetic energy upon exchanges.",
            },
        ),
        "swapfile": (
            InputValue,
            {
                "dtype": str,
                "default": "remd_idx",
                "help": "File to keep track of replica exchanges",
            },
        ),
        "repindex": (
            InputArray,
            {
                "dtype": int,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "List of current indices of the replicas compared to the starting indices",
            },
        ),
    }

    default_help = "Replica Exchange"
    default_label = "REMD"

    def store(self, remd):
        if remd == {}:
            return
        self.stride.store(remd.stride)
        self.repindex.store(remd.repindex)
        self.krescale.store(remd.rescalekin)
        self.swapfile.store(remd.swapfile)

    def fetch(self):
        rv = super(InputReplicaExchange, self).fetch()
        return rv

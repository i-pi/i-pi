"""Creates objects that deal with the alchemical exchanges."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

import numpy as np
from ipi.utils.inputvalue import InputDictionary, InputAttribute, InputValue, InputArray, input_default


__all__ = ['InputAtomSwap']


class InputAtomSwap(InputDictionary):

    """Alchemy input class.

    Handles generating the appropriate alchemical exchange class from the xml input file.

    Fields:
        spicesA/B: isotopes for exchanges.
        timestep: An optional float giving the size of the timestep in atomic
            units. Defaults to 1.0.
    """
    attribs = {"mode": (InputAttribute, {"dtype": str,
                                         "default": 'dummy',
                                         "help": " ",
                                         "options": ['dummy']})}

    fields = {
        "names": (InputArray, {"dtype": str,
                               "default": input_default(factory=np.zeros, args=(0,), kwargs={'dtype': np.dtype('|S6')}),
                               "help": "The names of the atoms to be to exchanged, in the format [name1, name2, ... ]."}),
        "nxc": (InputValue, {"dtype": float,
                             "default": 1,
                             "help": "The average number of exchanges per step to be attempted "}),
        "ealc": (InputValue, {"dtype": float,
                              "default": 0.0,
                              "help": "The contribution to the conserved quantity for the atom swapper"})
    }

    dynamic = {}

    default_help = "Holds all the information for doing Monte Carlo atom swap moves. "
    default_label = "ATOMSWAP"

    def store(self, alc):
        """Takes an alchemical exchange instance and stores a minimal representation of it.

        Args:
            alc: An alchemy object.
        """

        if alc == {}:
            return

        self.names.store(alc.names)
        self.nxc.store(alc.nxc)
        self.ealc.store(alc.ealc)

    def fetch(self):
        """Creates an ensemble object.

        Returns:
            An ensemble object of the appropriate mode and with the appropriate
            objects given the attributes of the InputEnsemble object.
        """

        rv = super(InputAtomSwap, self).fetch()
        rv["mode"] = self.mode.fetch()
        return rv

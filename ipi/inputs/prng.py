"""Creates objects that deal with random number generation.

Generates a random number generator either from a seed number, or from a
state vector.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np

from ipi.utils.prng import *
from ipi.utils.inputvalue import *


__all__ = ["InputRandom"]


class InputRandom(Input):

    """Random input class.

    Handles generating the appropriate random number class from the xml
    input file, and generating the xml checkpoint tags and data from an
    instance of the object.

    Attributes:
       seed: An optional integer giving a seed to initialise the random number
          generator from. Defaults to 123456.
       state: An optional array giving the state of the random number generator.
          Defaults to an empty array.
       has_gauss: An optional integer giving whether there is a stored
          Gaussian number or not. Defaults to 0.
       gauss: An optional float giving the stored Gaussian number. Defaults to
          0.0.
       set_pos: An optional integer giving the position in the state array
          that is being read from. Defaults to 0.
    """

    fields = {
        "seed": (
            InputValue,
            {
                "dtype": int,
                "default": 123456,
                "help": "This is the seed number used to generate the initial state of the random number generator.",
            },
        ),
        "state": (
            InputArray,
            {
                "dtype": np.uint,
                "default": input_default(
                    factory=np.zeros, kwargs={"shape": (0,), "dtype": np.uint}
                ),
                "help": "Gives the state vector for the random number generator. Avoid directly modifying this unless you are very familiar with the inner workings of the algorithm used.",
            },
        ),
        "has_gauss": (
            InputValue,
            {
                "dtype": int,
                "default": 0,
                "help": "Determines whether there is a stored gaussian number or not. A value of 0 means there is none stored.",
            },
        ),
        "gauss": (
            InputValue,
            {"dtype": float, "default": 0.00, "help": "The stored Gaussian number."},
        ),
        "set_pos": (
            InputValue,
            {
                "dtype": int,
                "default": 0,
                "help": "Gives the position in the state array that the random number generator is reading from.",
            },
        ),
    }

    default_help = "Deals with the pseudo-random number generator."
    default_label = "PRNG"

    def store(self, prng):
        """Takes a random number instance and stores a minimal
        representation of it.

        Args:
           prng: A random number object from which to initialise from.
        """

        super(InputRandom, self).store(prng)
        self.seed.store(prng.seed)
        gstate = prng.state
        self.state.store(gstate[1])
        self.set_pos.store(gstate[2])
        self.has_gauss.store(gstate[3])
        self.gauss.store(gstate[4])

    def fetch(self):
        """Creates a random number object.

        Returns:
           An random number object of the appropriate type and with the
           appropriate properties given the attributes of the InputRandom
           object.
        """

        super(InputRandom, self).fetch()
        if not self.state._explicit:
            return Random(seed=self.seed.fetch())
        else:
            return Random(
                state=(
                    "MT19937",
                    self.state.fetch(),
                    self.set_pos.fetch(),
                    self.has_gauss.fetch(),
                    self.gauss.fetch(),
                )
            )

"""Creates objects that deal with random number generation.

Generates a random number generator either from a seed number, or from a
state vector.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import json

from ipi.utils.prng import *
from ipi.utils.inputvalue import *
from ipi.utils.io import NumpyEncoder


__all__ = ["InputRandom"]


class InputRandom(Input):
    """Random input class.

    Handles generating the appropriate random number class from the xml
    input file, and generating the xml checkpoint tags and data from an
    instance of the object.

    Attributes:
       seed: An optional integer giving a seed to initialise the random number
          generator from. Defaults to 123456.
       state: An optional string giving the state of the random number generator.s
       n_threads: An optional integer giving whether to use threaded evaluation. Defaults to 1.
    """

    attribs = {
        "n_threads": (
            InputValue,
            {
                "dtype": int,
                "default": 1,
                "help": """
Use parallel PRNG generator. Will make trajectories less reproducible and is only faster
if the arrays are very large.
""",
            },
        ),
    }
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
            InputValue,
            {
                "dtype": str,
                "default": "",
                "help": "Gives the state vector for the random number generator. Avoid directly modifying this unless you are very familiar with the inner workings of the algorithm used.",
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
        self.state.store(json.dumps(prng.state, cls=NumpyEncoder))
        self.n_threads.store(prng.n_threads)

    def fetch(self):
        """Creates a random number object.

        Returns:
           An random number object of the appropriate type and with the
           appropriate properties given the attributes of the InputRandom
           object.
        """

        super(InputRandom, self).fetch()
        if not self.state._explicit:
            return Random(seed=self.seed.fetch(), n_threads=self.n_threads.fetch())
        else:
            return Random(
                state=json.loads(self.state.fetch()), n_threads=self.n_threads.fetch()
            )

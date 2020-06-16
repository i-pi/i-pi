"""Classes used to generate pseudo-random numbers.

Allows the user to specify a seed for the random number generator.
These are used in initialising the velocities and in stochastic thermostats.
The state of the random number generator is kept track of, so that the if the
simulation is restarted from a checkpoint, we will see the same dynamics as if
it had not been stopped.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np


__all__ = ["Random"]


class Random(object):

    """Class to interface with the standard pseudo-random number generator.

    Initialises the standard numpy pseudo-random number generator from a seed
    at the beginning of the simulation, and keeps track of the state so that
    it can be output to the checkpoint files throughout the simulation.

    Attributes:
        rng: The random number generator to be used.
        seed: The seed number to start the generator.
        state: A tuple of five objects giving the current state of the random
            number generator. The first is the type of random number generator,
            here 'MT19937', the second is an array of 624 integers, the third
            is the current position in the array that is being read from, the
            fourth gives whether it has a gaussian random number stored, and
            the fifth is this stored Gaussian random number, or else the last
            Gaussian random number returned.
    """

    def __init__(self, seed=12345, state=None):
        """Initialises Random.

        Args:
            seed: An optional seed giving an integer to initialise the state with.
            state: An optional state tuple to initialise the state with.
        """

        self.rng = np.random.mtrand.RandomState(seed=seed)
        self.seed = seed
        if state is None:
            self.rng.seed(seed)
        else:
            self.state = state

    def get_state(self):
        """Interface to the standard get_state() function."""

        return self.rng.get_state()

    def set_state(self, value):
        """Interface to the standard set_state() function.

        Should only be used with states generated from another similar random
        number generator, such as one from a previous run.
        """

        return self.rng.set_state(value)

    state = property(get_state, set_state)

    @property
    def u(self):
        """Interface to the standard random_sample() function.

        Returns:
            A pseudo-random number from a uniform distribution from 0-1.
        """

        return self.rng.random_sample()

    @property
    def g(self):
        """Interface to the standard standard_normal() function.

        Returns:
            A pseudo-random number from a normal Gaussian distribution.
        """

        return self.rng.standard_normal()

    def gamma(self, k, theta=1.0):
        """Interface to the standard gamma() function.

        Args:
            k: Shape parameter for the gamma distribution.
            theta: Mean of the distribution.

        Returns:
            A random number from a gamma distribution with a shape k and a
            mean value theta.
        """

        return self.rng.gamma(k, theta)

    def gvec(self, shape):
        """Interface to the standard_normal array function.

        Args:
            shape: The shape of the array to be returned.

        Returns:
            An array with the required shape where each element is taken from
            a normal Gaussian distribution.
        """

        return self.rng.standard_normal(shape)

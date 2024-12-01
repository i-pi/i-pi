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
import concurrent.futures
import time

__all__ = ["Random"]

_MIN_STEP_THREADED = 100  # minimum stride to use multithreaded prng


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

    def __init__(self, seed=-1, state=None, n_threads=1):
        """Initialises Random.

        Args:
            seed: An optional seed giving an integer to initialise the state with.
            state: An optional state tuple to initialise the state with.
            n_threads: Whether to use multi-threaded generation (only useful if
                    arrays to be filled are large!)
        """

        # Here we make the seed random if it has not been specified in the input file
        if seed == -1:
            seed = int(time.time() * 1000)

        self.seed = seed

        self.rng = [
            np.random.Generator(np.random.MT19937(s))
            for s in np.random.SeedSequence(seed).spawn(n_threads)
        ]

        self.n_threads = n_threads
        if self.n_threads == 1:
            self.gvec = self.gvec_serial
            self.gfill = self.gfill_serial
        else:
            self.executor = concurrent.futures.ThreadPoolExecutor(n_threads)
            self.gvec = self.gvec_threaded
            self.gfill = self.gfill_threaded

        if state is not None:
            self.state = state

    def get_state(self):
        """Interface to the standard get_state() function."""

        return [r.bit_generator.state for r in self.rng]

    def set_state(self, value):
        """Interface to the standard set_state() function.

        Should only be used with states generated from another similar random
        number generator, such as one from a previous run.
        """

        for r, s in zip(self.rng, value):
            r.bit_generator.state = s

    state = property(get_state, set_state)

    @property
    def u(self):
        """Interface to the standard random_sample() function.

        Returns:
            A pseudo-random number from a uniform distribution from 0-1.
        """

        return self.rng[0].uniform()

    @property
    def g(self):
        """Interface to the standard standard_normal() function.

        Returns:
            A pseudo-random number from a normal Gaussian distribution.
        """

        return self.rng[0].standard_normal()

    def gamma(self, k, theta=1.0, size=None):
        """Interface to the standard gamma() function.

        Args:
            k: Shape parameter for the gamma distribution.
            theta: Mean of the distribution.

        Returns:
            A random number from a gamma distribution with a shape k and a
            mean value theta.
        """

        return self.rng[0].gamma(k, theta, size)

    def poisson(self, lam=1.0, size=None):
        """Interface to the standard poisson() function.

        Args:
            lam: Mean of the Poisson distribution

        Returns:
            A random number from a Poisson distribution
        """

        return self.rng[0].poisson(lam, size)

    def uniform(self, low=0.0, high=1.0, size=None):
        """Interface to the standard uniform() function.

        Args:
            Same as numpy.Generator.uniform

        Returns:
            Uniform random reals in the prescribed interval
        """

        return self.rng[0].uniform(low, high, size)

    def integers(self, low, high=None, size=None, dtype=np.int64, endpoint=False):
        """Interface to the standard integers() function.

        Args:
            Same as numpy.Generator.integers

        Returns:
            Random integers in the prescribed interval
        """

        return self.rng[0].integers(low, high, size, dtype, endpoint)

    def shuffle(self, x, axis=0):
        """Interface to the standard shuffle() function.

        Args:
            Same as numpy.Generator.shuffle

        Returns:
            None
        """

        self.rng[0].shuffle(x, axis)

    def gfill_serial(self, out):
        """Fills a pre-allocated array serially

        Args:
            out: The array to be filled.
        """

        self.rng[0].standard_normal(out=out)

    def gfill_threaded(self, out):
        """Fills a pre-allocated array in parallel

        Args:
            out: The array to be filled.
        """

        out_flat = out.flatten()
        step_size = np.ceil(len(out_flat) / self.n_threads).astype(int)
        if step_size < _MIN_STEP_THREADED:
            # falls back to serial if the vector is too small
            self.gfill_serial(out)
        else:
            # threaded execution
            futures = {}
            for i in range(self.n_threads):
                futures[
                    self.executor.submit(
                        self.rng[i].standard_normal,
                        out=out_flat[step_size * i : step_size * (i + 1)],
                    )
                ] = i
            concurrent.futures.wait(futures)

    def gvec_serial(self, shape):
        """Interface to the standard_normal array function.

        Args:
            shape: The shape of the array to be returned.

        Returns:
            An array with the required shape where each element is taken from
            a normal Gaussian distribution.
        """

        rvec = np.empty(shape=shape)
        self.gfill_serial(rvec)
        return rvec

    def gvec_threaded(self, shape):
        """Interface to the standard_normal array function.

        Args:
            shape: The shape of the array to be returned.

        Returns:
            An array with the required shape where each element is taken from
            a normal Gaussian distribution.
        """

        rvec = np.empty(shape=shape)
        self.gfill_threaded(rvec)

        return rvec

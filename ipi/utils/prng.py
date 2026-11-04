"""Classes used to generate pseudo-random numbers.

Allows the user to specify a seed for the random number generator.
These are used in initialising the velocities and in stochastic thermostats.
The state of the random number generator is kept track of, so that the if the
simulation is restarted from a checkpoint, we will see the same dynamics as if
it had not been stopped.

Scalars and host-side arrays are drawn with numpy generators. Arrays of
Gaussian numbers (`gvec`, `gfill`) are drawn with the generator of the array
backend, i.e. on the device when using torch. With JAX, whose arrays cannot be
filled in place, only `gvec` is available.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np
import concurrent.futures
import functools
import time

from ipi.utils.array_backend import backend, device
from ipi.utils.messages import verbosity, warning

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
        self._torch_rng = None
        self._jax_key = None
        if backend() == "jax":
            import jax

            # Gaussian arrays are drawn with a JAX key, split at each draw and
            # seeded independently of the numpy generators
            self._jax_key = jax.random.key(
                int(np.random.SeedSequence(seed).generate_state(1)[0])
            )
            self.gvec = self.gvec_jax
            self.gammavec = self.gammavec_jax
        elif backend() == "torch":
            import torch

            # Gaussian arrays are drawn on the device by a torch generator,
            # seeded independently of the numpy generators
            self._torch_rng = torch.Generator(device=device())
            self._torch_rng.manual_seed(
                int(np.random.SeedSequence(seed).generate_state(1, np.uint64)[0])
            )
            self.gvec = self.gvec_torch
            self.gfill = self.gfill_torch
            self.gammavec = self.gammavec_torch
        elif self.n_threads == 1:
            self.gvec = self.gvec_serial
            self.gfill = self.gfill_serial
        else:
            self.executor = concurrent.futures.ThreadPoolExecutor(n_threads)
            self.gvec = self.gvec_threaded
            self.gfill = self.gfill_threaded

        if state is not None:
            self.state = state

    @property
    def torch_generator(self):
        """The torch generator that draws the Gaussian arrays on the device, or
        None if they are drawn with numpy."""

        return self._torch_rng

    def get_state(self):
        """Interface to the standard get_state() function.

        Returns the list of the states of the numpy generators, or a dictionary
        that also contains the state of the torch generator or of the JAX key,
        if present.
        """

        states = [r.bit_generator.state for r in self.rng]
        if self._jax_key is not None:
            import jax

            return {"numpy": states, "jax": jax.random.key_data(self._jax_key).tolist()}
        if self._torch_rng is None:
            return states
        return {
            "numpy": states,
            "torch": {
                "device": str(self._torch_rng.device),
                "state": self._torch_rng.get_state().tolist(),
            },
        }

    def set_state(self, value):
        """Interface to the standard set_state() function.

        Should only be used with states generated from another similar random
        number generator, such as one from a previous run. If the state of the
        torch generator is missing, or refers to a different device, the torch
        generator keeps the state obtained from the seed.
        """

        torch_state = jax_state = None
        if isinstance(value, dict):
            torch_state = value.get("torch")
            jax_state = value.get("jax")
            value = value["numpy"]
        for r, s in zip(self.rng, value):
            r.bit_generator.state = s
        if self._jax_key is not None:
            import jax

            if jax_state is not None:
                self._jax_key = jax.random.wrap_key_data(
                    jax.numpy.asarray(jax_state, dtype=jax.numpy.uint32)
                )
            else:
                warning(
                    "No state for the JAX random key: it is initialized from the seed",
                    verbosity.low,
                )
        if self._torch_rng is not None:
            if torch_state is not None and torch_state["device"] == str(
                self._torch_rng.device
            ):
                import torch

                self._torch_rng.set_state(
                    torch.tensor(torch_state["state"], dtype=torch.uint8)
                )
            else:
                warning(
                    "No state for the torch random number generator on this device: "
                    "it is initialized from the seed",
                    verbosity.low,
                )

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

    def gammavec(self, k, shape):
        """Returns an array of random numbers from a gamma distribution with
        shape `k` and scale 1, drawn with the generator of the backend.

        Args:
            k: Shape parameter for the gamma distribution.
            shape: The shape of the array to be returned.
        """

        return self.rng[0].gamma(k, 1.0, shape)

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

        out_flat = out.reshape(-1)  # a view of out whenever out is contiguous
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
            if not np.may_share_memory(out_flat, out):
                out[...] = out_flat.reshape(out.shape)

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

    def gfill_torch(self, out):
        """Fills a pre-allocated torch array with the torch generator.

        Args:
            out: The array to be filled.
        """

        out.normal_(generator=self._torch_rng)

    def gvec_torch(self, shape):
        """Returns an array of Gaussian random numbers drawn with the torch
        generator, on its device.

        Args:
            shape: The shape of the array to be returned.
        """

        import torch

        return torch.randn(
            shape,
            generator=self._torch_rng,
            dtype=torch.float64,
            device=self._torch_rng.device,
        )

    def gvec_jax(self, shape):
        """Returns an array of Gaussian random numbers drawn with a new key
        split from the JAX key.

        Args:
            shape: The shape of the array to be returned.
        """

        self._jax_key, rvec = _jax_normal()(self._jax_key, shape)
        return rvec

    def gammavec_torch(self, k, shape):
        """Returns an array of gamma-distributed random numbers drawn with the
        torch generator, on its device (see `gammavec`)."""

        import torch

        alpha = torch.empty(shape, dtype=torch.float64, device=self._torch_rng.device)
        return torch._standard_gamma(alpha.fill_(k), generator=self._torch_rng)

    def gammavec_jax(self, k, shape):
        """Returns an array of gamma-distributed random numbers drawn with a new
        key split from the JAX key (see `gammavec`)."""

        self._jax_key, rvec = _jax_gamma()(self._jax_key, k, shape)
        return rvec


@functools.cache
def _jax_normal():
    """Returns a jitted function that splits `key` and draws an array of shape
    `shape` with one of the two keys, returning the other key and the array
    (a single call, as eager JAX operations are slow)."""

    import jax

    def normal(key, shape):
        key, subkey = jax.random.split(key)
        return key, jax.random.normal(subkey, shape, dtype=jax.numpy.float64)

    return jax.jit(normal, static_argnums=1)


@functools.cache
def _jax_gamma():
    """Returns a jitted function that splits `key` and draws an array of shape
    `shape` from a gamma distribution with shape `k`, with one of the two keys,
    returning the other key and the array."""

    import jax

    def gamma(key, k, shape):
        key, subkey = jax.random.split(key)
        return key, jax.random.gamma(subkey, k, shape, dtype=jax.numpy.float64)

    return jax.jit(gamma, static_argnums=(1, 2))

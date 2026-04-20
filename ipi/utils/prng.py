"""Classes used to generate pseudo-random numbers.

Allows the user to specify a seed for the random number generator.
These are used in initialising the velocities and in stochastic thermostats.
The state of the random number generator is kept track of, so that the if the
simulation is restarted from a checkpoint, we will see the same dynamics as if
it had not been stopped.

A single stream drives every draw, so the full trajectory is reproducible
from the seed. Under the torch array backend the stream is a
``torch.Generator``; under numpy it is ``numpy.random.Generator`` (MT19937).
Checkpoints tag the state with the active backend — restart requires the
same backend as the run that produced it.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import math
import time
import concurrent.futures

import numpy as np
from array_api_compat import is_torch_namespace

from ipi.utils.array_backend import xp, device as _xp_device
from ipi.utils.messages import warning

__all__ = ["Random"]

_MIN_STEP_THREADED = 100  # minimum stride to use multithreaded prng
_IS_TORCH_BACKEND = is_torch_namespace(xp)


class Random(object):
    """Class to interface with a single-stream pseudo-random generator.

    Attributes:
        seed: The seed number to start the generator.
        rng: Under the numpy backend, the list of ``numpy.random.Generator``
            streams (one per thread). ``None`` under torch.
        n_threads: Number of threads used for bulk gaussian fills. Forced to
            1 under the torch backend (single stream is mandatory for
            reproducibility).
    """

    def __init__(self, seed=-1, state=None, n_threads=1):
        """Initialises Random.

        Args:
            seed: An optional seed giving an integer to initialise the state with.
            state: An optional state to initialise the state with.
            n_threads: Number of threads for parallel gaussian fills (numpy
                backend only). Ignored (with a warning) under torch.
        """

        if seed == -1:
            seed = int(time.time() * 1000)

        self.seed = seed

        if _IS_TORCH_BACKEND:
            import torch

            if n_threads != 1:
                warning(
                    "torch backend uses a single PRNG stream; "
                    "threaded PRNG is disabled (n_threads forced to 1)."
                )
                n_threads = 1
            # Generator must live on the device where sampled tensors
            # are allocated (torch's default device). Mismatched devices
            # raise at sample time.
            self._torch_gen = torch.Generator(device=_xp_device)
            self._torch_gen.manual_seed(int(seed))
            self.rng = None
        else:
            self.rng = [
                np.random.Generator(np.random.MT19937(s))
                for s in np.random.SeedSequence(seed).spawn(n_threads)
            ]
            self._torch_gen = None

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
        """Current state of the random stream.

        Tagged with the active backend so cross-backend restarts can be
        detected and rejected.
        """

        if _IS_TORCH_BACKEND:
            data = self._torch_gen.get_state().cpu().numpy().tolist()
            return {"backend": "torch", "data": data}
        return [r.bit_generator.state for r in self.rng]

    def set_state(self, value):
        """Restore a previously saved state.

        Raises ``ValueError`` if the saved state was produced by a different
        array backend — those states are not interchangeable.
        """

        saved_is_torch = isinstance(value, dict) and value.get("backend") == "torch"
        if _IS_TORCH_BACKEND:
            if not saved_is_torch:
                raise ValueError(
                    "PRNG state was saved under the numpy backend; restart "
                    "requires the same array backend as the original run."
                )
            import torch

            data = torch.tensor(value["data"], dtype=torch.uint8)
            self._torch_gen.set_state(data)
            return
        if saved_is_torch:
            raise ValueError(
                "PRNG state was saved under the torch backend; restart "
                "requires the same array backend as the original run."
            )
        for r, s in zip(self.rng, value):
            r.bit_generator.state = s

    state = property(get_state, set_state)

    @property
    def u(self):
        """A scalar draw from U(0, 1)."""

        if _IS_TORCH_BACKEND:
            import torch

            return torch.rand((), generator=self._torch_gen).item()
        return self.rng[0].uniform()

    @property
    def g(self):
        """A scalar draw from the standard normal."""

        if _IS_TORCH_BACKEND:
            import torch

            return torch.randn((), generator=self._torch_gen).item()
        return self.rng[0].standard_normal()

    def gamma(self, k, theta=1.0, size=None):
        """Draw from a gamma distribution (shape ``k``, scale ``theta``).

        Torch's ``distributions.Gamma.sample`` does not accept a generator,
        so under torch we implement Marsaglia-Tsang manually against the
        single stream. Only scalar sampling (``size=None``) is currently
        needed by i-PI under torch.
        """

        if _IS_TORCH_BACKEND:
            if size is not None:
                raise NotImplementedError(
                    "prng.gamma(size=...) is not implemented under the torch "
                    "backend; only scalar sampling is supported."
                )
            return float(theta) * self._gamma_scalar_torch(float(k))
        return self.rng[0].gamma(k, theta, size)

    def _gamma_scalar_torch(self, k):
        """Scalar gamma(k, 1) via Marsaglia-Tsang, consuming self._torch_gen."""

        import torch

        if k < 1.0:
            # boost: X ~ Gamma(k) <=> X = Y * U**(1/k) with Y ~ Gamma(k+1)
            y = self._gamma_scalar_torch(k + 1.0)
            u = torch.rand((), generator=self._torch_gen).item()
            return y * (u ** (1.0 / k))
        d = k - 1.0 / 3.0
        c = 1.0 / math.sqrt(9.0 * d)
        while True:
            x = torch.randn((), generator=self._torch_gen).item()
            v = (1.0 + c * x) ** 3
            if v <= 0.0:
                continue
            u = torch.rand((), generator=self._torch_gen).item()
            x2 = x * x
            if u < 1.0 - 0.0331 * x2 * x2:
                return d * v
            if math.log(u) < 0.5 * x2 + d * (1.0 - v + math.log(v)):
                return d * v

    def poisson(self, lam=1.0, size=None):
        """Draw from a Poisson distribution with mean ``lam``."""

        if _IS_TORCH_BACKEND:
            import torch

            shape = _normalise_size(size)
            rates = torch.full(shape, float(lam)) if shape else torch.tensor(float(lam))
            out = torch.poisson(rates, generator=self._torch_gen)
            if size is None:
                return int(out.item())
            return out.to(torch.int64).cpu().numpy()
        return self.rng[0].poisson(lam, size)

    def uniform(self, low=0.0, high=1.0, size=None):
        """Uniform reals on [low, high)."""

        if _IS_TORCH_BACKEND:
            import torch

            if size is None:
                return (
                    low
                    + (high - low) * torch.rand((), generator=self._torch_gen).item()
                )
            shape = _normalise_size(size)
            out = torch.empty(shape)
            out.uniform_(low, high, generator=self._torch_gen)
            return out.cpu().numpy()
        return self.rng[0].uniform(low, high, size)

    def integers(self, low, high=None, size=None, dtype=np.int64, endpoint=False):
        """Random integers in the prescribed interval."""

        if _IS_TORCH_BACKEND:
            import torch

            # numpy convention: integers(high) draws from [0, high); integers(low, high)
            # draws from [low, high). Mirror that.
            if high is None:
                lo, hi = 0, low
            else:
                lo, hi = low, high
            if endpoint:
                hi = hi + 1
            shape = _normalise_size(size)
            out = torch.randint(
                int(lo), int(hi), shape, generator=self._torch_gen, dtype=torch.int64
            )
            if size is None:
                return int(out.item())
            return out.cpu().numpy().astype(dtype, copy=False)
        return self.rng[0].integers(low, high, size, dtype, endpoint)

    def shuffle(self, x, axis=0):
        """In-place shuffle along ``axis``."""

        if _IS_TORCH_BACKEND:
            import torch

            n = x.shape[axis]
            perm = torch.randperm(n, generator=self._torch_gen).cpu().numpy()
            x[...] = np.take(x, perm, axis=axis)
            return
        self.rng[0].shuffle(x, axis)

    def gfill_serial(self, out):
        """Fills a pre-allocated array serially.

        Args:
            out: The array (numpy) to be filled.
        """

        if _IS_TORCH_BACKEND:
            import torch

            buf = torch.empty(tuple(out.shape))
            buf.normal_(generator=self._torch_gen)
            out[...] = buf.cpu().numpy()
            return
        self.rng[0].standard_normal(out=out)

    def gfill_threaded(self, out):
        """Fills a pre-allocated array in parallel (numpy backend only)."""

        out_flat = out.flatten()
        step_size = np.ceil(len(out_flat) / self.n_threads).astype(int)
        if step_size < _MIN_STEP_THREADED:
            # falls back to serial if the vector is too small
            self.gfill_serial(out)
        else:
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
        """Gaussian array in the active backend (serial)."""

        if _IS_TORCH_BACKEND:
            import torch

            out = torch.empty(tuple(shape) if not isinstance(shape, int) else (shape,))
            out.normal_(generator=self._torch_gen)
            return out
        rvec = np.empty(shape=shape)
        self.gfill_serial(rvec)
        return xp.asarray(rvec)

    def gvec_threaded(self, shape):
        """Gaussian array in the active backend (numpy-threaded)."""

        rvec = np.empty(shape=shape)
        self.gfill_threaded(rvec)
        return xp.asarray(rvec)


def _normalise_size(size):
    """numpy-style ``size`` argument to a tuple shape (``()`` for scalar)."""

    if size is None:
        return ()
    if isinstance(size, int):
        return (size,)
    return tuple(size)

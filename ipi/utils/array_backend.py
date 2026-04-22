"""Array-API backend for i-PI.

Exposes a single, process-global namespace `xp` that the numeric modules
use in place of `numpy`. Consumers do:

    from ipi.utils.array_backend import xp

    y = xp.asarray(...)

The selected backend is fixed for the lifetime of the process. It is
resolved once, in this order:

1. the value of the environment variable ``IPI_ARRAY_BACKEND``
   (``numpy`` / ``torch`` / ``jax``), if set
2. otherwise, numpy

``set_array_backend(...)`` can override the choice programmatically, but
only if it runs *before* any numeric module has done
``from ipi.utils.array_backend import xp`` — otherwise those modules
have already captured the previous binding. The ``i-pi`` CLI honours a
``--array-backend`` flag by setting ``IPI_ARRAY_BACKEND`` in the
environment before any ipi import, so this restriction never bites
end users of the command-line tool.

Mid-run backend swaps are not supported: cached matrices, RNG state and
similar are bound to whatever backend was active at construction.

I/O and socket code paths keep using numpy directly and convert with
``to_numpy()`` at the boundary.
"""

import os

import numpy as np

from array_api_compat import size as xp_size

__all__ = [
    "xp",
    "xp_size",
    "device",
    "dtype",
    "set_array_backend",
    "array_namespace",
    "to_numpy",
]


_CREATION_FNS = frozenset(
    {
        "asarray",
        "zeros",
        "zeros_like",
        "ones",
        "ones_like",
        "empty",
        "empty_like",
        "full",
        "full_like",
        "arange",
        "linspace",
        "eye",
    }
)


class _XpWithDevice:
    """Thin proxy that injects ``device=`` into array-creation functions.

    Avoids calling ``torch.set_default_device(...)``, which installs a
    global ``__torch_function__`` mode that fires on every aten op and
    adds ~3 us of Python dispatch overhead per kernel. Non-creation
    attributes are forwarded as-is; creation functions get a thin
    closure that sets ``device=`` if the caller didn't pass one.
    """

    __slots__ = ("_xp", "_device", "_cache")

    def __init__(self, xp, device):
        self._xp = xp
        self._device = device
        # Cache attribute lookups so hot-path code that reads `xp.foo`
        # repeatedly doesn't keep repaying the wrapping cost.
        self._cache = {}

    def __getattr__(self, name):
        cache = self._cache
        hit = cache.get(name)
        if hit is not None:
            return hit
        attr = getattr(self._xp, name)
        if name in _CREATION_FNS and callable(attr):
            device = self._device

            def wrapped(*args, **kwargs):
                kwargs.setdefault("device", device)
                return attr(*args, **kwargs)

            cache[name] = wrapped
            return wrapped
        cache[name] = attr
        return attr


def _load_backend(name):
    """Load backend namespace and return (namespace, device, dtype).

    `device` and `dtype` are the resolved defaults (strings / backend
    dtype objects, or None where the backend has no notion of them).
    """
    if name == "numpy":
        import array_api_compat.numpy as _numpy_backend

        return _numpy_backend, "cpu", None
    if name == "torch":
        import torch
        import array_api_compat.torch as _torch_backend

        # IPI_DEVICE / IPI_DTYPE apply only to torch. We avoid
        # `torch.set_default_device`: it installs a global
        # ``__torch_function__`` mode that fires on *every* aten op
        # and adds measurable per-kernel dispatch overhead. Instead,
        # wrap the backend namespace so that the handful of creation
        # functions (``asarray``, ``zeros``, ``empty``, ...) receive
        # ``device=`` transparently.
        dev = os.environ.get("IPI_DEVICE", "cpu")
        # i-PI is double-precision throughout (matches numpy's default for
        # float operations); torch otherwise defaults to float32 which
        # would cause dtype mismatches when a float64-typed array meets
        # a prng-drawn / default-dtype array. `IPI_DTYPE=float32` opts in
        # explicitly.
        dt_name = os.environ.get("IPI_DTYPE", "float64")
        torch.set_default_dtype(getattr(torch, dt_name))
        # Prime the CUDA primary context eagerly. cuFFT (and other
        # lazy-init kernels) otherwise create the context from
        # whichever worker thread hits them first and emit a
        # UserWarning about it. Allocating any tensor on the target
        # device forces primary-context creation now; subsequent
        # threads inherit it implicitly via `cudaSetDevice`.
        if dev != "cpu" and dev.startswith("cuda") and torch.cuda.is_available():
            torch.cuda.init()
            torch.empty(0, device=dev)
        return _XpWithDevice(_torch_backend, dev), dev, torch.get_default_dtype()
    if name == "jax":
        import array_api_compat.jax.numpy as _jax_backend

        return _jax_backend, None, None
    raise ValueError(f"Unknown array backend: {name!r}")


xp, device, dtype = _load_backend(os.environ.get("IPI_ARRAY_BACKEND", "numpy"))


def set_array_backend(backend):
    """Rebind the process-global array backend.

    ``backend`` is either a backend name (``"numpy"``/``"torch"``/``"jax"``)
    or an already-imported array-API-compatible module. torch/jax are
    only imported if explicitly requested, so users who don't opt in
    never pay the dependency.

    Must be called before any numeric module has imported ``xp``; see
    the module docstring.
    """
    global xp, device, dtype
    if isinstance(backend, str):
        xp, device, dtype = _load_backend(backend)
    else:
        xp = backend


def array_namespace(*xs):
    """Namespace of one or more arrays (thin wrapper over array_api_compat)."""
    from array_api_compat import array_namespace as _ns

    return _ns(*xs)


def to_numpy(x):
    """Convert an array to a plain ``numpy.ndarray``.

    No-op if ``x`` is already a numpy array or a plain scalar. Used at
    I/O and socket boundaries where numpy is expected.
    """
    # Unwrap depend_array so array_api_compat sees the raw backend
    # array and can dispatch on its real type.
    if hasattr(x, "_value"):
        x = x._value
    if isinstance(x, np.ndarray):
        return x
    try:
        return np.asarray(x)
    except TypeError:
        from array_api_compat import to_device

        return np.asarray(to_device(x, "cpu"))

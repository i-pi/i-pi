"""Array-API backend resolution for i-PI.

Exposes a single, process-global namespace `xp` that the numeric modules
use instead of `numpy` directly. Defaults to numpy (via
`array_api_compat.numpy`) so that existing behavior is unchanged.

Switch backends with `set_default_backend("torch")` (or pass the module
itself). Optional backends (torch, jax) are imported lazily on demand,
so they are not required unless the user asks for them.

I/O and socket code paths keep using numpy directly and convert with
`to_numpy()` at the boundary.
"""

import numpy as np


__all__ = [
    "get_xp",
    "set_default_backend",
    "array_namespace",
    "to_numpy",
    "asarray",
]


_default_xp = None  # resolved lazily on first get_xp() call


def _load_numpy_backend():
    import array_api_compat.numpy as xp_numpy

    return xp_numpy


def get_xp():
    """Return the current default array-API namespace."""
    global _default_xp
    if _default_xp is None:
        _default_xp = _load_numpy_backend()
    return _default_xp


def set_default_backend(backend):
    """Set the process-global array backend.

    `backend` can be the string "numpy"/"torch"/"jax" or an already-
    imported array-api-compatible module. torch/jax are only imported
    if explicitly requested here.
    """
    global _default_xp
    if isinstance(backend, str):
        if backend == "numpy":
            _default_xp = _load_numpy_backend()
        elif backend == "torch":
            import array_api_compat.torch as xp_torch

            _default_xp = xp_torch
        elif backend == "jax":
            import array_api_compat.jax.numpy as xp_jax

            _default_xp = xp_jax
        else:
            raise ValueError(f"Unknown backend: {backend!r}")
    else:
        _default_xp = backend


def array_namespace(*xs):
    """Namespace of one or more arrays (thin wrapper over array_api_compat)."""
    from array_api_compat import array_namespace as _ns

    return _ns(*xs)


def to_numpy(x):
    """Convert an array to a plain `numpy.ndarray`.

    No-op if `x` is already a numpy array or a plain scalar. Used at
    I/O and socket boundaries where numpy is expected.
    """
    if isinstance(x, np.ndarray):
        return x
    try:
        return np.asarray(x)
    except TypeError:
        from array_api_compat import to_device

        return np.asarray(to_device(x, "cpu"))


def asarray(x, dtype=None):
    """Convert an object to an array in the default backend."""
    xp = get_xp()
    if dtype is None:
        return xp.asarray(x)
    return xp.asarray(x, dtype=dtype)

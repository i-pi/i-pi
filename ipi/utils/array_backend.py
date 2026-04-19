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

__all__ = [
    "xp",
    "set_array_backend",
    "array_namespace",
    "to_numpy",
]


def _load_backend(name):
    if name == "numpy":
        import array_api_compat.numpy as _numpy_backend

        return _numpy_backend
    if name == "torch":
        import torch
        import array_api_compat.torch as _torch_backend

        # IPI_DEVICE / IPI_DTYPE apply only to torch. Setting the defaults
        # here makes subsequent xp.zeros/xp.asarray land on the right
        # device and dtype without plumbing them through every call site.
        device = os.environ.get("IPI_DEVICE")
        if device:
            torch.set_default_device(device)
        dtype = os.environ.get("IPI_DTYPE")
        if dtype:
            torch.set_default_dtype(getattr(torch, dtype))
        return _torch_backend
    if name == "jax":
        import array_api_compat.jax.numpy as _jax_backend

        return _jax_backend
    raise ValueError(f"Unknown array backend: {name!r}")


xp = _load_backend(os.environ.get("IPI_ARRAY_BACKEND", "numpy"))


def set_array_backend(backend):
    """Rebind the process-global array backend.

    ``backend`` is either a backend name (``"numpy"``/``"torch"``/``"jax"``)
    or an already-imported array-API-compatible module. torch/jax are
    only imported if explicitly requested, so users who don't opt in
    never pay the dependency.

    Must be called before any numeric module has imported ``xp``; see
    the module docstring.
    """
    global xp
    if isinstance(backend, str):
        xp = _load_backend(backend)
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
    if isinstance(x, np.ndarray):
        return x
    try:
        return np.asarray(x)
    except TypeError:
        from array_api_compat import to_device

        return np.asarray(to_device(x, "cpu"))

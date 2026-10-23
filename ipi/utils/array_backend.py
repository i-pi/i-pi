"""Selection of the array library used to store and process i-PI data.

The functions of the selected library are collected in the namespace `xp`,
used in place of numpy for operations on the arrays that hold the state of
a simulation:

    from ipi.utils.array_backend import xp

    v = xp.zeros(3)

`set_backend` selects the library and the device. It updates `xp` in place,
so it also affects the modules that have already imported it, but it does not
convert existing arrays: it must be called before a simulation is created.

Supported backends are "numpy" (the default) and "torch". For numpy, `xp`
contains the numpy namespace. For torch, it contains the array-API compatible
namespace from `array_api_compat`, with the array-creation functions
allocating on the selected device, and choosing the type of the arrays like
numpy does, i.e. with float64 as the floating-point type. The global default
type of torch is left unchanged.

Input and output, sockets and forcefields use numpy arrays: `to_numpy` converts
backend arrays at these boundaries.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

import types

import numpy as np
import array_api_compat

__all__ = [
    "xp",
    "xp_size",
    "set_backend",
    "backend",
    "device",
    "to_numpy",
    "synchronize",
]

xp = types.SimpleNamespace()
xp_size = array_api_compat.size

_current = {"backend": None, "device": None}


def _creation_functions(namespace, device, dtype):
    """Returns the array-creation functions of `namespace`, modified to allocate
    on `device` and to choose the type of the arrays like numpy, using `dtype`
    as the floating-point type. Functions `*_like` follow their argument."""

    def always(args, kwargs):
        return True

    def any_float(args, kwargs):
        return any(isinstance(v, float) for v in (*args, *kwargs.values()))

    def with_defaults(func, floating):
        def wrapped(*args, **kwargs):
            kwargs.setdefault("device", device)
            if kwargs.get("dtype") is None and floating(args, kwargs):
                kwargs["dtype"] = dtype
            return func(*args, **kwargs)

        return wrapped

    def asarray(obj, *args, **kwargs):
        # Python scalars and sequences are converted by numpy, to infer their type
        if kwargs.get("dtype") is None and not array_api_compat.is_array_api_obj(obj):
            obj = np.asarray(obj)
        kwargs.setdefault("device", device)
        return namespace.asarray(obj, *args, **kwargs)

    functions = {"asarray": asarray}
    for name in ["zeros", "ones", "empty", "eye", "linspace"]:
        functions[name] = with_defaults(getattr(namespace, name), always)
    for name in ["full", "arange"]:
        functions[name] = with_defaults(getattr(namespace, name), any_float)
    return functions


def set_backend(name="numpy", device="cpu"):
    """Selects the array library and the device used to store arrays.

    Args:
        name: "numpy" or "torch".
        device: The device on which arrays are allocated, e.g. "cpu" or
            "cuda". The numpy backend only supports "cpu".
    """

    if name == "numpy":
        if device != "cpu":
            raise ValueError("The numpy backend only supports the 'cpu' device")
        members = dict(vars(np))
    elif name == "torch":
        import torch
        import array_api_compat.torch as namespace

        members = dict(vars(namespace))
        members.update(_creation_functions(namespace, device, torch.float64))
        if device.startswith("cuda"):
            # creates the CUDA context and the lazily-initialized linear
            # algebra kernels from the main thread, as initializing them
            # concurrently from several threads fails
            torch.cuda.init()
            torch.linalg.eigh(torch.eye(2, device=device))
    else:
        raise ValueError(f"Unknown array backend '{name}'")

    xp.__dict__.clear()
    xp.__dict__.update(members)
    _current["backend"] = name
    _current["device"] = device


def backend():
    """Name of the selected array library."""

    return _current["backend"]


def device():
    """Device on which arrays are allocated."""

    return _current["device"]


def to_numpy(x):
    """Converts an array (or a scalar) to numpy, copying it to the host if needed."""

    if isinstance(x, np.ndarray):
        return x
    if array_api_compat.is_torch_array(x):
        # explicit conversion, which does not go through Tensor.__array__
        return x.cpu().numpy()
    if array_api_compat.is_array_api_obj(x):
        x = array_api_compat.to_device(x, "cpu")
    return np.asarray(x)


def synchronize():
    """Waits for the completion of the operations queued on the device."""

    if _current["backend"] == "torch" and _current["device"] != "cpu":
        import torch

        torch.cuda.synchronize()


set_backend("numpy")

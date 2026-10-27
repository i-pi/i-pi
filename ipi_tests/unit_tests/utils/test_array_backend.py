"""Tests for the selection of the array backend."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

import numpy as np
import pytest

import ipi.utils.array_backend as ab
import ipi.utils.depend as dp
from ipi.utils.array_backend import xp

try:
    import torch

    cuda_available = torch.cuda.is_available()
except ImportError:
    cuda_available = False


@pytest.fixture
def torch_backend():
    """Selects the torch backend on CPU, and restores numpy afterwards."""

    torch = pytest.importorskip("torch")
    ab.set_backend("torch")
    yield torch
    ab.set_backend("numpy")


def test_numpy_default():
    """The default backend is numpy, and conversions are no-ops"""

    assert ab.backend() == "numpy" and ab.device() == "cpu"
    a = xp.zeros(3)
    assert type(a) is np.ndarray
    assert ab.to_numpy(a) is a
    with pytest.raises(ValueError):
        ab.set_backend("numpy", device="cuda")
    with pytest.raises(ValueError):
        ab.set_backend("unknown")


def test_switch_updates_imported_namespace(torch_backend):
    """Modules that imported `xp` see the selected backend"""

    assert xp is ab.xp
    a = xp.zeros(3)
    assert isinstance(a, torch_backend.Tensor)
    assert a.dtype == torch_backend.float64
    assert type(ab.to_numpy(a)) is np.ndarray


def test_types_like_numpy(torch_backend):
    """Array types are chosen like numpy, without changing torch's default"""

    torch = torch_backend
    default = torch.get_default_dtype()
    for created, reference in [
        (xp.zeros(2), np.zeros(2)),
        (xp.eye(2), np.eye(2)),
        (xp.linspace(0, 1, 3), np.linspace(0, 1, 3)),
        (xp.full(2, 1), np.full(2, 1)),
        (xp.full(2, 1.0), np.full(2, 1.0)),
        (xp.arange(3), np.arange(3)),
        (xp.arange(3.0), np.arange(3.0)),
        (xp.asarray(0.5), np.asarray(0.5)),
        (xp.asarray([1, 2]), np.asarray([1, 2])),
        (xp.asarray(np.zeros(2, np.float32)), np.zeros(2, np.float32)),
    ]:
        assert ab.to_numpy(created).dtype == reference.dtype
    assert torch.get_default_dtype() == default


def test_depend_storage(torch_backend):
    """depend_array stores arrays with the selected backend, unless `host` is set"""

    class Holder:
        def __init__(self):
            self.x_ = dp.depend_array(name="x", value=np.zeros(4))
            self.y_ = dp.depend_array(
                name="y",
                value=np.zeros(4),
                func=lambda: 2 * self.x,
                dependencies=[self.x_],
            )
            self.names_ = dp.depend_array(
                name="names", value=np.zeros(4, np.dtype("|U6")), host=True
            )

    dp.dproperties(Holder, ["x", "y", "names"])

    h = Holder()
    assert isinstance(h.x, torch_backend.Tensor)
    assert type(h.names) is np.ndarray
    h.x_[1:3] += 1.0
    h.x_[0] = 3.0
    assert np.array_equal(ab.to_numpy(h.y), [6.0, 2.0, 2.0, 0.0])
    h.x = xp.ones(4)
    assert np.array_equal(ab.to_numpy(h.y), 2 * np.ones(4))
    # numpy arrays are not converted when writing to an existing array
    with pytest.raises(TypeError):
        h.x = np.zeros(4)


@pytest.mark.skipif(not cuda_available, reason="CUDA not available")
def test_cuda_device(torch_backend):
    """Arrays are allocated on the selected device"""

    ab.set_backend("torch", device="cuda")
    a = dp.depend_array(name="a", value=np.arange(3.0))
    assert a.get().device.type == "cuda"
    assert xp.zeros(2).device.type == "cuda"
    assert np.array_equal(ab.to_numpy(a.get()), np.arange(3.0))

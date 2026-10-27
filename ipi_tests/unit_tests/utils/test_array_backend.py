"""Tests for the selection of the array backend."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

import pickle
from copy import deepcopy

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


def check_types_like_numpy():
    """Checks that the arrays created by `xp` have the types numpy would give"""

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


def test_types_like_numpy(torch_backend):
    """Array types are chosen like numpy, without changing torch's default"""

    default = torch_backend.get_default_dtype()
    check_types_like_numpy()
    assert torch_backend.get_default_dtype() == default


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


@pytest.fixture
def jax_backend():
    """Selects the JAX backend on CPU, and restores numpy afterwards."""

    jax = pytest.importorskip("jax")
    ab.set_backend("jax")
    yield jax
    ab.set_backend("numpy")


def test_jax_types_like_numpy(jax_backend):
    """With JAX, arrays are created with the types numpy would give"""

    assert isinstance(xp.zeros(3), jax_backend.Array)
    check_types_like_numpy()


def test_jax_depend_writes(jax_backend):
    """depend_arrays holding JAX arrays are written by replacing the array, and
    the arrays read from them cannot be modified"""

    a = dp.depend_array(name="a", value=np.zeros(4))
    b = dp.depend_array(
        name="b", value=np.zeros(4), func=lambda: 2 * a.get(), dependencies=[a]
    )
    assert type(a) is dp.depend_immutable
    a[0] = 3.0
    a[1:3] += 1.0
    a[xp.asarray([True, False, False, True])] *= 2.0
    # masks select a copy, as for the arrays themselves
    assert not isinstance(a[xp.asarray([True, False, False, True])], dp.depend_array)
    assert np.array_equal(ab.to_numpy(b.get()), [12.0, 2.0, 2.0, 0.0])
    a.set(1.0)
    assert np.array_equal(ab.to_numpy(b.get()), 2 * np.ones(4))
    with pytest.raises(TypeError):
        a.get()[0] = 0.0


def test_jax_views_follow_parent(jax_backend):
    """Views of depend_arrays holding JAX arrays read the current value of the
    parent, and writing them updates the parent and its dependants"""

    a = dp.depend_array(name="a", value=np.zeros((2, 3)))
    b = dp.depend_array(
        name="b", value=np.zeros((2, 3)), func=lambda: 2 * a.get(), dependencies=[a]
    )
    row, flat = a[1], a[1].flatten()
    a += 1.0
    assert np.array_equal(ab.to_numpy(flat.get()), np.ones(3))
    row += 1.0
    flat[0] = 5.0
    assert np.array_equal(ab.to_numpy(a.get()), [[1.0, 1.0, 1.0], [5.0, 2.0, 2.0]])
    assert np.array_equal(ab.to_numpy(b.get()), 2 * ab.to_numpy(a.get()))
    # a view of one of two synchronized arrays
    c = dp.depend_array(name="c", value=np.zeros((2, 3)))
    dp.dsynchro(c, a, lambda: 2 * a.get(), lambda: c.get() / 2)
    row[2] = 4.0
    assert np.array_equal(ab.to_numpy(c.get()[1]), [10.0, 4.0, 8.0])
    c.set(np.zeros((2, 3)))
    assert np.array_equal(ab.to_numpy(row.get()), np.zeros(3))


def test_jax_depend_copies(jax_backend):
    """depend_arrays holding JAX arrays can be deep-copied and pickled, and their
    value cannot be replaced by assigning it"""

    a = dp.depend_array(name="a", value=np.arange(3.0))
    for b in (deepcopy(a), pickle.loads(pickle.dumps(a))):
        assert type(b) is dp.depend_immutable
        b[0] = 5.0
        assert np.array_equal(ab.to_numpy(b.get()), [5.0, 1.0, 2.0])
    assert np.array_equal(ab.to_numpy(a.get()), [0.0, 1.0, 2.0])
    with pytest.raises(AttributeError):
        a._value = xp.zeros(3)


def test_jax_random_state(jax_backend):
    """With JAX, Gaussian arrays are drawn with a JAX key, which is restored with
    the state of the generator"""

    from ipi.utils.prng import Random

    prng = Random(seed=12345)
    first = prng.gvec((2, 3))
    assert isinstance(first, jax_backend.Array) and first.shape == (2, 3)
    state = prng.get_state()
    second = prng.gvec(6)
    assert not np.allclose(ab.to_numpy(first).reshape(6), ab.to_numpy(second))
    prng.set_state(state)
    assert np.array_equal(ab.to_numpy(prng.gvec(6)), ab.to_numpy(second))

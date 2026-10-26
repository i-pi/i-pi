"""Used to test the depend array view mechanism."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import pickle
from copy import deepcopy

import numpy as np
import pytest

import ipi.engine.atoms
import ipi.utils.array_backend as ab
import ipi.utils.depend as dp


class Holder:
    """Exposes depend members through dproperties."""

    def __init__(self):
        self._x = dp.depend_array(name="x", value=np.zeros(4))
        self._y = dp.depend_array(
            name="y", value=np.zeros(4), func=lambda: 2 * self.x, dependencies=[self._x]
        )


dp.dproperties(Holder, ["x", "y"])


def test_slicing():
    """Depend: indexing a depend_array returns live views or raw scalars"""
    a = dp.depend_array(name="a", value=np.zeros((2, 2), float))
    d1 = dp.depend_array(name="d1", value=np.zeros(6))
    assert isinstance(a[0], dp.depend_array)
    assert isinstance(d1[2:5], dp.depend_array)
    assert isinstance(a[...], dp.depend_array)
    assert isinstance(a[0, :], dp.depend_array)
    assert isinstance(d1[[0, 2, 4]], dp.depend_array)
    assert not isinstance(a[0, 0], dp.depend_array)
    assert not isinstance(d1[2], dp.depend_array)


def test_property_returns_readonly_array():
    """Depend: properties return non-writeable arrays"""
    h = Holder()
    assert type(h.x) is np.ndarray
    with pytest.raises(ValueError):
        h.x[0] = 1.0
    with pytest.raises(ValueError):
        h.x += 1.0


def test_handle_is_not_an_array():
    """Depend: depend objects cannot be used as arrays"""
    h = Holder()
    with pytest.raises(TypeError):
        np.asarray(h._x)


def test_writes_taint_dependants():
    """Depend: writes through the handle, its views and the property taint dependants"""
    h = Holder()
    assert np.array_equal(h.y, np.zeros(4))
    h._x[1] = 1.0
    assert h._y.tainted()
    assert np.array_equal(h.y, [0.0, 2.0, 0.0, 0.0])
    h._x[2:][0] = 3.0
    assert np.array_equal(h.y, [0.0, 2.0, 6.0, 0.0])
    h.x = np.ones(4)
    assert np.array_equal(h.y, 2 * np.ones(4))


def test_inplace_operators():
    """Depend: in-place operators modify the buffer and taint dependants"""
    h = Holder()
    buffer = h.x
    h._x += 1.0
    h._x *= 3.0
    assert np.array_equal(h.y, 6 * np.ones(4))
    # the array is updated in place
    assert np.array_equal(buffer, 3 * np.ones(4))
    h._x[0] -= 3.0
    assert np.array_equal(h.y, [0.0, 6.0, 6.0, 6.0])


def test_indexed_inplace_operators():
    """Depend: `obj._x[i] op= v` updates the array and taints dependants"""
    h = Holder()
    h._x[1:3] += 2.0
    assert np.array_equal(h.y, [0.0, 4.0, 4.0, 0.0])
    a = dp.depend_array(name="a", value=np.zeros((2, 3)))
    a[1] += 1.0
    # a view obtained with a different index is copied
    a[0] = a[1]
    assert np.array_equal(a.get(), np.ones((2, 3)))
    with pytest.raises(NameError):
        h._y[0:2] += 1.0


@pytest.mark.parametrize("backend", ["numpy", "torch"])
@pytest.mark.parametrize("index", [[True, False, True, False], [0, 2]])
def test_advanced_index_inplace_operators(backend, index):
    """Depend: `obj._x[i] op= v` with a mask or an array of indices, which
    select a copy rather than a view, updates the array and taints dependants"""
    if backend == "torch":
        pytest.importorskip("torch")
    ab.set_backend(backend)
    try:
        index = ab.xp.asarray(index)
        h = Holder()
        h._x[index] += 2.0
        assert np.array_equal(ab.to_numpy(h.y), [4.0, 0.0, 4.0, 0.0])
        a = dp.depend_array(name="a", value=ab.xp.ones((2, 4)))
        a[:, index] *= 3.0
        assert np.array_equal(ab.to_numpy(a.get()), [[3.0, 1.0, 3.0, 1.0]] * 2)
    finally:
        ab.set_backend("numpy")


def test_computed_cannot_be_set():
    """Depend: setting a computed quantity raises"""
    h = Holder()
    with pytest.raises(NameError):
        h.y = np.ones(4)


def test_readonly():
    """Depend: assigning to an Atoms property"""
    atoms = ipi.engine.atoms.Atoms(2)
    atoms.q = np.zeros(2 * 3)


def test_reshape_signatures():
    """Depend: reshape accepts both tuple and varargs forms."""
    da = dp.depend_array(name="r", value=np.arange(6.0))
    r1 = da.reshape((2, 3))
    r2 = da.reshape(2, 3)
    assert isinstance(r1, dp.depend_array) and r1.get().shape == (2, 3)
    assert isinstance(r2, dp.depend_array) and r2.get().shape == (2, 3)
    # Reshaped view shares memory with the parent.
    r1[0, 0] = 99.0
    assert da[0] == 99.0


def test_pickle_roundtrip():
    """Depend: pickle/unpickle preserves values and recreates the lock."""
    da = dp.depend_array(name="p", value=np.array([1.0, 2.0, 3.0]))
    rt = pickle.loads(pickle.dumps(da))
    assert isinstance(rt, dp.depend_array)
    assert np.array_equal(rt.get(), da.get())
    with rt._status.lock:
        pass
    rt[0] = 5.0
    assert rt[0] == 5.0


def test_deepcopy():
    """Depend: deep copies have independent storage and status."""
    da = dp.depend_array(name="c", value=np.array([1.0, 2.0, 3.0]))
    cp = deepcopy(da)
    cp[0] = 5.0
    assert cp[0] == 5.0
    assert da[0] == 1.0
    assert cp._status is not da._status


def test_strided_slice_refresh_via_parent():
    """Depend: a strided 1D view of an auto-computed array refreshes
    through its parent, so first access through the view does not crash
    with a shape mismatch."""
    natoms = 4
    counter = [0]

    def get_full():
        counter[0] += 1
        return np.arange(3 * natoms, dtype=float) + 100.0 * counter[0]

    trigger = dp.depend_value(name="trig", value=0)
    full = dp.depend_array(
        name="full", value=np.zeros(3 * natoms), func=get_full, dependencies=[trigger]
    )
    sx = full._view(full._value[0::3], full._wvalue[0::3])

    # first access through the view delegates refresh to the parent
    sx_first = sx.get().copy()
    assert counter[0] == 1
    assert np.array_equal(sx_first, full.get()[0::3])

    # tainting `trigger` taints the view through the shared status
    trigger.taint()
    assert sx.tainted()
    sx_new = sx.get().copy()
    assert counter[0] == 2
    assert np.array_equal(sx_new, full.get()[0::3])


def test_synchronizer():
    """Depend: writing a synchronized peer makes the other stale"""
    a = dp.depend_array(name="a", value=np.ones(3))
    b = dp.depend_array(name="b", value=np.zeros(3))
    b_sum = dp.depend_value(name="bsum", func=lambda: b.get().sum(), dependencies=[b])
    dp.dsynchro(a, b, lambda: 0.5 * b.get(), lambda: 2 * a.get())
    # `a` holds the current data when the two are linked
    assert not a.tainted() and b.tainted()
    assert np.array_equal(b.get(), 2 * np.ones(3))
    b[:] = 4.0
    assert a.tainted()
    assert np.array_equal(a.get(), 2 * np.ones(3))
    # a view shares the synchronizer of its parent
    a[1:][0] = 0.0
    assert b.tainted() and b_sum.tainted()
    assert np.array_equal(b.get(), [4.0, 0.0, 4.0])
    assert b_sum.get() == 8.0
    # in-place updates also propagate
    b += 1.0
    assert np.array_equal(a.get(), [2.5, 0.5, 2.5])
    with pytest.raises(AssertionError):
        dp.dsynchro(a, dp.depend_array(name="c", value=np.zeros(3)), None, None)

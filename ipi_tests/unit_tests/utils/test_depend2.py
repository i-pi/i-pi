"""Used to test the depend array view mechanism."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np

import ipi.engine.atoms
import ipi.utils.depend as dp
from ipi.utils.array_backend import xp

# Test arrays are built with `xp` so their storage matches the active
# array backend (numpy when IPI_ARRAY_BACKEND is unset, torch when set).
# Using `np.zeros(...)` directly would mix numpy storage with a torch-
# backed depend_array and fail at assignment time.
a = dp.depend_array(name="a", value=xp.zeros((2, 2)))
b = dp.depend_array(name="b", value=xp.zeros((2, 2)))


def test_slicing():
    """Depend: __getitem__ returns raw backend views (no depend wrapping)."""
    raw_t = type(xp.zeros((1,)))
    # Integer index into a 2D array returns a backend row view.
    c = a[0]
    assert isinstance(c, raw_t)
    assert not isinstance(c, dp.depend_array)
    # 1D slice returns a backend view, not a depend_array.
    d1 = dp.depend_array(name="d1", value=xp.zeros(6))
    s = d1[2:5]
    assert isinstance(s, raw_t)
    assert not isinstance(s, dp.depend_array)
    # Ellipsis and tuple-with-slice return backend views.
    assert isinstance(a[...], raw_t) and not isinstance(a[...], dp.depend_array)
    assert isinstance(a[0, :], raw_t) and not isinstance(a[0, :], dp.depend_array)
    # Full multidim integer index returns a backend scalar.
    assert not isinstance(a[0, 0], dp.depend_array)
    # Integer index on a 1D array returns a backend scalar.
    assert not isinstance(d1[2], dp.depend_array)


def test_dslice():
    """Depend: dslice returns a depend-aware slice that shares _status."""
    parent = dp.depend_array(name="parent", value=xp.zeros((4, 3)))
    child = parent.dslice((1, slice(None)))
    assert isinstance(child, dp.depend_array)
    # Status is shared by reference (taint propagation works both ways).
    assert child._status is parent._status
    # Refresh delegation: child knows its parent.
    assert child._parent is parent


def test_value():
    """Depend: .value returns the raw backend tensor and triggers refresh."""
    raw_t = type(xp.zeros((1,)))
    da = dp.depend_array(name="x", value=xp.zeros((3, 2)))
    v = da.value
    assert isinstance(v, raw_t)
    # Mutating through the parent's __setitem__ retaints downstream
    # consumers; .value reflects the latest data on next read.
    da[0, 0] = 1.0
    v2 = da.value
    assert float(v2[0, 0]) == 1.0


def test_addition():
    """Depend: Addition test"""
    c = a + b
    print((type(c)))
    # `a + b` returns the backend-native array type (ndarray under numpy,
    # torch.Tensor under torch).
    assert isinstance(c, type(xp.zeros((1,))))


def test_increment():
    """Depend: Increment test"""
    c = xp.zeros((2, 2))
    c += a
    print((type(c)))
    assert isinstance(c, type(xp.zeros((1,))))


def test_dot():
    """Depend: Dot test"""
    # Under torch, aten ops reject depend_array wrappers — extract raw via .value first.
    c = xp.tensordot(a.value, b.value, axes=([-1], [-2]))
    print((type(c)))
    assert isinstance(c, type(xp.zeros((1,))))


def test_dotf():
    """Depend: `@` on depend_array operands returns a raw backend array."""
    c = a @ b
    assert isinstance(c, type(xp.zeros((1,))))


def test_readonly():
    """Depend: read-only flag"""
    atoms = ipi.engine.atoms.Atoms(2)
    atoms.q = xp.zeros(2 * 3)

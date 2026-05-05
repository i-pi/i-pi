"""Used to test the depend array view mechanism."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import pickle

import numpy as np

import ipi.engine.atoms
import ipi.utils.depend as dp


a = dp.depend_array(name="a", value=np.zeros((2, 2), float))
b = dp.depend_array(name="b", value=np.zeros((2, 2), float))


def test_slicing():
    """Depend: Slicing test"""
    # Integer index into a 2D array returns a row (sub-array), not a scalar
    c = a[0]
    print((type(c)))
    assert isinstance(c, dp.depend_array)
    # Slice of a 1D array must also return a depend_array, not a plain ndarray
    d1 = dp.depend_array(name="d1", value=np.zeros(6))
    s = d1[2:5]
    assert isinstance(
        s, dp.depend_array
    ), f"1D slice returned {type(s)}, expected depend_array"
    # Ellipsis on a 2D array returns a depend_array view, not a scalar
    assert isinstance(a[...], dp.depend_array)
    # Tuple-with-slice on a 2D array returns a depend_array view
    assert isinstance(a[0, :], dp.depend_array)
    # Fancy indexing on a 1D array returns a depend_array
    assert isinstance(d1[[0, 2, 4]], dp.depend_array)
    # Full multidim integer index on a 2D array returns a raw scalar
    assert not isinstance(a[0, 0], dp.depend_array)
    # Integer index on a 1D array returns a raw scalar
    assert not isinstance(d1[2], dp.depend_array)


def test_addition():
    """Depend: Addition test"""
    c = a + b
    print((type(c)))
    assert isinstance(c, np.ndarray)


def test_increment():
    """Depend: Increment test"""
    c = np.zeros((2, 2))
    c += a
    print((type(c)))
    assert isinstance(c, np.ndarray)


def test_dot():
    """Depend: Dot test"""
    c = np.dot(a, b)
    print((type(c)))
    assert isinstance(c, np.ndarray)


def test_dotf():
    """Depend: Dot-f test"""

    rdot = np.dot

    def fdot(aa, bb):
        return rdot(aa, bb).view(np.ndarray)

    np.dot = fdot

    c = np.dot(a, b)
    assert isinstance(c, np.ndarray)
    np.dot = rdot


def test_readonly():
    """Depend: read-only flag"""
    atoms = ipi.engine.atoms.Atoms(2)
    atoms.q = np.zeros(2 * 3)


def test_reshape_signatures():
    """Depend: reshape accepts both tuple and varargs forms."""
    da = dp.depend_array(name="r", value=np.arange(6.0))
    r1 = da.reshape((2, 3))
    r2 = da.reshape(2, 3)
    assert isinstance(r1, dp.depend_array) and r1.shape == (2, 3)
    assert isinstance(r2, dp.depend_array) and r2.shape == (2, 3)
    # Reshaped view shares memory with the parent.
    r1[0, 0] = 99.0
    assert da[0] == 99.0


def test_pickle_roundtrip():
    """Depend: pickle/unpickle preserves values and drops the threadlock cleanly."""
    da = dp.depend_array(name="p", value=np.array([1.0, 2.0, 3.0]))
    blob = pickle.dumps(da)
    rt = pickle.loads(blob)
    assert isinstance(rt, dp.depend_array)
    assert np.array_equal(np.asarray(rt), np.asarray(da))
    # Threadlock is recreated on unpickle.
    with rt._threadlock:
        pass


def test_asarray_zerocopy():
    """Depend: np.asarray on a depend_array exposes the underlying buffer."""
    da = dp.depend_array(name="z", value=np.zeros(4))
    arr = np.asarray(da)
    arr[1] = 7.0
    assert da[1] == 7.0


def test_slice_writes_taint_parent_dependant():
    """Depend: writing through a slice taints downstream consumers of the parent."""
    parent = dp.depend_array(name="par", value=np.zeros(4))
    child = dp.depend_array(
        name="ch", value=np.zeros(4), func=lambda: 2 * np.asarray(parent)
    )
    child.add_dependency(parent)
    # Force a refresh so child becomes untainted.
    _ = np.asarray(child)
    assert not child.tainted()
    # Writing through a slice of the parent must taint the child.
    parent[1:3] = np.array([5.0, 6.0])
    assert child.tainted()
    # And the recompute uses the updated parent values.
    assert np.array_equal(np.asarray(child), np.array([0.0, 10.0, 12.0, 0.0]))


def test_iteration():
    """Depend: iterating over a 2D depend_array yields the rows."""
    da = dp.depend_array(name="it", value=np.arange(6.0).reshape(2, 3))
    rows = list(da)
    assert len(rows) == 2
    assert np.array_equal(rows[0], np.array([0.0, 1.0, 2.0]))
    assert np.array_equal(rows[1], np.array([3.0, 4.0, 5.0]))


def test_strided_slice_refresh_via_parent():
    """Depend: a strided 1D slice of an auto-computed array refreshes
    through its parent (not through its own `_func`), so first-access
    through the slice does not crash with a shape mismatch.

    This is the pattern that was previously used for ForceBead's per-axis
    force shortcuts: a parent depend_array with a `_func` returning the
    full-shape buffer, plus strided 1D slices that share `_tainted` and
    delegate refresh upward via `_parent`.
    """
    natoms = 4
    fbase = np.zeros(3 * natoms)
    counter = [0]

    def get_full():
        counter[0] += 1
        return np.arange(3 * natoms, dtype=float) + 100.0 * counter[0]

    trigger = dp.depend_value(name="trig", value=0)
    full = dp.depend_array(
        name="full", value=fbase, func=get_full, dependencies=[trigger]
    )
    # Strided 1D view of the parent's buffer, with refresh delegation.
    sx = dp.depend_array(name="sx", value=fbase[0 : 3 * natoms : 3], parent=full)
    dp.dcopy(full, sx)

    # First-access through the slice (before `full` is touched) delegates
    # refresh to the parent — no shape mismatch.
    sx_first = np.asarray(sx).copy()
    assert counter[0] == 1
    assert np.array_equal(sx_first, np.asarray(full)[0::3])

    # Tainting `trigger` taints the slice via the shared `_tainted`.
    trigger.taint()
    assert sx.tainted()
    sx_new = np.asarray(sx).copy()
    assert counter[0] == 2
    assert np.array_equal(sx_new, np.asarray(full)[0::3])

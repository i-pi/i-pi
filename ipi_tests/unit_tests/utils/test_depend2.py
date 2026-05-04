"""Used to test the depend array view mechanism."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


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
    assert isinstance(s, dp.depend_array), (
        f"1D slice returned {type(s)}, expected depend_array"
    )


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

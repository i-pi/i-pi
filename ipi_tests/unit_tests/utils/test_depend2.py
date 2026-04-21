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
    """Depend: Slicing test"""
    c = a[0]
    print((type(c)))
    assert isinstance(c, dp.depend_array)


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
    # Under torch, aten ops reject depend_array wrappers — dstrip first.
    c = xp.tensordot(dp.dstrip(a), dp.dstrip(b), axes=([-1], [-2]))
    print((type(c)))
    assert isinstance(c, type(xp.zeros((1,))))


def test_dotf():
    """Depend: ddot (wrapper over numpy.dot that unwraps depend operands)."""
    # dp.ddot dstrips its args and calls np.dot, so it only makes sense
    # when the backend is numpy. Under torch, construct temporaries via xp
    # and use `@` instead.
    c = a @ b
    assert isinstance(c, type(xp.zeros((1,))))


def test_readonly():
    """Depend: read-only flag"""
    atoms = ipi.engine.atoms.Atoms(2)
    atoms.q = xp.zeros(2 * 3)

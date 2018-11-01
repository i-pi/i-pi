"""Tests ring polymer contraction."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np
from numpy.testing import assert_almost_equal as assert_equals

from ipi.utils import nmtransform


def check_up_and_down_scaling(n, q):
    """Check if q expanding and then contracting a ring polymer is a no-op.

    Args:
       n: The number of beads in the scaled ring polymer.
       q: The original position array.
    """

    rescale = nmtransform.nm_rescale(q.shape[0], n)
    print "Initial position of the beads:"
    print q, q.shape, (q.shape[0], n)

    # rescale up to the n beads
    beads_n = rescale.b1tob2(q)
    print "Upscaled to %d beads:" % n
    print beads_n, beads_n.shape

    beads_final = rescale.b2tob1(beads_n)
    print "Final position of the beads:"
    print beads_final

    assert_equals(q, beads_final)
    return beads_n


def check_rpc_consistency(n, q):
    """Check if q expanding and then contracting a ring polymer is a no-op.

    Args:
       n: The number of beads in the scaled ring polymer.
       q: The original position array.
    """

    rescale1 = nmtransform.nm_rescale(q.shape[0], n)
    rescale2 = nmtransform.nm_rescale(n, q.shape[0])

    beads_n = rescale1.b1tob2(q)
    beads_1 = rescale1.b2tob1(beads_n)
    beads_2 = rescale2.b1tob2(beads_n)

    assert_equals(beads_1, beads_2)


def check_centroid_pos(n, q):
    """Check if expanding and then contracting a ring polymer
    maintains the centroid.

    Args:
       n: The number of beads in the scaled ring polymer.
       q: The original position array.
    """

    beads_big = check_up_and_down_scaling(n, q)
    rescale_big = nmtransform.mk_rs_matrix(n, 1)
    rescale_q = nmtransform.mk_rs_matrix(q.shape[0], 1)

    centroid_big = np.dot(rescale_big, beads_big)
    centroid_q = np.dot(rescale_q, q)

    assert_equals(centroid_q, centroid_big)


numbers_to_check = range(10, 56, 9)


def test_1_to_n():
    """Contraction with one bead."""

    for n in numbers_to_check:
        q = np.array([[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]])
        yield check_up_and_down_scaling, n, q
        yield check_rpc_consistency, n, q
        yield check_centroid_pos, n, q


def test_2_to_n():
    """Contraction with two beads."""

    for n in numbers_to_check:
        q = np.array([[0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                      [0.0, 0.1, 0.0, 1.0, 0.1, 0.0]])
        yield check_up_and_down_scaling, n, q
        yield check_rpc_consistency, n, q
        yield check_centroid_pos, n, q


def test_3_to_n():
    """Contraction with tree beads."""

    for n in numbers_to_check:
        q = np.array([[0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                      [0.0, 0.1, 0.0, 1.0, 0.1, 0.0],
                      [0.0, -0.1, 0.0, 1.0, -0.1, 0.0]])
        yield check_up_and_down_scaling, n, q
        yield check_rpc_consistency, n, q
        yield check_centroid_pos, n, q


def test_4_to_n():
    """Contraction with four beads."""

    for n in numbers_to_check:
        q = np.array([[0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                      [0.0, 0.1, 0.0, 1.0, 0.1, 0.0],
                      [0.0, 0.2, 0.0, 1.0, 0.2, 0.0],
                      [0.0, -0.1, 0.0, 1.0, -0.1, 0.0]])
        yield check_up_and_down_scaling, n, q
        yield check_rpc_consistency, n, q
        yield check_centroid_pos, n, q

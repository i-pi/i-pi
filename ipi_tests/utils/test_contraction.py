"""Tests ring polymer contraction."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import pytest
import numpy as np
from numpy.testing import assert_almost_equal

from ipi.utils import nmtransform


def check_up_and_down_scaling(n, q):
    """Check if q expanding and then contracting a ring polymer is a no-op.

    Args:
       n: The number of beads in the scaled ring polymer.
       q: The original position array.
    """

    rescale = nmtransform.nm_rescale(q.shape[0], n)
    print("Initial position of the beads:")
    print(q, q.shape, (q.shape[0], n))

    # rescale up to the n beads
    beads_n = rescale.b1tob2(q)
    print("Upscaled to %d beads:" % n)
    print(beads_n, beads_n.shape)

    beads_final = rescale.b2tob1(beads_n)
    print("Final position of the beads:")
    print(beads_final)

    assert_almost_equal(q, beads_final)
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

    assert_almost_equal(beads_1, beads_2)


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

    assert_almost_equal(centroid_q, centroid_big)


numbers_to_check = list(range(10, 56, 9))


qs_to_check = [
    np.array([[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]]),
    np.array([[0.0, 0.0, 0.0, 1.0, 0.0, 0.0], [0.0, 0.1, 0.0, 1.0, 0.1, 0.0]]),
    np.array(
        [
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.1, 0.0, 1.0, 0.1, 0.0],
            [0.0, -0.1, 0.0, 1.0, -0.1, 0.0],
        ]
    ),
    np.array(
        [
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.1, 0.0, 1.0, 0.1, 0.0],
            [0.0, 0.2, 0.0, 1.0, 0.2, 0.0],
            [0.0, -0.1, 0.0, 1.0, -0.1, 0.0],
        ]
    ),
]


@pytest.mark.parametrize("n", numbers_to_check)
@pytest.mark.parametrize("q", qs_to_check)
def test_contraction(n, q):
    """Test contraction to `n` replicas with original positions `q`."""

    check_up_and_down_scaling(n, q)
    check_rpc_consistency(n, q)
    check_centroid_pos(n, q)

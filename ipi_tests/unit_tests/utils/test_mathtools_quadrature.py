"""Unit tests for the Legendre roots and the rotation quadrature, which must be real."""

import numpy as np

from ipi.utils.mathtools import (
    roots_legendre,
    get_rotation_quadrature_legendre,
)


def test_roots_legendre_real():
    # legroots returns complex on numpy>=2.5; roots/weights must stay real
    roots, weights = roots_legendre(2)
    assert roots.dtype.kind == "f"
    assert weights.dtype.kind == "f"
    assert np.allclose(sorted(roots), [-0.5773502692, 0.5773502692])


def test_rotation_quadrature_legendre_real():
    # complex roots would propagate to complex rotation matrices and corrupt
    # the data sent to the driver over the socket
    for mat, weight, angles in get_rotation_quadrature_legendre(2):
        assert mat.dtype.kind == "f"
        assert np.isrealobj(np.asarray(angles))

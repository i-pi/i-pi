"""Unit tests for the sinh(x)/x helper, including complex arguments.

Complex eigenvalues from a non-symmetric barostat momentum tensor reach _sinch
(via np.linalg.eig), so the small-argument test must use abs(); these checks
guard that path together with the real-valued behaviour.
"""

import numpy as np

from ipi.utils.mathtools import (
    _sinch,
    sinch,
    roots_legendre,
    get_rotation_quadrature_legendre,
)


def test_sinch_real_matches_definition():
    for x in (0.5, 2.0, -1.3):
        assert np.isclose(_sinch(x), np.sinh(x) / x)


def test_sinch_small_real_uses_taylor():
    # near zero sinh(x)/x -> 1 (the direct ratio would be 0/0)
    assert np.isfinite(_sinch(1e-9))
    assert np.isclose(_sinch(1e-9), 1.0)


def test_sinch_complex():
    # a genuinely complex argument (would raise on the old `x2 < 1e-12` test)
    for x in (1j, 0.3 + 0.7j, 1e-9j):
        assert np.isclose(_sinch(x), np.sinh(x) / x if abs(x) > 1e-6 else 1.0)


def test_sinch_vectorized_over_complex_eigenvalues():
    v = np.array([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 0.1]])
    eigvals = np.linalg.eig(v)[0]
    assert np.iscomplexobj(eigvals)
    out = sinch(0.5 * eigvals)
    assert np.allclose(out, np.sinh(0.5 * eigvals) / (0.5 * eigvals))


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

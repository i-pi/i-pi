"""Identities for matrix functions and signed log sums used by dynamics."""

import math

import numpy as np
import pytest

from ipi.utils.mathtools import (
    logsumlog,
    matrix_exp,
    root_herm,
    stab_cholesky,
)


@pytest.mark.parametrize(
    "a,b",
    [
        (2.0, 3.0),
        (-2.0, 3.0),
        (2.0, -3.0),
        (-2.0, -3.0),
        (3.0, 2.0),
        (-3.0, 2.0),
        (3.0, -2.0),
        (-3.0, -2.0),
    ],
)
def test_signed_log_sum_matches_direct_addition(a, b):
    """All sign and magnitude orderings must represent the ordinary sum."""
    magnitude, sign = logsumlog(
        (math.log(abs(a)), np.sign(a)), (math.log(abs(b)), np.sign(b))
    )
    assert sign * math.exp(magnitude) == pytest.approx(a + b)


@pytest.mark.parametrize("diagonal", [0.5, -0.5])
def test_matrix_exponential_of_jordan_block(diagonal):
    """For N cubed equal to zero, exp(a I + N) = exp(a) (I + N + N²/2)."""
    nilpotent = np.array([[0.0, 2.0, 0.0], [0.0, 0.0, 3.0], [0.0, 0.0, 0.0]])
    matrix = diagonal * np.eye(3) + nilpotent
    expected = math.exp(diagonal) * (np.eye(3) + nilpotent + nilpotent @ nilpotent / 2)
    np.testing.assert_allclose(matrix_exp(matrix), expected, rtol=2e-12, atol=2e-12)


@pytest.mark.parametrize("factor", [root_herm, stab_cholesky])
@pytest.mark.parametrize("singular", [False, True])
def test_covariance_factor_reconstructs_positive_matrix(factor, singular):
    """A covariance factor must reproduce correlated noise, including a zero mode."""
    basis = np.array(
        [[1.0, 0.0, 0.0], [2.0, 3.0, 0.0], [-1.0, 0.5, 0.0 if singular else 2.0]]
    )
    covariance = basis @ basis.T
    result = factor(covariance)
    np.testing.assert_allclose(result @ result.T, covariance, rtol=1e-12, atol=1e-12)


def test_symmetric_root_rejects_nonsymmetric_input():
    """An invalid covariance matrix must not silently produce a noise factor."""
    with pytest.raises(ValueError, match="Non-Hermitian"):
        root_herm(np.array([[1.0, 2.0], [0.0, 1.0]]))

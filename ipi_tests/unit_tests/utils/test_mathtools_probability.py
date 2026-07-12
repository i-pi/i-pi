"""Tests for probability-related mathematical helpers."""

import pytest

from ipi.utils.mathtools import gaussian_inv


@pytest.mark.parametrize(
    "probability, expected",
    [
        (0.001, -3.0902323062),
        (0.025, -1.9599639845),
        (0.5, 0.0),
        (0.975, 1.9599639845),
        (0.999, 3.0902323062),
    ],
)
def test_gaussian_inv_matches_known_quantiles(probability, expected):
    """Checks central and tail values of the inverse normal distribution."""

    assert gaussian_inv(probability) == pytest.approx(expected, abs=1.0e-6)

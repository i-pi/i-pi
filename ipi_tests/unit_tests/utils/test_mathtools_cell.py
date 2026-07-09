"""Tests for the cell-geometry helpers in ipi.utils.mathtools."""

import math

import numpy as np
import pytest

from ipi.utils.mathtools import abc2h, h2abc, det_ut3x3, invert_ut3x3

# (a, b, c, alpha, beta, gamma) in radians
CELLS = [
    (1.0, 1.0, 1.0, math.pi / 2, math.pi / 2, math.pi / 2),  # cubic
    (3.0, 4.0, 5.0, math.pi / 2, math.pi / 2, math.pi / 2),  # orthorhombic
    (2.0, 2.0, 3.0, math.pi / 2, math.pi / 2, 2 * math.pi / 3),  # hexagonal-ish
    (2.5, 3.5, 4.5, 1.2, 1.3, 1.1),  # triclinic
]


@pytest.mark.parametrize("abc", CELLS)
def test_abc2h_h2abc_roundtrip(abc):
    h = abc2h(*abc)
    assert np.allclose(np.array(h2abc(h)), np.array(abc))


@pytest.mark.parametrize("abc", CELLS)
def test_det_equals_diagonal_product(abc):
    h = abc2h(*abc)
    assert det_ut3x3(h) == pytest.approx(np.linalg.det(h))


@pytest.mark.parametrize("abc", CELLS)
def test_invert_ut3x3_matches_numpy(abc):
    h = abc2h(*abc)
    assert np.allclose(invert_ut3x3(h), np.linalg.inv(h))

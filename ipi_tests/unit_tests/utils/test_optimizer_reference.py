"""Optimization steps checked by stationarity, symmetry and known quadratics."""

import numpy as np
import pytest

from ipi.utils.array_backend import xp
from ipi.utils.mintools import (
    Damped_BFGS,
    FIRE,
    _max_atomic_displacement,
    min_trm,
    nichols,
    Powell,
)
from ipi.utils.softexit import softexit


def test_stuck_fire_requests_a_soft_exit():
    """FIRE stops with a soft exit when it has gone uphill for more than Nmax steps"""

    softexit.reset()
    try:
        # the force opposes the velocity for more than Nmax steps
        with pytest.raises(SystemExit):
            FIRE(xp.zeros(3), None, (0.0, xp.ones(3)), v=xp.ones(3), N_up=5, Nmax=5)
        assert softexit.triggered
    finally:
        softexit.reset()


@pytest.mark.parametrize("radius", [0.1, 10.0])
@pytest.mark.parametrize("eigenvalues", [[1.0, 2.0, 4.0], [-1.0, 2.0, 4.0]])
def test_trust_region_step_satisfies_stationarity_and_radius(radius, eigenvalues):
    h = np.diag(eigenvalues)
    force = np.array([0.2, 0.4, -0.3])
    step = min_trm(force, h, radius)
    assert np.linalg.norm(step) <= radius * (1 + 1e-4)
    residual = force - h @ step
    multiplier = np.dot(residual, step) / np.dot(step, step)
    np.testing.assert_allclose(residual, multiplier * step, atol=1e-7)
    assert multiplier >= -1e-7
    assert np.min(np.linalg.eigvalsh(h + multiplier * np.eye(3))) >= -1e-7
    assert 0.5 * step @ h @ step - force @ step <= 0


@pytest.mark.parametrize(
    "diagonal", [[1.0, 2.0, 3.0], [-1.0, 2.0, 3.0], [-3.0, -2.0, 1.0]]
)
def test_nichols_step_is_covariant_under_rotation(diagonal):
    rotation = np.linalg.qr(
        np.array([[1.0, 2.0, 3.0], [3.0, 1.0, 2.0], [2.0, 4.0, 1.0]])
    )[0]
    force = np.array([0.1, -0.3, 0.2])
    d = np.array(diagonal)
    base = nichols(force, d, np.eye(3), np.ones(3), 0.5)
    rotated = nichols(rotation @ force, d, rotation, np.ones(3), 0.5)
    np.testing.assert_allclose(rotated, rotation @ base, atol=1e-12)


def test_powell_update_satisfies_secant_equation_and_symmetry():
    h = np.eye(3)
    step = np.array([0.2, 0.1, -0.3])
    change = np.array([0.5, -0.2, 0.8])
    updated = Powell(step, change, h)
    np.testing.assert_allclose(updated @ step, change, atol=1e-12)
    np.testing.assert_allclose(updated, updated.T, atol=1e-12)


def test_damped_bfgs_step_is_measured_on_atoms():
    """With fixed coordinates, the step of Damped_BFGS is limited by the
    displacements of the atoms, not of triplets of free coordinates, also
    with the free coordinates of several beads"""

    # z of the first atom, x and y of the second
    mask = np.array([False, False, True, True, True, False])
    assert _max_atomic_displacement(np.array([3.0, 4.0, 0.0]), mask) == 4.0
    assert _max_atomic_displacement(np.array([3.0, 4.0, 0.0]), np.ones(3, bool)) == 5.0
    dx = np.array([3.0, 4.0, 0.0, 1.0, 0.0, 2.0])
    assert _max_atomic_displacement(dx, mask) == 4.0


def test_damped_bfgs_at_a_stationary_point():
    """At a zero gradient, Damped_BFGS stays at the same point and leaves the
    Hessian unchanged"""

    positions = []

    def fdf(x):
        positions.append(np.array(x))
        return 0.0, np.zeros(3)

    x0 = np.array([[0.1, 0.2, 0.3]])
    hessian = 2 * np.eye(3)
    quality = Damped_BFGS(x0, fdf, (0.0, np.zeros(3)), hessian, 0.5, np.ones(3, bool))
    assert np.isfinite(quality)
    np.testing.assert_array_equal(hessian, 2 * np.eye(3))
    np.testing.assert_array_equal(positions[-1], x0)

"""Optimization steps checked by stationarity, symmetry and known quadratics."""

import numpy as np
import pytest

from ipi.utils.mintools import L_BFGS_nls, min_brent_neb, min_trm, nichols, Powell


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


@pytest.mark.parametrize("center", [-2.0, 1.5])
@pytest.mark.parametrize("quartic", [False, True])
def test_neb_scalar_search_finds_known_minimum(center, quartic):

    def objective(x):
        d = x - center
        return 0.0, d * d + (0.2 * d**4 if quartic else 0.0)

    position, value = min_brent_neb(
        objective, fdf0=objective(0.0)[1], x0=0.0, tol=1e-9, itmax=200, init_step=0.1
    )
    assert position == pytest.approx(center, abs=1e-6)
    assert value < 1e-12


def test_lbfgs_without_line_search_on_quadratic():

    h = np.diag([1.0, 2.0, 3.0])
    x = np.array([[0.2, -0.3, 0.4]])

    def objective(q):
        return 0.5 * q.ravel() @ h @ q.ravel(), (h @ q.ravel()).reshape(q.shape)

    direction = -objective(x)[1]
    qlist = np.zeros((3, 3))
    glist = np.zeros((3, 3))
    for k in range(30):
        if np.linalg.norm(objective(x)[1]) < 1e-8:
            break
        x, value, direction, qlist, glist = L_BFGS_nls(
            x,
            direction,
            objective,
            qlist,
            glist,
            big_step=0.2,
            init_step=0.05,
            m=3,
            k=k,
        )
        assert value == pytest.approx(objective(x)[0])
    np.testing.assert_allclose(x, 0, atol=1e-6)

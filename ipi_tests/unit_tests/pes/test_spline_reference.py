"""Spline driver evaluated on a tabulated quadratic."""

import json

import numpy as np
import pytest

from ipi.pes import spline


@pytest.fixture
def quadratic(tmp_path):
    x = np.linspace(-2, 2, 21)
    path = tmp_path / "quadratic.dat"
    np.savetxt(path, np.column_stack([x, -2 * x, x * 0, x * 0, x * x]))
    return spline.Spline_driver(str(path))


def test_tabulated_quadratic_energy_and_transverse_forces(quadratic):
    pos = np.array([[0.3, -0.2, 0.1]])
    energy, force, virial, extra = quadratic.compute_structure(np.eye(3), pos)
    assert energy == pytest.approx(
        pos[0, 0] ** 2 + quadratic.k / 2 * np.sum(pos[0, 1:] ** 2)
    )
    np.testing.assert_allclose(force[0, 1:], -quadratic.k * pos[0, 1:])
    np.testing.assert_array_equal(virial, np.zeros((3, 3)))
    assert extra == ""


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="Spline force omits the coordinate-scale factor in the energy derivative",
)
def test_tabulated_quadratic_longitudinal_force(quadratic):
    pos = np.array([[0.3, 0.0, 0.0]])
    assert quadratic.get_forces(pos)[0, 0] == pytest.approx(-2 * pos[0, 0], abs=1e-10)


def test_constant_friction_is_serialized(quadratic, monkeypatch):
    monkeypatch.setattr(spline, "friction", True)
    _, _, _, extra = quadratic.compute_structure(np.eye(3), np.zeros((1, 3)))
    np.testing.assert_allclose(
        json.loads(extra)["friction"], np.eye(3) * spline.fric_value
    )

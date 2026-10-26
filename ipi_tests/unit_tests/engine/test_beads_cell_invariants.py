"""Physical invariants for bead state, spring forces, and periodic cells."""

import numpy as np
import pytest

from ipi.engine.atoms import Atoms
from ipi.engine.beads import Beads
from ipi.engine.cell import Cell, GenericCell


@pytest.fixture
def beads():
    """A nonuniform ring polymer with unequal atomic masses."""
    result = Beads(2, 4)
    result._names[:] = ["H", "O"]
    result._m[:] = [2.0, 5.0]
    result._q[:] = np.random.default_rng(82).normal(size=(4, 6))
    result._p[:] = np.arange(24).reshape(4, 6) / 10
    return result


def test_spring_force_is_negative_energy_gradient(beads):
    """A displacement's energy change and force must agree for every bead."""
    original = beads.q.copy()
    force = beads.fpath.copy()
    numerical = np.zeros_like(original)
    step = 1e-5
    for index in np.ndindex(original.shape):
        beads._q[index] = original[index] + step
        forward = beads.vpath
        beads._q[index] = original[index] - step
        backward = beads.vpath
        beads._q[index] = original[index]
        numerical[index] = -(forward - backward) / (2 * step)
    np.testing.assert_allclose(force, numerical, rtol=1e-8, atol=1e-8)
    np.testing.assert_allclose(force.sum(axis=0), 0, atol=1e-12)


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="Beads.fpath is not invalidated when masses change",
)
def test_spring_force_updates_after_mass_change(beads):
    """Changing isotope masses must invalidate already evaluated spring forces."""
    original_force = beads.fpath.copy()
    beads._m[:] *= 2
    np.testing.assert_allclose(beads.fpath, 2 * original_force, atol=1e-12)


def test_bead_views_and_clones_preserve_state_without_aliasing(beads):
    """Per-bead edits update the centroid; cloning produces independent storage."""
    original_centroid = beads.qc.copy()
    beads[1]._q[:] += 4
    np.testing.assert_allclose(beads.qc, original_centroid + 1)
    clone = beads.clone(2)
    np.testing.assert_array_equal(clone.q, beads.q[:2])
    np.testing.assert_array_equal(clone.p, beads.p[:2])
    np.testing.assert_array_equal(clone.names, beads.names)
    clone._q[:] = -100
    clone._m[:] = 9
    assert np.all(beads.q != -100)
    np.testing.assert_array_equal(beads.m, [2, 5])
    with pytest.raises(ValueError, match="larger number of beads"):
        beads.clone(5)


@pytest.mark.parametrize("cell_type", [Cell, GenericCell])
def test_periodic_wrapping_preserves_fractional_coordinates(cell_type):
    """Lattice translations leave the wrapped position and pair separation unchanged."""
    h = np.array([[3.0, 0.3, -0.2], [0.0, 4.0, 0.4], [0.0, 0.0, 5.0]])
    if cell_type is GenericCell:
        h = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]) @ h
    cell = cell_type(h.copy())
    fractional = np.array([[0.1, -0.2, 0.3], [-0.2, 0.1, -0.1]])
    shifts = np.array([[3, -2, 1], [-1, 4, -3]])
    atoms = Atoms(2)
    atoms._q[:] = ((fractional + shifts) @ h.T).ravel()
    expected = fractional @ h.T
    np.testing.assert_allclose(cell.apply_pbc(atoms[0]), expected[0], atol=1e-12)
    np.testing.assert_allclose(
        cell.minimum_distance(atoms[0], atoms[1]), expected[0] - expected[1], atol=1e-12
    )
    q = cell.array_pbc(np.array(atoms.q))
    np.testing.assert_allclose(q.reshape(-1, 3), expected, atol=1e-12)
    clone = cell.clone()
    old_volume = cell.V
    cell._h[:] *= 2
    assert cell.V == pytest.approx(8 * old_volume)
    np.testing.assert_allclose(cell.ih @ cell.h, np.eye(3), atol=1e-12)
    np.testing.assert_allclose(clone.h, h, atol=1e-12)


@pytest.mark.xfail(
    strict=True,
    raises=AttributeError,
    reason="Beads total kinetic energy refers to missing kins attribute",
)
def test_total_bead_kinetic_energy_sums_independent_replicas():
    beads = Beads(2, 4)
    beads._m[:] = [2.0, 3.0]
    beads._p[:] = np.arange(24).reshape(4, 6) / 10
    expected = np.sum(beads.p**2 / np.repeat([2.0, 3.0], 3)) / 2
    assert beads.kin == pytest.approx(expected)

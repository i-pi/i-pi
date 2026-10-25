"""Analytic energies and forces for isotropic, anisotropic, and batched oscillators."""

import numpy as np
import pytest

from ipi.pes.harmonic import Harmonic_driver


@pytest.mark.parametrize("constants", [(2.0,), (2.0, 3.0, 5.0)])
@pytest.mark.parametrize("natoms", [1, 3])
def test_harmonic_forces_match_energy_gradient(constants, natoms):
    """The force is the negative energy derivative along each Cartesian component."""
    driver = Harmonic_driver(*constants)
    positions = np.random.default_rng(20).normal(size=(natoms, 3))
    cell = np.eye(3) * 10
    energy, force, virial, _ = driver.compute(cell, positions)
    expected_energy = 0.5 * np.sum(np.asarray(constants) * positions**2)
    assert energy == pytest.approx(expected_energy)
    delta = 1e-5
    for index in np.ndindex(positions.shape):
        forward, backward = positions.copy(), positions.copy()
        forward[index] += delta
        backward[index] -= delta
        derivative = (
            driver.compute(cell, forward)[0] - driver.compute(cell, backward)[0]
        ) / (2 * delta)
        assert force[index] == pytest.approx(-derivative, rel=1e-8, abs=1e-8)
    np.testing.assert_array_equal(virial, np.zeros((3, 3)))


@pytest.mark.parametrize("constants", [(2.0,), (2.0, 3.0, 5.0)])
def test_harmonic_batch_matches_independent_structures(constants):
    """Batching preserves the energy and force for every structure in order."""
    driver = Harmonic_driver(*constants)
    cells = [np.eye(3) * 10 for _ in range(3)]
    positions = [np.random.default_rng(seed).normal(size=(2, 3)) for seed in range(3)]
    batch = driver.compute(cells, positions)
    assert len(batch) == 3
    for result, cell, position in zip(batch, cells, positions):
        expected = driver.compute(cell, position)
        for actual, reference in zip(result[:3], expected[:3]):
            np.testing.assert_allclose(actual, reference, rtol=1e-14, atol=1e-14)
        assert result[3] == expected[3]


def test_harmonic_batch_rejects_missing_structures():
    """A mismatched batch must fail instead of silently dropping an input."""
    with pytest.raises(ValueError, match="lists"):
        Harmonic_driver(1.0).compute([np.eye(3), np.eye(3)], [np.zeros((1, 3))])

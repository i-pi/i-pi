"""Unit tests for the ideal-gas Python driver."""

import numpy as np
import pytest

from ipi.pes.gas import Gas_driver


def test_gas_driver_single_structure_returns_zero_observables():
    """Checks that the single-structure gas driver is a no-op."""

    driver = Gas_driver()
    cell = np.eye(3)
    pos = np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]])

    pot, force, virial, extra = driver.compute(cell, pos)

    assert pot == 0.0
    assert np.allclose(force, np.zeros_like(pos))
    assert np.allclose(virial, np.zeros((3, 3)))
    assert extra == ""


def test_gas_driver_batched_structures_return_zero_observables():
    """Checks that batched gas-driver evaluations return one zero result each."""

    driver = Gas_driver()
    cells = [np.eye(3), np.eye(3) * 2.0]
    positions = [
        np.array([[0.1, 0.2, 0.3]]),
        np.array([[0.4, 0.5, 0.6]]),
    ]

    results = driver.compute(cells, positions)

    assert len(results) == len(positions)
    for pos, (pot, force, virial, extra) in zip(positions, results):
        assert pot == 0.0
        assert np.allclose(force, np.zeros_like(pos))
        assert np.allclose(virial, np.zeros((3, 3)))
        assert extra == ""


def test_gas_driver_batched_requires_matching_lists():
    """Checks that batched gas-driver inputs must be parallel lists."""

    driver = Gas_driver()

    with pytest.raises(ValueError, match="Both position and cell"):
        driver.compute([np.eye(3)], np.array([[0.1, 0.2, 0.3]]))

"""Unit tests for the double-well Python driver."""

import numpy as np

from ipi.pes.doublewell import DoubleWell_driver


def expected_observables(driver, positions):
    """Returns the analytic energy and forces for a set of positions."""

    energy = (
        driver.A * (positions[:, 0] - driver.delta) ** 2
        + driver.B * positions[:, 0] ** 4
        + 0.5 * driver.k * positions[:, 1] ** 2
        + 0.5 * driver.k * positions[:, 2] ** 2
    ).sum()
    forces = np.empty_like(positions)
    forces[:, 0] = (
        -2.0 * driver.A * (positions[:, 0] - driver.delta)
        - 4.0 * driver.B * positions[:, 0] ** 3
    )
    forces[:, 1:] = -driver.k * positions[:, 1:]
    return energy, forces


def test_doublewell_driver_returns_scalar_energy_for_single_and_batched_inputs():
    """Checks multi-atom and batched calls return one scalar energy each."""

    driver = DoubleWell_driver(w_b=500.0, v0=300.0, m=1837.36223469, delta=0.0)
    cell = np.eye(3)
    positions = [
        np.array([[0.25, 0.10, -0.20], [-0.40, 0.30, 0.15]]),
        np.array([[0.55, -0.10, 0.05]]),
    ]

    for pos, (energy, forces, virial, extra) in zip(
        positions, driver.compute([cell, 2.0 * cell], positions)
    ):
        expected_energy, expected_forces = expected_observables(driver, pos)
        assert np.isscalar(energy)
        assert np.isclose(energy, expected_energy)
        assert np.allclose(forces, expected_forces)
        assert np.allclose(virial, np.zeros((3, 3)))
        assert extra == "empty"

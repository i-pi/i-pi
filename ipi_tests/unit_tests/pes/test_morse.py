"""Tests for the Morse-harmonic Python driver."""

import numpy as np

from ipi.pes.morse import MorseHarmonic_driver
from ipi.utils.units import unit_to_internal


def test_morse_driver_converts_user_parameters_to_internal_units():
    """Checks conversion of the four user-facing parameters."""

    driver = MorseHarmonic_driver(De=2.0, a=3.0, z0=1.5, k=4.0)
    length = unit_to_internal("length", "angstrom", 1.0)
    energy = unit_to_internal("energy", "electronvolt", 1.0)

    assert driver.De == 2.0 * energy
    assert driver.a == 3.0 / length
    assert driver.z0 == 1.5 * length
    assert driver.k == 4.0 * energy / length**2


def test_morse_driver_force_matches_total_energy_gradient():
    """Checks returned forces against small changes in total energy."""

    driver = MorseHarmonic_driver(De=2.0, a=3.0, z0=1.5, k=4.0)
    positions = np.array([[0.2, -0.3, driver.z0 + 0.1], [-0.4, 0.5, driver.z0 - 0.2]])
    cell = np.eye(3)
    _, forces, virial, extra = driver.compute_structure(cell, positions)
    step = 1.0e-6
    numerical_forces = np.zeros_like(positions)

    for index in np.ndindex(positions.shape):
        forward = positions.copy()
        backward = positions.copy()
        forward[index] += step
        backward[index] -= step
        forward_energy = driver.compute_structure(cell, forward)[0]
        backward_energy = driver.compute_structure(cell, backward)[0]
        numerical_forces[index] = -(forward_energy - backward_energy) / (2.0 * step)

    assert np.allclose(forces, numerical_forces, rtol=1.0e-6, atol=1.0e-8)
    assert np.array_equal(virial, np.zeros((3, 3)))
    assert extra == "empty"

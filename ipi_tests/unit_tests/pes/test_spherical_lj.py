"""Tests for the spherical Lennard-Jones Python driver."""

import numpy as np

from ipi.pes.spherical_LJ import Spherical_LJ_driver, compute_energy_and_forces
from ipi.utils.units import unit_to_internal


def _instructions():
    return {
        "center": [0.0, 0.0, 0.0],
        "radius": 5.0,
        "sigma": 1.2,
        "epsilon": 0.7,
        "symbols": ["O"],
        "first_power": 9,
        "second_power": 3,
    }


def test_spherical_lj_forces_match_energy_gradient():
    """Checks analytical forces against small changes in energy."""

    instructions = _instructions()
    positions = np.array([[1.0, 0.2, -0.1], [-0.4, 1.3, 0.2]])
    _, forces = compute_energy_and_forces(positions, instructions)
    step = 1.0e-6
    numerical_forces = np.zeros_like(positions)

    for index in np.ndindex(positions.shape):
        forward = positions.copy()
        backward = positions.copy()
        forward[index] += step
        backward[index] -= step
        forward_energy = compute_energy_and_forces(forward, instructions)[0]
        backward_energy = compute_energy_and_forces(backward, instructions)[0]
        numerical_forces[index] = -(forward_energy - backward_energy) / (2.0 * step)

    assert np.allclose(forces, numerical_forces, rtol=1.0e-6, atol=1.0e-8)


def test_spherical_lj_driver_applies_only_to_selected_species():
    """Checks that unselected atoms receive no potential force."""

    instructions = _instructions()
    driver = Spherical_LJ_driver(
        template="", instructions=instructions, symbols=["H", "O"]
    )
    positions = np.array([[1.0, 0.0, 0.0], [0.0, 1.5, 0.0]])

    energy, forces, virial, _ = driver.compute_structure(np.eye(3), positions)
    expected_energy, expected_force = compute_energy_and_forces(
        positions[1:], instructions
    )

    assert energy == expected_energy
    assert np.array_equal(forces[0], np.zeros(3))
    assert np.allclose(forces[1], expected_force[0])
    assert np.array_equal(virial, np.zeros((3, 3)))


def test_spherical_lj_driver_converts_instruction_units():
    """Checks conversion of length and energy instructions."""

    instructions = _instructions()
    instructions.update(
        {
            "center_unit": "angstrom",
            "radius_unit": "angstrom",
            "sigma_unit": "angstrom",
            "epsilon_unit": "electronvolt",
        }
    )
    original = {
        key: np.asarray(instructions[key]).copy()
        for key in ("center", "radius", "sigma", "epsilon")
    }

    driver = Spherical_LJ_driver(template="", instructions=instructions, symbols=["O"])
    length = unit_to_internal("length", "angstrom", 1.0)
    energy = unit_to_internal("energy", "electronvolt", 1.0)

    assert np.allclose(driver.instructions["center"], original["center"] * length)
    assert driver.instructions["radius"] == original["radius"] * length
    assert driver.instructions["sigma"] == original["sigma"] * length
    assert driver.instructions["epsilon"] == original["epsilon"] * energy
    assert not any(key.endswith("_unit") for key in driver.instructions)

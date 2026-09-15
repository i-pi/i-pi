"""Deals with testing the Atoms object."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np

from ..common.folder import local

from ipi.engine.atoms import Atoms
from ipi.utils.io import read_file


def get_atoms(fin):
    """Reads atoms object from file @fin."""

    with open(local(fin), "r") as f:
        ret = read_file("xyz", f)
    return ret["atoms"]


def test_names():
    """Tests names of Atoms object."""
    atoms = get_atoms("test.pos_0.xyz")
    expected = ["O", "H", "H"]
    assert len(atoms.names) == 3
    for i, name in enumerate(atoms.names):
        print(atoms[i])
        print(atoms[i].q)
        assert name == expected[i]
        assert name == atoms[i].name

    # Same test with iterator instead
    for i, atom in enumerate(atoms):
        assert atom.name == expected[i]
        assert atom.name == atoms.names[i]


def test_kinetic_energy_and_stress_match_manual_calculation():
    """Checks kinetic observables against their defining expressions."""

    atoms = Atoms(2)
    masses = np.array([2.0, 4.0])
    momenta = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    atoms.m[:] = masses
    atoms.p[:] = momenta.flatten()

    expected_kinetic = np.sum(momenta**2, axis=1) / (2.0 * masses)
    expected_atom_stress = np.array(
        [
            np.triu(np.outer(momentum, momentum) / mass)
            for momentum, mass in zip(momenta, masses)
        ]
    )

    assert np.allclose([atom.kin.item() for atom in atoms], expected_kinetic)
    assert np.allclose([atom.kstress for atom in atoms], expected_atom_stress)
    assert np.isclose(atoms.kin, np.sum(expected_kinetic))
    assert np.allclose(atoms.kstress, np.sum(expected_atom_stress, axis=0))


def test_clone_copies_values_without_sharing_storage():
    """Checks that changing a clone does not change the original atoms."""

    atoms = Atoms(2)
    positions = np.arange(6, dtype=float)
    momenta = positions + 10.0
    masses = np.array([1.0, 2.0])
    names = np.array(["H", "O"])
    atoms.q[:] = positions
    atoms.p[:] = momenta
    atoms.m[:] = masses
    atoms.names[:] = names

    cloned = atoms.clone()

    assert np.array_equal(cloned.q, positions)
    assert np.array_equal(cloned.p, momenta)
    assert np.array_equal(cloned.m, masses)
    assert np.array_equal(cloned.names, names)

    cloned.q[0] = -1.0
    cloned.p[1] = -2.0
    cloned.m[0] = 3.0
    cloned.names[1] = "C"

    assert np.array_equal(atoms.q, positions)
    assert np.array_equal(atoms.p, momenta)
    assert np.array_equal(atoms.m, masses)
    assert np.array_equal(atoms.names, names)

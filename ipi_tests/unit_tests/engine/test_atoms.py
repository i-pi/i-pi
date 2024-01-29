"""Deals with testing the Atoms object."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


from ..common.folder import local

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

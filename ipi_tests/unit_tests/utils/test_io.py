"""Deals with testing the io system.

Note that this will only run if you have Python version 2.5 or later.
Otherwise, replace all the with statements with f = filestream.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import filecmp
import os
import re

import numpy as np
from numpy.testing import assert_equal

from ..common.folder import local
from ipi.utils.units import unit_to_internal
from ipi.utils.io import iter_file, read_file, print_file
from ipi.utils.parsing import read_trajectory


pos_xyz = np.array([i for i in range(3 * 3)])
pos_pdb = unit_to_internal("length", "angstrom", pos_xyz)

# Regular expression to match the cell parameters
abcABC = re.compile(r"CELL[\(\[\{]abcABC[\)\]\}]: ([-+0-9\.Ee ]*)\s*")


def read_cell_parameters_from_extxyz(filename):
    """
    Reads the cell parameters from an external XYZ file using a regular expression.
    """
    # Open the file and read the relevant line
    with open(filename, "r") as file:
        file.readline()  # Skip the first line (number of atoms)
        cell_line = file.readline()  # Read the second line (cell parameters)

    # Use the regular expression to search for the cell parameters
    _match = abcABC.search(cell_line)
    if _match:
        # Extract the cell parameters and convert them to a numpy array of floats
        cell_params = np.array(_match.group(1).split(), dtype=float)
        return cell_params
    else:
        raise ValueError("Cell parameters not found in the given file format.")


def read_positions_from_extxyz(filename):
    """
    Reads atomic positions from an external XYZ file, skipping the header.
    """
    return np.genfromtxt(
        filename, skip_header=2, dtype=None, encoding=None, usecols=(1, 2, 3)
    )


def test_read_pos_au():
    """
    Tests reading atomic positions in atomic units (bohr) and converts them to angstroms.
    """
    # positions saved in `atomic_unit`=`bohr`
    file = local("test.positions_au.xyz")

    # positions converted within `read_trajectory` into `angstrom`
    atoms = read_trajectory(file)
    assert len(atoms) == 1, "Wrong number of snapshots"
    atoms = atoms[0]

    # conversion factor from `bohr` to `angstrom`
    ang2bohr = unit_to_internal("length", "angstrom", 1)

    # extract the positions (in bohr) with brute force
    positions = np.genfromtxt(
        file, skip_header=2, dtype=None, encoding=None, usecols=(1, 2, 3)
    )

    # check that `read_trajectory` actually converts to  `angstrom`
    assert np.allclose(atoms.positions * ang2bohr, positions), "Wrong conversion"


def test_read_cell_au():
    """
    Tests reading cell parameters in atomic units (bohr) and converts them to angstroms.
    """
    # positions saved in `atomic_unit`=`bohr`
    file = local("test.positions_au.xyz")

    # positions converted within `read_trajectory` into `angstrom`
    atoms = read_trajectory(file)
    assert len(atoms) == 1, "Wrong number of snapshots"
    atoms = atoms[0]

    # conversion factor from `bohr` to `angstrom`
    ang2bohr = unit_to_internal("length", "angstrom", 1)
    atoms_cell = atoms.get_cell().cellpar()
    atoms_cell[:3] *= ang2bohr

    # extract the positions (in bohr) with brute force
    cell = read_cell_parameters_from_extxyz(file)
    print(cell)

    # check that `read_trajectory` actually converts to  `angstrom`
    assert np.allclose(atoms_cell, cell), "Wrong conversion"


def test_read_pos_ang():
    """
    Tests reading atomic positions in angstroms without conversion.
    """
    # positions saved in `angstrom`
    file = local("test.positions_ang.xyz")

    # positions are not converted
    atoms = read_trajectory(file)
    assert len(atoms) == 1, "Wrong number of snapshots"
    atoms = atoms[0]

    # extract the positions (in angstrom) with brute force
    positions = np.genfromtxt(
        file, skip_header=2, dtype=None, encoding=None, usecols=(1, 2, 3)
    )

    # check that `read_trajectory` actually does not convert the positions
    assert np.allclose(atoms.positions, positions), "Wrong conversion"


def test_read_cell_ang():
    """
    Tests reading cell parameters in angstroms without conversion.
    """
    # positions saved in `atomic_unit`=`bohr`
    file = local("test.positions_ang.xyz")

    # positions converted within `read_trajectory` into `angstrom`
    atoms = read_trajectory(file)
    assert len(atoms) == 1, "Wrong number of snapshots"
    atoms = atoms[0]
    atoms_cell = atoms.get_cell().cellpar()

    # extract the positions (in bohr) with brute force
    cell = read_cell_parameters_from_extxyz(file)

    # check that `read_trajectory` actually converts to  `angstrom`
    assert np.allclose(atoms_cell, cell), "Wrong conversion"


def test_read_xyz():
    """Tests that xyz files are read correctly."""

    with open(local("test.pos_0.xyz"), "r") as f:
        ret = read_file("xyz", f)
        atoms = ret["atoms"]
        assert len(atoms) == 3
        assert_equal(pos_xyz, atoms.q)


def test_iter_xyz():
    """Tests that xyz files with multiple frames are read correctly."""

    with open(local("test.pos_0.xyz"), "r") as f:
        for num, ret in enumerate(iter_file("xyz", f)):
            atoms = ret["atoms"]
            assert len(atoms) == 3
            assert_equal(pos_xyz * (num + 1), atoms.q)


def test_read_pdb():
    """Tests that pdb files are read correctly."""

    with open(local("test.pos_0.pdb"), "r") as f:
        ret = read_file("pdb", f)
        atoms = ret["atoms"]
        assert len(atoms) == 3
        assert_equal(pos_pdb, atoms.q)
        # TODO: test cell


def test_iter_pdb():
    """Tests that pdb files with multiple frames are read correctly."""

    with open(local("test.pos_0.pdb"), "r") as f:
        for num, ret in enumerate(iter_file("pdb", f)):
            atoms = ret["atoms"]
            assert len(atoms) == 3
            assert_equal(pos_pdb * (num + 1), atoms.q)


def test_print_pdb():
    """Tests that pdb files are printed correctly."""

    with open(local("test.pos_0.pdb"), "r") as f:
        with open(local("test.pos_1.pdb"), "w") as out:
            for num, ret in enumerate(iter_file("pdb", f)):
                atoms = ret["atoms"]
                assert len(atoms) == 3
                assert_equal(pos_pdb * (num + 1), atoms.q)
                print_file("pdb", atoms, ret["cell"], filedesc=out)

    assert filecmp.cmp(local("test.pos_0.pdb"), local("test.pos_1.pdb"))
    os.unlink(local("test.pos_1.pdb"))


def test_print_xyz():
    """Tests that xyz files are printed correctly."""

    with open(local("test.pos_0.xyz"), "r") as f:
        with open(local("test.pos_1.xyz"), "w") as out:
            for num, ret in enumerate(iter_file("xyz", f)):
                atoms = ret["atoms"]
                assert len(atoms) == 3
                assert_equal(pos_xyz * (num + 1), atoms.q)
                print_file("xyz", atoms, ret["cell"], filedesc=out)

    assert filecmp.cmp(local("test.pos_0.xyz"), local("test.pos_1.xyz"))
    os.unlink(local("test.pos_1.xyz"))


def test_read_xyz2():
    """Tests that mode/xyz files are read correctly."""

    with open(local("test.pos_0.xyz"), "r") as f:
        ret = read_file("xyz", f)
        atoms = ret["atoms"]
        assert len(atoms) == 3
        assert_equal(pos_xyz, atoms.q)


def test_iter_xyz2():
    """Tests that mode/xyz files with multiple frames are read correctly."""

    with open(local("test.pos_0.xyz"), "r") as f:
        for num, ret in enumerate(iter_file("xyz", f)):
            atoms = ret["atoms"]
            assert len(atoms) == 3
            assert_equal(pos_xyz * (num + 1), atoms.q)


def test_read_pdb2():
    """Tests that mode/pdb files are read correctly."""

    with open(local("test.pos_0.pdb"), "r") as f:
        ret = read_file("pdb", f)
        atoms = ret["atoms"]
        assert len(atoms) == 3
        assert_equal(pos_pdb, atoms.q)
        # TODO: test cell


def test_iter_pdb2():
    """Tests that mode/pdb files with multiple frames are read correctly."""

    with open(local("test.pos_0.pdb"), "r") as f:
        for num, ret in enumerate(iter_file("pdb", f)):
            atoms = ret["atoms"]
            assert len(atoms) == 3
            assert_equal(pos_pdb * (num + 1), atoms.q)


def test_print_pdb2():
    """Tests that mode/pdb files are printed correctly."""

    with open(local("test.pos_0.pdb"), "r") as f:
        with open(local("test.pos_1.pdb"), "w") as out:
            for num, ret in enumerate(iter_file("pdb", f)):
                atoms = ret["atoms"]
                assert len(atoms) == 3
                assert_equal(pos_pdb * (num + 1), atoms.q)
                print_file("pdb", atoms, ret["cell"], filedesc=out)

    assert filecmp.cmp(local("test.pos_0.pdb"), local("test.pos_1.pdb"))
    os.unlink(local("test.pos_1.pdb"))


def test_print_xyz2():
    """Tests that mode/xyz files are printed correctly."""

    with open(local("test.pos_0.xyz"), "r") as f:
        with open(local("test.pos_1.xyz"), "w") as out:
            for num, ret in enumerate(iter_file("xyz", f)):
                atoms = ret["atoms"]
                assert len(atoms) == 3
                assert_equal(pos_xyz * (num + 1), atoms.q)
                print_file("xyz", atoms, ret["cell"], filedesc=out)

    assert filecmp.cmp(local("test.pos_0.xyz"), local("test.pos_1.xyz"))
    os.unlink(local("test.pos_1.xyz"))

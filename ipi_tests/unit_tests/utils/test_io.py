"""Deals with testing the io system.

Note that this will only run if you have Python version 2.5 or later.
Otherwise, replace all the with statements with f = filestream.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import filecmp
import os

import numpy as np
from numpy.testing import assert_equal

from ..common.folder import local
from ipi.utils.units import unit_to_internal
from ipi.utils.io import iter_file, read_file, print_file
from ipi.utils.parsing import read_output


pos_xyz = np.array([i for i in range(3 * 3)])
pos_pdb = unit_to_internal("length", "angstrom", pos_xyz)


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


def test_read_output_low_verbosity():
    """Test the read_output function with low verbosity output."""
    # run locally with: python -m pytest test_io.py::test_read_output_low_verbosity -v -s
    file = local("prop.low.out")
    a, b = read_output(file)

    assert a is not None, "Expected non-None result for 'a'"
    assert b is not None, "Expected non-None result for 'b'"
    assert a.keys() == b.keys(), "Keys should be the same"
    keys = list(a.keys())
    assert keys == [
        "step",
        "time",
        "conserved",
        "kinetic_md",
        "potential",
        "Efield",
        "Eenvelope",
        "dipole",
    ], "wrong keys"
    for k, v in a.items():
        assert len(v) == 3, f"`{k}` has the wrong size: {len(v)}"


def test_read_output_high_verbosity():
    """Test the read_output function with high verbosity output."""
    # run locally with: python -m pytest test_io.py::test_read_output_high_verbosity -v -s
    file = local("prop.high.out")
    a, b = read_output(file)

    assert a is not None, "Expected non-None result for 'a'"
    assert b is not None, "Expected non-None result for 'b'"
    assert a.keys() == b.keys(), "Keys should be the same"
    keys = list(a.keys())
    assert keys == [
        "step",
        "time",
        "conserved",
        "kinetic_md",
        "potential",
        "Efield",
        "Eenvelope",
        "dipole(0)",
        "dipole(1)",
    ], "wrong keys"
    for k, v in a.items():
        assert len(v) == 3, f"`{k}` has the wrong size: {len(v)}"


if __name__ == "__main__":
    test_read_output_low_verbosity()
    test_read_output_low_verbosity()

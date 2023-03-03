#!/usr/bin/env python3
# pylint: disable=C0111,W0621,R0914,C0301
# +easier to find important problems

import re
import tempfile as tmp
import filecmp
import os

import pytest
import numpy as np
import numpy.testing as npt

from ....common import xyz_generator as xyz_gen
import ipi.utils.io.backends.io_xyz as io_xyz
import ipi.utils.mathtools as mt

from ipi.engine.atoms import Atoms
from ipi.engine.cell import Cell
from ipi.utils.units import Elements

#
# Testing reading xyz #
#

deg2rad = np.pi / 180.0
cell_string = " ".join(
    [
        str(x)
        for x in mt.abc2h(
            5.1, 5.2, 5.0, 91 * deg2rad, 89 * deg2rad, 90 * deg2rad
        ).flatten()
    ]
)
# default_cell_mat = mt.abc2h(-1.0, -1.0, -1.0, np.pi/2.0, np.pi/2.0, np.pi/2.0) # After changing the input standard
default_cell_mat = np.eye(3) * -1


# natoms, frames, comment, expected_cell, precision
tests_read_xyz = [
    (1, 1, "just a string plus few numbers: 1.10 2 .1", default_cell_mat, 5),
    (2, 1, "another random comment", default_cell_mat, 5),
    (1, 2, "random comment", default_cell_mat, 5),
    (2, 2, "random comment", default_cell_mat, 5),
    (10, 2, "random comment", default_cell_mat, 5),
    (10, 10, "random comment", default_cell_mat, 5),
    (
        2,
        3,
        "# CELL(abcABC): 5.1 5.2 5.0 91.0  89  90",
        mt.abc2h(5.1, 5.2, 5.0, 91 * deg2rad, 89 * deg2rad, 90 * deg2rad),
        5,
    ),
    (
        2,
        3,
        "# CELL[GENH]: " + cell_string,
        mt.abc2h(5.1, 5.2, 5.0, 91 * deg2rad, 89 * deg2rad, 90 * deg2rad),
        5,
    ),
    (
        2,
        3,
        "# CELL{H}: " + cell_string,
        mt.abc2h(
            *mt.genh2abc(
                mt.abc2h(5.1, 5.2, 5.0, 91 * deg2rad, 89 * deg2rad, 90 * deg2rad)
            )
        ),
        5,
    ),
    (
        2,
        3,
        "100 aaa # CELL(abcABC): 5.1 5.2 5.0 91.0  89  90 100 aaa",
        mt.abc2h(5.1, 5.2, 5.0, 91 * deg2rad, 89 * deg2rad, 90 * deg2rad),
        5,
    ),
    (
        2,
        3,
        "100 aaa # CELL[GENH]: " + cell_string + " 100 aaa",
        mt.abc2h(5.1, 5.2, 5.0, 91 * deg2rad, 89 * deg2rad, 90 * deg2rad),
        5,
    ),
    (
        2,
        3,
        "# CELL{H}: " + cell_string + " 100 aaa",
        mt.abc2h(
            *mt.genh2abc(
                mt.abc2h(5.1, 5.2, 5.0, 91 * deg2rad, 89 * deg2rad, 90 * deg2rad)
            )
        ),
        5,
    ),
]


@pytest.fixture(params=tests_read_xyz)
def create_random_xyz_traj_to_read(request):
    natoms, frames, comment, expected_cell, precision = request.param

    filedesc, xyz, atom_names = xyz_gen.xyz_traj_filedesc(natoms, frames, comment)
    filedesc.seek(0)

    check_comment = re.match(
        r"# CELL[\(\[\{]H[\)\]\}]: ([-0-9\.?Ee ]*)\s*", comment
    )  # pylint: disable=W1401

    if check_comment:
        genh = np.array(check_comment.group(1).split()[:9], float).reshape((3, 3))
        invgenh = np.linalg.inv(genh)
        for _ui in range(natoms * frames):
            _uu = np.array([xyz[3 * _ui], xyz[3 * _ui + 1], xyz[3 * _ui + 2]])
            _us = np.dot(_uu, invgenh)
            _uu = np.dot(expected_cell, _us)
            xyz[3 * _ui], xyz[3 * _ui + 1], xyz[3 * _ui + 2] = _uu

    return (
        filedesc,
        xyz,
        atom_names,
        natoms,
        frames,
        comment,
        expected_cell,
        precision,
    )


@pytest.mark.skip(reason="This needs to be updated to match current code.")
def test_read_xyz(create_random_xyz_traj_to_read):
    (
        filedesc,
        xyz,
        atom_names,
        natoms,
        frames,
        comment,
        expected_cell,
        precision,
    ) = create_random_xyz_traj_to_read

    for _fr in range(frames):
        tcomment, tcell, tqatoms, tnames, tmasses = io_xyz.read_xyz(filedesc)

        assert tcomment.strip() == comment
        npt.assert_array_almost_equal(
            tqatoms,
            xyz[_fr * natoms * 3 : _fr * 3 * natoms + 3 * natoms],
            decimal=precision,
        )
        npt.assert_array_equal(
            np.array(atom_names[_fr * natoms : _fr * natoms + natoms], dtype="|S4"),
            tnames,
        )
        npt.assert_array_almost_equal(tcell, expected_cell, decimal=precision)


@pytest.mark.skipif(True, reason="iter_xyz has been removed")
def test_iter_xyz(create_random_xyz_traj_to_read):
    (
        filedesc,
        xyz,
        atom_names,
        natoms,
        junk,
        comment,
        expected_cell,
        precision,
    ) = create_random_xyz_traj_to_read

    _fr = 0
    for _io in io_xyz.iter_xyz(filedesc):
        tcomment, tcell, tqatoms, tnames, tmasses = _io
        assert tcomment.strip() == comment
        npt.assert_array_almost_equal(
            tqatoms,
            xyz[_fr * natoms * 3 : _fr * 3 * natoms + 3 * natoms],
            decimal=precision,
        )
        npt.assert_array_equal(
            np.array(atom_names[_fr * natoms : _fr * natoms + natoms], dtype="|S4"),
            tnames,
        )
        npt.assert_array_almost_equal(tcell, expected_cell, decimal=precision)
        _fr += 1


#
# Testing writing xyz #
#

write_test_xyz = [
    (1, 1, "just a string plus few numbers: 1.10 2 .1", default_cell_mat, 10),
    (2, 1, "another random comment", default_cell_mat, 10),
    (1, 2, "random comment", default_cell_mat, 10),
    (2, 2, "random comment", default_cell_mat, 10),
    (10, 2, "random comment", default_cell_mat, 10),
    (10, 10, "random comment", default_cell_mat, 10),
    (
        2,
        3,
        "# CELL(abcABC): 5.1 5.2 5.0 91.0  89  90",
        mt.abc2h(5.1, 5.2, 5.0, 91 * deg2rad, 89 * deg2rad, 90 * deg2rad),
        10,
    ),
    (
        2,
        3,
        "# CELL[GENH]: " + cell_string,
        mt.abc2h(5.1, 5.2, 5.0, 91 * deg2rad, 89 * deg2rad, 90 * deg2rad),
        10,
    ),
    (
        2,
        3,
        "# CELL{H}: " + cell_string,
        mt.abc2h(
            *mt.genh2abc(
                mt.abc2h(5.1, 5.2, 5.0, 91 * deg2rad, 89 * deg2rad, 90 * deg2rad)
            )
        ),
        10,
    ),
]


@pytest.fixture(params=write_test_xyz)
def create_random_xyz_traj_to_write(request):
    natoms, frames, comment, expected_cell, precision = request.param

    a, b, c, alpha, beta, gamma = mt.h2abc_deg(expected_cell)

    fmt_header = "# CELL(abcABC): %10.5f  %10.5f  %10.5f  %10.5f  %10.5f  %10.5f  %s"

    comment = fmt_header % (a, b, c, alpha, beta, gamma, comment)

    filedesc, xyz, atom_names = xyz_gen.xyz_traj_filedesc(natoms, frames, comment)
    filedesc.seek(0)

    masses = [Elements.mass(_am) for _am in atom_names]

    cell_list = []
    atoms_list = []

    for _fr in range(frames):
        cell = Cell(expected_cell)
        atoms = Atoms(natoms)
        atoms.q[:] = xyz[_fr * natoms * 3 : (_fr + 1) * natoms * 3]
        atoms.names = atom_names[_fr * natoms : (_fr + 1) * natoms]
        atoms.m[:] = masses[_fr * natoms : (_fr + 1) * natoms]
        atoms_list.append(atoms)
        cell_list.append(cell)

    return (filedesc, atoms_list, cell_list, comment, precision)


def test_print_xyz(create_random_xyz_traj_to_write):
    filedesc, atoms_list, cell_list, title, precision = create_random_xyz_traj_to_write

    filedesc_orig = tmp.NamedTemporaryFile(
        mode="w", prefix="ipi_testing-tmp", delete=False
    )
    filedesc_test = tmp.NamedTemporaryFile(
        mode="w", prefix="ipi_testing-tmp", delete=False
    )

    filedesc_orig.write(filedesc.read())
    filedesc.close()

    filedesc_orig.flush()
    filedesc_test.flush()

    for atoms, cell in zip(atoms_list, cell_list):
        io_xyz.print_xyz(atoms, cell, filedesc=filedesc_test, title=title[88:])

    filedesc_orig.close()
    filedesc_test.close()

    assert filecmp.cmp(filedesc_orig.name, filedesc_test.name)

    os.remove(filedesc_orig.name)
    os.remove(filedesc_test.name)


# def test_print_xyz(atoms, cell, filedesc=sys.stdout, title="")


# if __name__ == '__main__':

#     for st in test_cell:
#         string, float_n = st
#         floats = 100.0 * np.random.random_sample(float_n)

# Generate also a random number of spaces between numbers
#         spaces = ' ' * np.random.random_integers(1, 12)
#         string += spaces
#         string += spaces.join([str(x) for x in floats.tolist()])
#         print string

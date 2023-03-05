#!/usr/bin/env python3

import mock
import pytest
from pytest_mock import mocker

import numpy as np
import numpy.testing as npt

from ....common import xyz_generator as xyz_gen

import ipi.utils.io as io
import ipi.utils.mathtools as mt


# default_cell_mat = mt.abc2h(-1.0, -1.0, -1.0, np.pi/2.0, np.pi/2.0, np.pi/2.0) # After changing the "input standard"
default_cell_mat = np.eye(3) * -1.0
deg2rad = np.pi / 180.0


test_read_file_prms = [
    # natoms, nframes, comment, output_type, file_type, expected_cell, unit_conv_cell, unit_conv_q
    (1, 1, "asdasd", "objects", "xyz", default_cell_mat, 1, 1),
    (2, 1, "asdasd", "objects", "xyz", default_cell_mat, 1, 1),
    (1, 2, "asdasd", "objects", "xyz", default_cell_mat, 1, 1),
    (
        10,
        10,
        "cell{angstrom} positions{angstrom}",
        "objects",
        "xyz",
        default_cell_mat,
        1.8897261,
        1.8897261,
    ),
    (
        10,
        10,
        " positions{angstrom} cell{angstrom} 100 aaa # CELL(abcABC): 5.1 5.2 5.0 91.0  89  90 100 aaa",
        "objects",
        "xyz",
        mt.abc2h(5.1, 5.2, 5.0, 91 * deg2rad, 89 * deg2rad, 90 * deg2rad),
        1.8897261,
        1.8897261,
    ),
    (1, 1, "asdasd", "array", "xyz", default_cell_mat, 1, 1),
    (2, 1, "asdasd", "array", "xyz", default_cell_mat, 1, 1),
    (1, 2, "asdasd", "array", "xyz", default_cell_mat, 1, 1),
    (
        10,
        10,
        "cell{angstrom} positions{angstrom}",
        "array",
        "xyz",
        default_cell_mat,
        1.8897261,
        1.8897261,
    ),
    (
        10,
        10,
        " positions{angstrom} cell{angstrom} 100 aaa # CELL(abcABC): 5.1 5.2 5.0 91.0  89  90 100 aaa",
        "array",
        "xyz",
        mt.abc2h(5.1, 5.2, 5.0, 91 * deg2rad, 89 * deg2rad, 90 * deg2rad),
        1.8897261,
        1.8897261,
    ),
]


@pytest.fixture(params=test_read_file_prms)
def prepare_read_file(request):
    (
        natoms,
        nframes,
        comment,
        output_type,
        file_type,
        expected_cell,
        unit_conv_cell,
        unit_conv_q,
    ) = request.param

    filedesc, xyz, atom_names = xyz_gen.xyz_traj_filedesc(natoms, nframes, comment)

    expected_q = xyz.copy()
    expected_names = atom_names[:]

    # if comment.find('CELL') < 0:
    #     unit_conv_cell = 1.0

    return (
        file_type,
        filedesc,
        output_type,
        expected_q,
        expected_cell,
        expected_names,
        unit_conv_cell,
        unit_conv_q,
    )


@pytest.mark.skip(reason="This needs to be updated to match current code.")
def test_read_file(prepare_read_file):
    (
        file_type,
        filedesc,
        output_type,
        expected_q,
        expected_cell,
        expected_names,
        unit_conv_cell,
        unit_conv_q,
    ) = prepare_read_file

    returned_q = np.array([])
    returned_cell = np.array([])
    returned_names = []

    while True:
        try:
            output = io.read_file(file_type, filedesc, output=output_type)
        except EOFError:
            break

        if output_type.strip() != "objects":
            returned_q = np.concatenate([returned_q, output["data"]])
            returned_cell = output["cell"].flatten()

        elif output_type.strip() == "objects":
            returned_q = np.concatenate([returned_q, output["atoms"].q])
            returned_cell = output["cell"].h.flatten()
            returned_names += output["atoms"].names.tolist()

    npt.assert_almost_equal(expected_q * unit_conv_q, returned_q, 5)
    npt.assert_almost_equal(expected_cell.flatten() * unit_conv_cell, returned_cell, 5)
    if output_type.strip() == "objects":
        assert all(
            [
                expected_names[_ii] == returned_names[_ii]
                for _ii in range(len(expected_names))
            ]
        )

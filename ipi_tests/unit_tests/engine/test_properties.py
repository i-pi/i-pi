import mock
import tempfile
import re

import pytest

import numpy as np
import numpy.testing as npt

import ipi.engine.properties
from ..common import xyz_generator as xyz_gen

test_Trajectories_print_traj_prms = [
    # natoms, nbeads, bead, cell, cell_units, property_, property_units, format_, units_conversion
    (1, 1, 0, np.random.rand(9), "atomic_unit", "positions", "atomic_unit", "xyz", 1),
    (
        1,
        1,
        0,
        np.random.rand(9),
        "angstrom",
        "positions",
        "atomic_unit",
        "xyz",
        0.52917721,
    ),
    (1, 1, 0, np.random.rand(9), "atomic_unit", "forces", "atomic_unit", "xyz", 1),
    (
        1,
        1,
        0,
        np.random.rand(9),
        "angstrom",
        "forces",
        "atomic_unit",
        "xyz",
        0.52917721,
    ),
    (1, 1, 0, np.random.rand(9), "atomic_unit", "positions", "atomic_unit", "pdb", 1),
    (
        1,
        1,
        0,
        np.random.rand(9),
        "angstrom",
        "positions",
        "atomic_unit",
        "pdb",
        0.52917721,
    ),
    (1, 1, 0, np.random.rand(9), "atomic_unit", "forces", "atomic_unit", "pdb", 1),
    (
        1,
        1,
        0,
        np.random.rand(9),
        "angstrom",
        "forces",
        "atomic_unit",
        "pdb",
        0.52917721,
    ),
    (5, 1, 0, np.random.rand(9), "atomic_unit", "positions", "atomic_unit", "xyz", 1),
    (
        5,
        1,
        0,
        np.random.rand(9),
        "angstrom",
        "positions",
        "atomic_unit",
        "xyz",
        0.52917721,
    ),
    (5, 1, 0, np.random.rand(9), "atomic_unit", "forces", "atomic_unit", "xyz", 1),
    (
        5,
        1,
        0,
        np.random.rand(9),
        "angstrom",
        "forces",
        "atomic_unit",
        "xyz",
        0.52917721,
    ),
    (5, 1, 0, np.random.rand(9), "atomic_unit", "positions", "atomic_unit", "pdb", 1),
    (
        5,
        1,
        0,
        np.random.rand(9),
        "angstrom",
        "positions",
        "atomic_unit",
        "pdb",
        0.52917721,
    ),
    (5, 1, 0, np.random.rand(9), "atomic_unit", "forces", "atomic_unit", "pdb", 1),
    (
        5,
        1,
        0,
        np.random.rand(9),
        "angstrom",
        "forces",
        "atomic_unit",
        "pdb",
        0.52917721,
    ),
    (5, 10, 1, np.random.rand(9), "atomic_unit", "positions", "atomic_unit", "xyz", 1),
    (
        5,
        10,
        2,
        np.random.rand(9),
        "angstrom",
        "positions",
        "atomic_unit",
        "xyz",
        0.52917721,
    ),
    (5, 10, 3, np.random.rand(9), "atomic_unit", "forces", "atomic_unit", "xyz", 1),
    (
        5,
        10,
        4,
        np.random.rand(9),
        "angstrom",
        "forces",
        "atomic_unit",
        "xyz",
        0.52917721,
    ),
    (5, 10, 5, np.random.rand(9), "atomic_unit", "positions", "atomic_unit", "pdb", 1),
    (
        5,
        10,
        6,
        np.random.rand(9),
        "angstrom",
        "positions",
        "atomic_unit",
        "pdb",
        0.52917721,
    ),
    (5, 10, 7, np.random.rand(9), "atomic_unit", "forces", "atomic_unit", "pdb", 1),
    (
        5,
        10,
        8,
        np.random.rand(9),
        "angstrom",
        "forces",
        "atomic_unit",
        "pdb",
        0.52917721,
    ),
]


def create_a_fake_system_obj(natoms, q, forces, atom_names, cell):
    sys_attr = {
        "beads.natoms": natoms,
        "beads.q": q,
        "simul.step": 0,
        "beads.names": atom_names[:natoms],
        "cell.h": cell,
        "forces.f": forces,
    }

    sys = mock.Mock(**sys_attr)

    return sys


@pytest.fixture(params=test_Trajectories_print_traj_prms)
def prepare_Trajectories_print_traj(request):
    (
        natoms,
        nbeads,
        bead,
        cell,
        cell_units,
        property_,
        property_unit,
        format_,
        unit_conv,
    ) = request.param

    junk, xyz, atom_names = xyz_gen.xyz_traj(natoms, nbeads, "comment")

    # Using the same xyz_generator to generate random forces
    junk, forces, atom_names = xyz_gen.xyz_traj(natoms, nbeads, "comment")
    xyz = xyz.reshape((nbeads, natoms * 3))
    forces = forces.reshape((nbeads, natoms * 3))
    cell = cell.reshape((3, 3))

    system_mock = create_a_fake_system_obj(
        natoms, xyz, forces, atom_names[:natoms], cell
    )

    stream = tempfile.NamedTemporaryFile(mode="w")

    if property_ == "forces":
        expected_position = forces.copy()
    elif property_ == "positions":
        expected_position = xyz.copy()
    expected_cell = cell.copy()

    expected_comment = re.compile(
        r"^\s*\cell{%s\}\s+Traj\:\s*%s\{%s\}\s+Step\:\s+%i\s*Bead\:\s*%i\s*$"
        % (cell_units, property_, property_unit, 1, bead)
    )

    expected_names = atom_names[:]

    return (
        system_mock,
        stream,
        bead,
        expected_position,
        expected_cell,
        expected_comment,
        expected_names,
        format_,
        property_,
        cell_units,
        unit_conv,
    )


@pytest.mark.skip(reason="This needs to be updated to match current code.")
def test_Trajectories_print_traj(prepare_Trajectories_print_traj, mocker):
    mock_io = mocker.patch("ipi.engine.properties.io.print_file", autospec=True)

    (
        system,
        stream,
        bead,
        expected_position,
        expected_cell,
        expected_comment,
        expected_names,
        format_,
        property_,
        cell_units,
        unit_conv,
    ) = prepare_Trajectories_print_traj

    ipi.engine.properties.io.print_file = mock_io

    trj = ipi.engine.properties.Trajectories()
    trj.bind(system)

    # Function to test call
    trj.print_traj(
        property_, stream, b=bead, format=format_, cell_units=cell_units, flush=True
    )

    comment = mock_io.call_args[1]["title"]
    file_type, atoms, cell, junk = mock_io.call_args[0]

    assert file_type == format_
    assert expected_comment.match(comment)
    npt.assert_almost_equal(atoms.q, expected_position[bead], 5)
    npt.assert_equal(atoms.names, expected_names[: system.beads.natoms])
    npt.assert_almost_equal(cell.h, expected_cell * unit_conv)

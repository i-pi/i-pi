"""Configuration, trajectory, and property-file input/output checks."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

import filecmp
import io

import numpy as np
import pytest
from numpy.testing import assert_equal

from ..common.folder import local
from ipi.engine.beads import Beads
from ipi.engine.cell import Cell
from ipi.utils.units import unit_to_internal, unit_to_user
from ipi.utils.io import (
    iter_file,
    read_file,
    print_file,
    print_file_path,
    netstring_encoded_savez,
    netstring_encoded_loadz,
)
from ipi.scripting import read_output

pos_xyz = np.arange(9)
pos_pdb = unit_to_internal("length", "angstrom", pos_xyz)

FORMATS = [
    pytest.param("xyz", pos_xyz, id="xyz"),
    pytest.param("pdb", pos_pdb, id="pdb"),
]


@pytest.mark.parametrize("mode,positions", FORMATS)
def test_read_configuration(mode, positions):
    with open(local(f"test.pos_0.{mode}")) as stream:
        atoms = read_file(mode, stream)["atoms"]
        assert len(atoms) == 3
        assert_equal(positions, atoms.q)


@pytest.mark.parametrize("mode,positions", FORMATS)
def test_iter_configurations(mode, positions):
    with open(local(f"test.pos_0.{mode}")) as stream:
        for index, frame in enumerate(iter_file(mode, stream)):
            assert len(frame["atoms"]) == 3
            assert_equal(positions * (index + 1), frame["atoms"].q)


@pytest.mark.parametrize("mode,positions", FORMATS)
def test_print_configurations(tmp_path, mode, positions):
    source = local(f"test.pos_0.{mode}")
    destination = tmp_path / f"test.pos_1.{mode}"
    with open(source) as stream, destination.open("w") as output:
        for index, frame in enumerate(iter_file(mode, stream)):
            atoms = frame["atoms"]
            assert len(atoms) == 3
            assert_equal(positions * (index + 1), atoms.q)
            print_file(mode, atoms, frame["cell"], filedesc=output)
    assert filecmp.cmp(source, destination)


@pytest.mark.parametrize(
    "verbosity,dipoles",
    [
        pytest.param("low", ["dipole"], id="low"),
        pytest.param("high", ["dipole(0)", "dipole(1)"], id="high"),
    ],
)
def test_read_output(verbosity, dipoles):
    data, metadata = read_output(local(f"prop.{verbosity}.out"))
    assert data is not None, "Expected property data"
    assert metadata is not None, "Expected property metadata"
    assert data.keys() == metadata.keys(), "Keys should be the same"
    assert (
        list(data)
        == [
            "step",
            "time",
            "conserved",
            "kinetic_md",
            "potential",
            "Efield",
            "Eenvelope",
        ]
        + dipoles
    ), "wrong keys"
    for key, values in data.items():
        assert len(values) == 3, f"{key} has the wrong size: {len(values)}"


@pytest.mark.parametrize("mode", ["xyz", "pdb"])
@pytest.mark.parametrize("units", ["atomic_unit", "angstrom"])
def test_path_output_roundtrip(mode, units):
    beads = Beads(2, 3)
    beads.names_[:] = ["H", "O"]
    beads.q_[:] = np.arange(18).reshape(3, 6) / 10
    cell = Cell(np.diag([10.0, 12.0, 14.0]))
    stream = io.StringIO()
    print_file_path(
        mode, beads, cell, stream, key="positions", units=units, cell_units=units
    )
    stream.seek(0)
    if mode == "pdb":
        lines = stream.getvalue().splitlines()
        atoms = [line for line in lines if line.startswith("ATOM")]
        positions = np.array(
            [
                [float(line[a:b]) for a, b in [(30, 38), (38, 46), (46, 54)]]
                for line in atoms
            ]
        )
        np.testing.assert_allclose(
            positions,
            beads.q.reshape(-1, 3) * unit_to_user("length", units, 1),
            atol=0.00051,
        )
        connections = {
            tuple(map(int, line.split()[1:]))
            for line in lines
            if line.startswith("CONECT")
        }
        assert connections == {(1, 5), (2, 6), (1, 3), (2, 4), (3, 5), (4, 6)}
        return
    frames = list(
        iter_file(mode, stream, dimension="length", units=units, cell_units=units)
    )
    assert len(frames) == 3
    for index, frame in enumerate(frames):
        np.testing.assert_allclose(frame["atoms"].q, beads.q[index], atol=0.001)
        np.testing.assert_allclose(frame["cell"].h, cell.h, atol=0.001)
        np.testing.assert_array_equal(frame["atoms"].names, beads.names)


def test_pdb_path_rejects_force_dimensions():
    with pytest.raises(ValueError, match="positions"):
        print_file_path("pdb", Beads(1, 2), Cell(), io.StringIO(), dimension="force")


@pytest.mark.parametrize(
    "encoded,match",
    [("0:,", "empty"), ("3;abc,", "delimiter"), ("3:abc;", "delimiter")],
)
def test_corrupt_netstrings_are_rejected(encoded, match):
    with pytest.raises(ValueError, match=match):
        netstring_encoded_loadz(io.StringIO(encoded))


@pytest.mark.parametrize("compressed", [False, True])
def test_netstring_roundtrip_preserves_multiple_matrices(tmp_path, compressed):
    # a text file, as for the planetary output (newlines are translated)
    expected = np.arange(-13.0, 50.0).reshape(7, 9) / 7
    with open(tmp_path / "matrices", "w") as stream:
        netstring_encoded_savez(stream, compressed=compressed, matrix=expected)
        netstring_encoded_savez(stream, compressed=compressed, matrix=expected.T)
    with open(tmp_path / "matrices") as stream:
        first, second = (netstring_encoded_loadz(stream)["matrix"] for _ in range(2))
    np.testing.assert_array_equal(first, expected)
    np.testing.assert_array_equal(second, expected.T)

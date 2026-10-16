"""Path trajectories preserve bead order, coordinates, and cell units."""

import io
import numpy as np
import pytest
from ipi.engine.beads import Beads
from ipi.engine.cell import Cell
from ipi.utils.units import unit_to_user
from ipi.utils.io import (
    print_file_path,
    iter_file,
    netstring_encoded_savez,
    netstring_encoded_loadz,
)


@pytest.mark.parametrize("mode", ["xyz", "pdb"])
@pytest.mark.parametrize("units", ["atomic_unit", "angstrom"])
def test_path_output_roundtrip(mode, units):
    beads = Beads(2, 3)
    beads.names[:] = ["H", "O"]
    beads.q[:] = np.arange(18).reshape(3, 6) / 10
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
@pytest.mark.xfail(
    strict=True,
    raises=TypeError,
    reason="netstring NPZ writer uses text StringIO for binary NumPy data",
)
def test_netstring_roundtrip_preserves_multiple_matrices(compressed):
    stream = io.StringIO()
    expected = np.arange(9).reshape(3, 3)
    netstring_encoded_savez(stream, compressed=compressed, matrix=expected)
    netstring_encoded_savez(stream, compressed=compressed, matrix=expected.T)
    stream.seek(0)
    np.testing.assert_array_equal(netstring_encoded_loadz(stream)["matrix"], expected)
    np.testing.assert_array_equal(netstring_encoded_loadz(stream)["matrix"], expected.T)

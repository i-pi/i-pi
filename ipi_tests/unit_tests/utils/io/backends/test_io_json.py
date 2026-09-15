"""Tests for JSON trajectory input and output."""

import io

import numpy as np

from ipi.engine.atoms import Atoms
from ipi.engine.cell import Cell
from ipi.utils.io.backends.io_json import iter_json, print_json, read_json


def test_json_roundtrip_preserves_configuration():
    """Checks positions, labels, cell and title after a JSON round-trip."""

    atoms = Atoms(2)
    atoms.q[:] = np.arange(6, dtype=float)
    atoms.names[:] = ["H", "O"]
    cell = Cell(np.diag([2.0, 3.0, 4.0]))
    stream = io.StringIO()

    print_json(atoms, cell, stream, title="frame 1")
    stream.seek(0)
    title, restored_cell, positions, names, masses = read_json(stream)

    assert title == "frame 1"
    assert np.allclose(restored_cell, cell.h)
    assert np.array_equal(positions, atoms.q)
    assert names.astype(str).tolist() == atoms.names.tolist()
    assert np.array_equal(masses, np.zeros(2))


def test_iter_json_reads_multiple_frames_and_stops_at_eof():
    """Checks iteration over every frame in a JSON trajectory."""

    atoms = Atoms(1)
    atoms.names[:] = ["H"]
    cell = Cell(np.eye(3))
    stream = io.StringIO()
    for frame in range(2):
        atoms.q[:] = frame
        print_json(atoms, cell, stream, title=f"frame {frame}")
    stream.seek(0)

    frames = list(iter_json(stream))

    assert [frame[0] for frame in frames] == ["frame 0", "frame 1"]
    assert np.array_equal(frames[0][2], np.zeros(3))
    assert np.array_equal(frames[1][2], np.ones(3))

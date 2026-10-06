from io import StringIO

import numpy as np
import pytest

from ipi.engine.atoms import Atoms
from ipi.engine.cell import Cell
from ipi.inputs.outputs import InputTrajectory
from ipi.utils.io import iter_file, print_file
from ipi.utils.io.inputs.io_xml import xml_parse_string
from ipi.utils.units import unit_to_internal, unit_to_user


def test_xyz_high_precision_round_trip():
    coordinates_angstrom = np.array(
        [
            30.123456789012345,
            -22.987654321098765,
            0.1234567890123456,
            31.234567890123456,
            -21.876543210987654,
            1.2345678901234567,
        ]
    )
    cell_angstrom = np.diag(
        [40.123456789012345, 41.234567890123456, 42.345678901234567]
    )
    atoms = Atoms(2)
    atoms.names = np.array(["H", "O"])
    cell = Cell(unit_to_internal("length", "angstrom", cell_angstrom))
    output = StringIO()

    expected = []
    for frame, displacement in enumerate([0.0, 1.0e-12]):
        frame_coordinates = coordinates_angstrom.copy()
        frame_coordinates[0] += displacement
        atoms.q[:] = unit_to_internal("length", "angstrom", frame_coordinates)
        expected.append(atoms.q.copy())
        print_file(
            "xyz_high_precision",
            atoms,
            cell,
            filedesc=output,
            title="frame=%d " % frame,
            key="positions",
            dimension="length",
            units="angstrom",
            cell_units="angstrom",
        )

    output.seek(0)
    frames = list(iter_file("xyz_high_precision", output))

    assert len(frames) == 2
    assert "frame=0" in frames[0]["comment"]
    assert "frame=1" in frames[1]["comment"]
    np.testing.assert_array_equal(frames[0]["atoms"].names, ["H", "O"])
    np.testing.assert_allclose(frames[0]["cell"].h, cell.h, rtol=0.0, atol=1.0e-13)
    np.testing.assert_allclose(
        frames[0]["atoms"].q, expected[0], rtol=0.0, atol=1.0e-13
    )
    np.testing.assert_allclose(
        frames[1]["atoms"].q, expected[1], rtol=0.0, atol=1.0e-13
    )
    difference_angstrom = unit_to_user(
        "length", "angstrom", frames[1]["atoms"].q[0] - frames[0]["atoms"].q[0]
    )
    assert difference_angstrom == pytest.approx(1.0e-12, rel=5.0e-3)


def test_input_trajectory_accepts_xyz_high_precision():
    root = xml_parse_string(
        "<trajectory format='xyz_high_precision'>positions</trajectory>"
    )
    trajectory = InputTrajectory()
    trajectory.parse(root.fields[0][1])
    assert trajectory.format.fetch() == "xyz_high_precision"

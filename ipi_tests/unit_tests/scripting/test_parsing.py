"""Tests for the output/trajectory parsers in ipi.scripting.parsing."""

import numpy as np
import pytest

from ipi.scripting import read_output, read_trajectory

OUTPUT_FILE = """\
# column   1     --> step : The current simulation time step.
# column   2     --> potential{electronvolt} : The potential energy.
# cols.    3-5   --> dipole : The electric dipole (x,y,z components).
   0.00000000e+00   1.00000000e+00   1.0   2.0   3.0
   1.00000000e+00   2.00000000e+00   4.0   5.0   6.0
"""

# minimal i-PI xyz frame: cubic 5 Ang cell, two atoms, positions in angstrom
XYZ_FILE = """\
2
# CELL{H}: 5.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 5.0 positions{angstrom} cell{angstrom}
H 0.0 0.0 0.0
H 1.0 0.0 0.0
"""


def test_read_output_values_and_info(tmp_path):
    f = tmp_path / "sim.out"
    f.write_text(OUTPUT_FILE)
    values, info = read_output(str(f))

    assert set(values) == {"step", "potential", "dipole"}
    assert np.allclose(values["step"], [0.0, 1.0])
    assert np.allclose(values["potential"], [1.0, 2.0])
    # multi-column property keeps its (nframes, 3) shape
    assert values["dipole"].shape == (2, 3)
    assert np.allclose(values["dipole"][1], [4.0, 5.0, 6.0])

    # units are parsed from the header, empty when absent
    assert info["potential"][0] == "electronvolt"
    assert info["step"][0] is None
    assert "potential energy" in info["potential"][1].lower()


def test_read_trajectory_returns_ase_atoms(tmp_path):
    pytest.importorskip("ase")
    f = tmp_path / "sim.xyz"
    f.write_text(XYZ_FILE)
    frames = read_trajectory(str(f))

    assert len(frames) == 1
    atoms = frames[0]
    assert list(atoms.symbols) == ["H", "H"]
    assert np.allclose(atoms.get_cell(), 5.0 * np.eye(3))
    assert np.allclose(atoms.positions[1], [1.0, 0.0, 0.0])


def test_read_trajectory_rejects_unknown_format(tmp_path):
    pytest.importorskip("ase")
    f = tmp_path / "sim.weird"
    f.write_text(XYZ_FILE)
    with pytest.raises(ValueError, match="Unrecognized file format"):
        read_trajectory(str(f))

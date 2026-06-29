"""Tests for the output/trajectory parsers in ipi.scripting.parsing."""

import numpy as np
import pytest

from ipi.scripting import read_output, read_trajectory
from ipi.utils.units import unit_to_user

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

# i-PI abuses the xyz format to store per-atom quantities other than positions:
# the data sits in the position columns and the comment names the quantity. Here
# it's the off-diagonal kinetic-energy estimator, in atomic units.
KINETIC_OD_FILE = """\
2
# CELL(abcABC): 5.0 5.0 5.0 90.0 90.0 90.0 Step: 5 Bead: 0 kinetic_od{atomic_unit} cell{atomic_unit}
O 0.1 0.2 0.3
H -0.4 0.5 -0.6
"""

# extended-xyz frame with several per-atom properties: positions plus a generic
# vector column (lands in `arrays`) and a forces column (lands on the calculator)
EXTXYZ_FILE = """\
2
Lattice="5.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 5.0" Properties=species:S:1:pos:R:3:myvec:R:3:forces:R:3 pbc="T T T"
H 0.0 0.0 0.0 0.1 0.2 0.3 1.0 2.0 3.0
H 1.0 0.0 0.0 -0.1 -0.2 -0.3 -1.0 -2.0 -3.0
"""

# `extras` files hold extra data returned by the calculator, one block per step
# tagged with `#EXTRAS(name)# Step: N`. Numeric blocks are stacked into an array.
EXTRAS_FILE = """\
 #EXTRAS(dipole)# Step:           0
  1.0  2.0  3.0
 #EXTRAS(dipole)# Step:           2
  4.0  5.0  6.0
"""

# when a block is not numeric, the reader falls back to raw per-step strings
EXTRAS_RAW_FILE = """\
 #EXTRAS(info)# Step:           0
 {"a": 1}
 #EXTRAS(info)# Step:           4
 {"a": 2}
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


def test_read_trajectory_named_quantity_not_zeroed(tmp_path):
    """A non-position quantity stored in xyz must keep its data and be unit
    converted. Regression: ASE's positions could alias the data buffer, so
    zeroing the positions in place wiped the quantity (it read all zeros)."""
    pytest.importorskip("ase")
    f = tmp_path / "sim.kod.xyz"
    f.write_text(KINETIC_OD_FILE)
    atoms = read_trajectory(str(f))[0]

    factor = unit_to_user("energy", "ase")  # atomic units -> ASE (eV)
    expected = np.array([[0.1, 0.2, 0.3], [-0.4, 0.5, -0.6]]) * factor
    assert "kinetic_od" in atoms.arrays
    np.testing.assert_allclose(atoms.arrays["kinetic_od"], expected)
    assert np.abs(atoms.arrays["kinetic_od"]).max() > 0.0  # guard vs zeroing

    # positions carry no meaningful data for such files and are set to zero
    assert np.abs(atoms.positions).max() == 0.0
    assert atoms.info["step"] == 5


def test_read_trajectory_extended_xyz_multiple_properties(tmp_path):
    """An extended-xyz file with several per-atom properties round-trips: the
    cell, positions, an extra vector array and the forces are all read."""
    pytest.importorskip("ase")
    f = tmp_path / "sim.extxyz"
    f.write_text(EXTXYZ_FILE)
    frames = read_trajectory(str(f))

    assert len(frames) == 1
    atoms = frames[0]
    assert np.allclose(atoms.get_cell(), 5.0 * np.eye(3))
    assert np.allclose(atoms.positions, [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    assert np.allclose(atoms.arrays["myvec"], [[0.1, 0.2, 0.3], [-0.1, -0.2, -0.3]])
    assert np.allclose(atoms.get_forces(), [[1.0, 2.0, 3.0], [-1.0, -2.0, -3.0]])


def test_read_trajectory_extras_numeric(tmp_path):
    """An `extras` file with numeric blocks returns the steps and the data
    stacked into a single array, keyed by the property name."""
    f = tmp_path / "sim.dipole"
    f.write_text(EXTRAS_FILE)
    data = read_trajectory(str(f), format="extras")

    assert set(data) == {"step", "dipole"}
    assert np.array_equal(data["step"], [0, 2])
    assert data["step"].dtype.kind == "i"
    assert np.allclose(data["dipole"], [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])


def test_read_trajectory_extras_raw_strings(tmp_path):
    """A non-numeric `extras` block falls back to one raw string per step."""
    f = tmp_path / "sim.info"
    f.write_text(EXTRAS_RAW_FILE)
    data = read_trajectory(str(f), format="extras")

    assert set(data) == {"step", "info"}
    assert np.array_equal(data["step"], [0, 4])
    assert len(data["info"]) == 2
    assert '{"a": 1}' in data["info"][0]
    assert '{"a": 2}' in data["info"][1]

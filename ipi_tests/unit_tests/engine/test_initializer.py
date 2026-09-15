#!/usr/bin/env python3
# pylint: disable=C0111,W0621,R0914,C0301
# +easier to find important problems

import pytest
import os

import numpy as np
import tempfile as tmp
import numpy.testing as npt
from types import SimpleNamespace

import ipi.engine.initializer as initializer
import ipi.utils.mathtools as mt

from ipi.inputs.simulation import InputSimulation
from ipi.utils.io.inputs.io_xml import xml_parse_string
from ipi.utils.units import Elements, unit_to_internal
from ipi.engine.atoms import Atoms
from ipi.engine.beads import Beads
from ipi.engine.cell import Cell


from ..common import xyz_generator as xyz_gen

# default_cell_mat = mt.abc2h(-1.0, -1.0, -1.0, np.pi/2.0, np.pi/2.0, np.pi/2.0) # After changing the input standard
default_cell_mat = np.eye(3) * -1
deg2rad = np.pi / 180.0


test_init_file_prms = [
    (1, 1, "asdasd", default_cell_mat, 1.0, 1.0),
    (2, 1, "asdasd", default_cell_mat, 1.0, 1.0),
    (1, 2, "asdasd", default_cell_mat, 1.0, 1.0),
    (10, 10, "dsadsa", default_cell_mat, 1.0, 1.0),
    (
        5,
        1,
        "asdasd positions{angstrom}  100 aaa cell{angstrom} asdasd ",
        default_cell_mat,
        1.8897261,
        1.8897261,
    ),
    (
        5,
        1,
        "asdasd positions{angstrom}  100 aaa cell{angstrom} 100 aaa # CELL(abcABC): 5.1 5.2 5.0 91.0  89  90 100 aaa asdasd ",
        mt.abc2h(5.1, 5.2, 5.0, 91 * deg2rad, 89 * deg2rad, 90 * deg2rad),
        1.8897261,
        1.8897261,
    ),
]


@pytest.fixture(params=test_init_file_prms)
def create_xyz_sample_file(request):
    """Create a fake xyz file and build the atoms and cell object from it."""

    (
        natoms,
        frames,
        comment,
        expected_cell,
        units_conv_at,
        units_conv_cell,
    ) = request.param

    filedesc, xyz, atoms_names = xyz_gen.xyz_traj_filedesc(natoms, frames, comment)

    # init_file needs to read from a real file...
    tmp_file = tmp.NamedTemporaryFile(mode="w", prefix="ipi_testing-tmp", delete=False)
    tmp_file.write(filedesc.read())
    tmp_file.close()
    filedesc.close()

    masses = np.zeros(natoms * frames)
    for _ii, _at in enumerate(atoms_names):
        masses[_ii] = Elements.mass(_at)

    ratoms = []
    for _fr in range(frames):
        ratoms.append(Atoms(natoms))
        ratoms[-1].q = xyz[_fr * natoms * 3 : 3 * (_fr + 1) * natoms] * units_conv_at
        ratoms[-1].m = masses[_fr * natoms : (_fr + 1) * natoms]
        ratoms[-1].names = atoms_names[_fr * natoms : (_fr + 1) * natoms]

    cell = Cell(expected_cell * units_conv_cell)

    # remove temp file created during this test
    def delete_files_after_testing():
        if os.path.isfile(tmp_file.name):
            os.remove(tmp_file.name)

    request.addfinalizer(delete_files_after_testing)
    return tmp_file, cell, ratoms


def test_init_file(create_xyz_sample_file):
    tmp_file, expected_cell, expected_ratoms = create_xyz_sample_file

    ret = initializer.init_file("xyz", tmp_file.name)

    for _ii, atoms in enumerate(ret[0]):
        npt.assert_array_almost_equal(expected_ratoms[_ii].q, atoms.q, 5)
        npt.assert_array_equal(expected_ratoms[_ii].names, atoms.names)
        npt.assert_array_almost_equal(expected_ratoms[_ii].m, atoms.m, 5)

    npt.assert_array_almost_equal(expected_cell.h, ret[1].h, 5)


def test_init_chk_skips_checkpoint_forcefields(tmp_path):
    """Checkpoint initialization reads state without constructing its PES."""

    checkpoint = tmp_path / "obsolete-forcefield.chk"
    checkpoint.write_text("""
<simulation>
  <ffdirect name='obsolete'>
    <pes>custom</pes>
  </ffdirect>
  <system>
    <forces><force forcefield='obsolete'/></forces>
    <motion mode='dynamics'>
      <dynamics mode='nve'>
        <timestep units='femtosecond'>0.5</timestep>
      </dynamics>
    </motion>
    <beads natoms='2' nbeads='2'>
      <q shape='(2, 6)' units='angstrom'>
        [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1]
      </q>
      <p shape='(2, 6)'>
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
      </p>
      <m shape='(2)'>[1, 2]</m>
      <names shape='(2)'>[H, He]</names>
    </beads>
    <cell shape='(3, 3)' units='angstrom'>
      [2, 0, 0, 0, 3, 0, 0, 0, 4]
    </cell>
  </system>
</simulation>
""")

    # The same checkpoint remains invalid as a full restart: FFDirect would
    # construct the custom PES and reject the missing pes_path.
    xmlchk = xml_parse_string(checkpoint.read_text())
    full_restart = InputSimulation()
    full_restart.parse(xmlchk.fields[0][1])
    with pytest.raises(ValueError, match="pes_path"):
        full_restart.fetch()

    beads, cell, motion = initializer.init_chk(checkpoint)
    q = np.arange(12, dtype=float).reshape(2, 6) / 10.0
    npt.assert_allclose(beads.q, q * unit_to_internal("length", "angstrom"))
    npt.assert_allclose(beads.p, np.arange(1, 13, dtype=float).reshape(2, 6))
    npt.assert_allclose(
        cell.h,
        np.diag([2.0, 3.0, 4.0]) * unit_to_internal("length", "angstrom"),
    )
    assert motion.dt == pytest.approx(unit_to_internal("time", "femtosecond", 0.5))

    # Exercise the regular <initialize> checkpoint code path as well.
    target = SimpleNamespace(beads=Beads(0, 0), cell=Cell())
    init = initializer.Initializer(
        nbeads=2,
        queue=[
            ("positions", initializer.InitIndexed(value=str(checkpoint), mode="chk")),
            ("masses", initializer.InitIndexed(value=str(checkpoint), mode="chk")),
            ("labels", initializer.InitIndexed(value=str(checkpoint), mode="chk")),
            ("momenta", initializer.InitIndexed(value=str(checkpoint), mode="chk")),
            ("cell", initializer.InitFile(value=str(checkpoint), mode="chk")),
        ],
    )
    init.init_stage1(target)
    npt.assert_allclose(target.beads.q, beads.q)
    npt.assert_allclose(target.beads.p, beads.p)
    npt.assert_allclose(target.cell.h, cell.h)

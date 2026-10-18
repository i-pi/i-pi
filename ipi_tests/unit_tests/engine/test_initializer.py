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


@pytest.fixture
def initialized_system():
    """A valid two-atom path; individual tests replace only the requested data."""
    beads = Beads(2, 4)
    beads.q[:] = np.arange(24).reshape(4, 6) / 10
    beads.p[:] = np.arange(1, 25).reshape(4, 6) / 10
    beads.m[:] = [2, 8]
    beads.names[:] = ["H", "He"]
    return SimpleNamespace(beads=beads, cell=Cell(np.eye(3) * 5))


@pytest.mark.parametrize("index", [-1, 1])
def test_mass_replacement_preserves_kinetic_energy(initialized_system, index):
    system = initialized_system
    original_p = system.beads.p.copy()
    kinetic = np.sum(original_p**2 / (2 * system.beads.m3), axis=1)
    masses = [8, 32] if index < 0 else [32]
    init = initializer.Initializer(
        4,
        [
            (
                "masses",
                initializer.InitIndexed(
                    value=np.array(masses), mode="manual", index=index
                ),
            )
        ],
    )
    init.init_stage1(system)
    expected_p = original_p.copy()
    expected_p[:, 0 if index < 0 else 3 :] *= 2
    npt.assert_allclose(system.beads.p, expected_p)
    npt.assert_allclose(
        np.sum(system.beads.p**2 / (2 * system.beads.m3), axis=1), kinetic
    )


@pytest.mark.parametrize("kind", ["positions", "momenta", "velocities"])
@pytest.mark.parametrize("index,bead", [(-1, -1), (1, -1), (-1, 2), (1, 2)])
def test_manual_vector_updates_only_selected_atoms_and_beads(
    initialized_system, kind, index, bead
):
    system = initialized_system
    target = system.beads.q if kind == "positions" else system.beads.p
    expected = target.copy()
    nbeads = 4 if bead < 0 else 1
    natoms = 2 if index < 0 else 1
    values = np.arange(1, nbeads * natoms * 3 + 1, dtype=float).reshape(nbeads, -1)
    replacement = values.copy()
    if kind == "velocities":
        replacement *= np.repeat([2, 8] if index < 0 else [8], 3)
    if kind != "positions":
        replacement *= np.sqrt(4 / nbeads)
    expected[
        slice(None) if bead < 0 else slice(bead, bead + 1),
        slice(None) if index < 0 else slice(3, 6),
    ] = replacement
    init = initializer.Initializer(
        4,
        [
            (
                kind,
                initializer.InitIndexed(
                    value=values.ravel().copy(), mode="manual", index=index, bead=bead
                ),
            )
        ],
    )
    init.init_stage1(system)
    npt.assert_allclose(target, expected)


@pytest.mark.parametrize("source_beads,target_beads", [(1, 4), (4, 1), (4, 8)])
def test_path_resizing_preserves_centroid(source_beads, target_beads):
    # A constant path represents the same physical configuration at any bead count.
    source = np.tile([1.0, -2.0, 3.0, 4.0, -5.0, 6.0], (source_beads, 1))
    target = np.zeros((target_beads, 6))
    initializer.set_vector(initializer.InitIndexed(), target, source)
    npt.assert_allclose(target, np.tile(source[0], (target_beads, 1)), atol=1e-14)
    # A nonconstant path contracted to one bead must become its centroid.
    if target_beads == 1:
        source[:, 0] = [0, 2, 4, 6]
        initializer.set_vector(initializer.InitIndexed(), target, source)
        assert target[0, 0] == pytest.approx(3)


@pytest.mark.parametrize(
    "index,bead,shape,message",
    [
        (-1, -1, (1, 3), "different atom numbers"),
        (2, -1, (1, 3), "atom index"),
        (-1, 4, (1, 6), "bead index"),
    ],
)
def test_invalid_vector_selection_does_not_modify_state(index, bead, shape, message):
    target = np.ones((4, 6))
    with pytest.raises(ValueError, match=message):
        initializer.set_vector(
            initializer.InitIndexed(index=index, bead=bead), target, np.zeros(shape)
        )
    npt.assert_array_equal(target, 1)


@pytest.mark.parametrize("kind", ["masses", "labels", "positions", "momenta"])
def test_indexed_initialization_requires_known_system_size(kind):
    target = SimpleNamespace(beads=Beads(0, 0), cell=Cell())
    init = initializer.Initializer(
        1, [(kind, initializer.InitIndexed(value=np.ones(3), mode="manual", index=0))]
    )
    with pytest.raises(ValueError, match="size of the system"):
        init.init_stage1(target)


@pytest.mark.parametrize("kind", ["masses", "labels"])
def test_mass_and_label_cannot_differ_between_beads(initialized_system, kind):
    value = np.array([1.0, 2.0]) if kind == "masses" else np.array(["H", "He"])
    init = initializer.Initializer(
        4, [(kind, initializer.InitIndexed(value=value, mode="manual", bead=1))]
    )
    with pytest.raises(ValueError, match="single bead"):
        init.init_stage1(initialized_system)


@pytest.mark.parametrize(
    "missing,message",
    [
        ("positions", "atomic positions"),
        ("cell", "cell"),
        ("masses", "masses"),
        ("labels", "labels"),
    ],
)
def test_missing_initial_state_is_rejected(initialized_system, missing, message):
    system = initialized_system
    if missing == "positions":
        system.beads = Beads(0, 0)
    elif missing == "cell":
        system.cell.h[:] = 0
    elif missing == "masses":
        system.beads.m[0] = 0
    else:
        system.beads.names[0] = ""
    with pytest.raises(ValueError, match=message):
        initializer.Initializer(4).init_stage1(system)


def test_manual_initialization_constructs_complete_system():
    system = SimpleNamespace(beads=Beads(0, 0), cell=Cell())
    queue = [
        (key, initializer.InitIndexed(value=np.array(value), mode="manual"))
        for key, value in [
            ("positions", [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]),
            ("masses", [2.0, 8.0]),
            ("labels", ["H", "He"]),
            ("cell", np.diag([4.0, 5.0, 6.0]).ravel()),
        ]
    ]
    initializer.Initializer(1, queue).init_stage1(system)
    npt.assert_array_equal(system.beads.q, [[0, 1, 2, 3, 4, 5]])
    npt.assert_array_equal(system.beads.p, 0)
    npt.assert_array_equal(system.beads.m, [2, 8])
    npt.assert_array_equal(system.beads.names, ["H", "He"])
    assert system.cell.V == pytest.approx(120)


@pytest.mark.parametrize("index", [-1, 1])
@pytest.mark.parametrize("temperature", [0.0, 300.0])
@pytest.mark.parametrize("nbeads", [1, 4])
def test_thermal_momenta_have_correct_mass_and_temperature_scaling(
    initialized_system, index, temperature, nbeads
):
    from ipi.engine.ensembles import Ensemble
    from ipi.engine.normalmodes import NormalModes
    from ipi.utils.prng import Random
    from ipi.utils.units import Constants

    system = initialized_system
    system.beads.resize(2, nbeads)
    system.beads.m[:] = [2, 8]
    system.beads.names[:] = ["H", "He"]
    system.beads.p[:] = 7
    system.ensemble = Ensemble(temp=150.0)
    system.nm = NormalModes()
    system.prng = Random(seed=321)
    selected = slice(None) if index < 0 else slice(3, 6)
    shape = (nbeads, 6 if index < 0 else 3)
    gaussian = Random(seed=321).gvec(shape)
    initializer.Initializer(
        nbeads,
        [
            (
                "momenta",
                initializer.InitIndexed(value=temperature, mode="thermal", index=index),
            )
        ],
    ).init_stage1(system)
    # Orthogonal normal-mode transforms preserve the mass-scaled squared norm.
    kinetic = np.sum(
        system.beads.p[:, selected] ** 2 / (2 * system.beads.m3[:, selected])
    )
    expected = (
        0.5 * nbeads * (temperature or 150.0) * Constants.kb * np.sum(gaussian**2)
    )
    assert kinetic == pytest.approx(expected, rel=1e-12)
    if index >= 0:
        npt.assert_array_equal(system.beads.p[:, :3], 7)


@pytest.mark.parametrize(
    "cell_comment", ["", "# CELL(abcABC): 4 5 6 90 90 90 cell{atomic_unit}"]
)
def test_xyz_initializes_positions_masses_labels_and_optional_cell(
    tmp_path, initialized_system, cell_comment
):
    path = tmp_path / "start.xyz"
    path.write_text(f"2\n{cell_comment} positions{{atomic_unit}}\nH 1 2 3\nHe 4 5 6\n")
    queue = [
        (key, initializer.InitIndexed(value=str(path), mode="xyz", units="automatic"))
        for key in ["positions", "masses", "labels", "cell"]
    ]
    initializer.Initializer(4, queue).init_stage1(initialized_system)
    npt.assert_allclose(initialized_system.beads.q, np.tile([1, 2, 3, 4, 5, 6], (4, 1)))
    npt.assert_allclose(
        initialized_system.beads.m, [Elements.mass("H"), Elements.mass("He")]
    )
    npt.assert_array_equal(initialized_system.beads.names, ["H", "He"])
    npt.assert_allclose(
        initialized_system.cell.h,
        np.diag([4, 5, 6]) if cell_comment else np.eye(3) * 5,
        atol=1e-14,
    )


@pytest.mark.xfail(
    strict=True,
    raises=UnboundLocalError,
    reason="Empty initialization file accesses undefined last-frame data",
)
def test_empty_initialization_file_has_clear_error(tmp_path):
    path = tmp_path / "empty.xyz"
    path.write_text("")
    with pytest.raises(ValueError, match="[Ee]mpty|[Ff]rame|[Cc]onfiguration"):
        initializer.init_file("xyz", path)


def test_whole_beads_cannot_be_initialized_from_manual_vector():
    with pytest.raises(ValueError, match="Cannot initialize manually"):
        initializer.init_beads(initializer.InitIndexed(mode="manual"), 1)

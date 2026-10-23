"""Initialization references: SI units, Fourier projection and fixed Gaussians."""

from types import SimpleNamespace

import numpy as np
import pytest

from ipi.engine.beads import Beads
from ipi.engine.cell import Cell
from ipi.engine.ensembles import Ensemble
from ipi.engine.initializer import (
    Initializer,
    InitIndexed,
    InitFile,
    init_file,
    init_beads,
    init_vector,
    set_vector,
)
from ipi.engine.normalmodes import NormalModes
from ipi.utils.prng import Random
from ipi_tests.unit_tests.engine.mathematical_references import (
    fourier_matrix,
    mode_mass_factors,
    projection,
)
from ipi_tests.unit_tests.engine.test_property_mathematical_references import verifies
from ipi_tests.unit_tests.reference_units import (
    ANGSTROM_TO_BOHR,
    BOHR_IN_ANGSTROM,
    HARTREE_IN_JOULE,
    UNIT_RTOL,
)


@verifies("initializer:init_file", "initializer:init_beads", "initializer:init_vector")
@pytest.mark.parametrize("format", ["xyz", "pdb", "ase"])
def test_initialization_formats_have_known_geometry_and_independent_atomic_masses(
    tmp_path, format
):
    # Historical CIAAW standard weights, https://ciaaw.org/pubs/TSAW-1997.pdf;
    # CODATA 2022 unified atomic mass / electron mass = 1822.888486209.
    expected_mass = np.array([1.00794, 4.002602]) * 1822.888486209
    q = np.array([[0.1, -0.2, 0.3], [0.4, 0.5, -0.6]])
    if format == "pdb":
        file = tmp_path / "start.pdb"
        # Literal fixed-width PDB records; coordinates and orthogonal cell in Angstrom.
        file.write_text(
            "CRYST1   10.000   12.000   14.000  90.00  90.00  90.00 P 1           1\n"
            + "ATOM      1  H   GAS A   1       0.100  -0.200   0.300  1.00  0.00           H\n"
            + "ATOM      2  He  GAS A   1       0.400   0.500  -0.600  1.00  0.00          He\nEND\n"
        )
    else:
        file = tmp_path / "start.xyz"
        file.write_text(
            "2\n# CELL(abcABC): 10 12 14 90 90 90 positions{angstrom} cell{angstrom}\nH .1 -.2 .3\nHe .4 .5 -.6\n"
        )
    mode = "ase" if format == "ase" else format
    atoms, cell = init_file(mode, str(file), units="angstrom", cell_units="angstrom")
    assert len(atoms) == 1
    np.testing.assert_allclose(
        atoms[0].q, q.ravel() * ANGSTROM_TO_BOHR, rtol=UNIT_RTOL, atol=1e-12
    )
    # ASE supplies its own modern element masses, so check it separately below.
    if format != "ase":
        np.testing.assert_allclose(atoms[0].m, expected_mass, rtol=UNIT_RTOL)
    else:
        np.testing.assert_allclose(
            atoms[0].m, np.array([1.008, 4.002602]) * 1822.888486209, rtol=UNIT_RTOL
        )
    np.testing.assert_array_equal(atoms[0].names, ["H", "He"])
    if format != "ase":
        np.testing.assert_allclose(
            cell.h, np.diag([10, 12, 14]) * ANGSTROM_TO_BOHR, rtol=UNIT_RTOL, atol=1e-12
        )
    bead_data = init_beads(InitFile(value=str(file), mode=mode), 1, units="angstrom")
    np.testing.assert_allclose(
        bead_data.q[0], q.ravel() * ANGSTROM_TO_BOHR, rtol=UNIT_RTOL
    )
    vector = init_vector(InitIndexed(value=str(file), mode=mode), 1, units="angstrom")
    np.testing.assert_allclose(vector[0], q.ravel() * ANGSTROM_TO_BOHR, rtol=UNIT_RTOL)


@verifies(
    "initializer:init_vector",
    "initializer:set_vector",
    "initializer:Initializer.init_stage1",
)
@pytest.mark.parametrize("source,target", [(3, 5), (5, 3), (5, 7), (7, 5)])
@pytest.mark.parametrize("kind", ["positions", "momenta", "velocities"])
def test_nonconstant_path_initialization_matches_independent_fourier_projection(
    tmp_path, source, target, kind
):
    values = (
        0.2 * np.sin(np.arange(source * 6).reshape(source, 6) * 0.43)
        + np.arange(6) * 0.1
    )
    matrix = projection(source, target)
    expected = matrix @ values
    if kind == "positions":
        target_q = np.zeros((target, 6))
        set_vector(InitIndexed(), target_q, values)
        np.testing.assert_allclose(target_q, expected, atol=1e-13)
    system = SimpleNamespace(beads=Beads(2, target), cell=Cell(np.eye(3) * 10))
    system.beads.m[:] = [2, 5]
    system.beads.names[:] = ["H", "He"]
    system.beads.q[:] = 1
    if kind == "velocities":
        expected *= np.repeat([2.0, 5.0], 3)
    if kind != "positions":
        expected *= np.sqrt(target / source)
    file = tmp_path / "source.xyz"
    frames = []
    for frame in values:
        frames += ["2", kind + "{atomic_unit}"]
        for name, xyz in zip(["H", "He"], frame.reshape(2, 3)):
            frames.append(name + " " + " ".join(repr(float(x)) for x in xyz))
    file.write_text("\n".join(frames) + "\n")
    Initializer(
        target, [(kind, InitIndexed(mode="xyz", value=str(file), units="automatic"))]
    ).init_stage1(system)
    np.testing.assert_allclose(
        system.beads.q if kind == "positions" else system.beads.p, expected, atol=1e-12
    )


@verifies("initializer:Initializer.init_stage1")
@pytest.mark.parametrize(
    "style,frequencies",
    [("manual", [0.3, 0.5, 0.7]), ("pa-cmd", [0.4]), ("wmax-cmd", [0.8, 0.3])],
)
@pytest.mark.parametrize("index", [-1, 1])
@pytest.mark.parametrize("thermal_temperature", [0.3, 0.0, -0.2])
def test_thermal_initialization_with_nonphysical_mode_masses_uses_fixed_draws(
    monkeypatch, style, frequencies, index, thermal_temperature
):
    beads = Beads(2, 4)
    beads.m[:] = [2, 5]
    beads.names[:] = ["H", "He"]
    beads.q[:] = 0.1
    beads.p[:] = 7
    system = SimpleNamespace(
        beads=beads,
        cell=Cell(np.eye(3) * 10),
        ensemble=Ensemble(temp=0.4),
        nm=NormalModes(mode=style, freqs=frequencies, transform_method="matrix"),
        prng=Random(seed=1),
    )
    count = 6 if index < 0 else 3
    draws = np.linspace(-1.2, 1.4, 4 * count).reshape(4, count)

    def fixed(shape):
        assert shape == draws.shape
        return draws.copy()

    monkeypatch.setattr(system.prng, "gvec", fixed)
    factors = mode_mass_factors(4, 0.4, style, frequencies)
    masses = np.repeat([2, 5] if index < 0 else [5], 3)
    temperature = thermal_temperature if thermal_temperature > 0 else 0.4
    expected = fourier_matrix(4).T @ (
        draws * np.sqrt(4 * temperature * factors[:, None] * masses)
    )
    Initializer(
        4,
        [
            (
                "momenta",
                InitIndexed(value=thermal_temperature, mode="thermal", index=index),
            )
        ],
    ).init_stage1(system)
    selected = slice(None) if index < 0 else slice(3, 6)
    np.testing.assert_allclose(beads.p[:, selected], expected, atol=1e-12)
    if index >= 0:
        np.testing.assert_array_equal(beads.p[:, :3], 7)


@verifies(
    "initializer:Initializer.init_stage1",
    "initializer:init_file",
    "initializer:init_vector",
)
def test_initial_mass_cell_and_velocity_units_use_independent_si_anchors(tmp_path):
    # Atomic velocity a0*Eh/hbar, using exact SI h and independent CODATA anchors.
    velocity_factor = (6.62607015e-34 / (2 * np.pi)) / (
        HARTREE_IN_JOULE * BOHR_IN_ANGSTROM * 1e-10
    )
    values = np.array([1e5, -2e5, 3e5, -1e5, 2e5, -3e5])
    file = tmp_path / "velocity.xyz"
    file.write_text(
        "2\nvelocities{m/s}\nH 100000 -200000 300000\nHe -100000 200000 -300000\n"
    )
    system = SimpleNamespace(beads=Beads(2, 1), cell=Cell())
    system.beads.m[:] = [1, 1]
    system.beads.names[:] = ["H", "He"]
    system.beads.q[:] = 1
    queue = [
        ("masses", InitIndexed(mode="manual", value=np.array([2, 5]), units="dalton")),
        (
            "cell",
            InitIndexed(
                mode="manual",
                value=np.array([6, 1, 0.5, 0, 5, 0.7, 0, 0, 4]),
                units="angstrom",
            ),
        ),
        ("velocities", InitIndexed(mode="xyz", value=str(file), units="m/s")),
    ]
    Initializer(1, queue).init_stage1(system)
    mass = np.array([2, 5]) * 1822.888486209
    np.testing.assert_allclose(system.beads.m, mass, rtol=UNIT_RTOL)
    np.testing.assert_allclose(
        system.beads.p[0], values * velocity_factor * np.repeat(mass, 3), rtol=UNIT_RTOL
    )
    np.testing.assert_allclose(
        system.cell.h,
        np.array([[6, 1, 0.5], [0, 5, 0.7], [0, 0, 4]]) * ANGSTROM_TO_BOHR,
        rtol=UNIT_RTOL,
    )


@verifies("initializer:init_vector")
def test_momentum_file_preserves_literal_atomic_momenta(tmp_path):
    file = tmp_path / "momentum.xyz"
    file.write_text("2\nmomenta{atomic_unit}\nH .1 -.2 .3\nHe -.4 .5 -.6\n")
    vector = init_vector(
        InitIndexed(value=str(file), mode="xyz"),
        1,
        momenta=True,
        dimension="momentum",
        units="atomic_unit",
    )
    np.testing.assert_array_equal(vector, [[0.1, -0.2, 0.3, -0.4, 0.5, -0.6]])

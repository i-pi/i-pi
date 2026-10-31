"""Real PLUMED restraint energies, forces and update accounting."""

import sys
from types import SimpleNamespace

import numpy as np
import pytest

import ipi.utils.array_backend as ab
from ipi.engine.forcefields import FFPlumed, ForceRequest
from ipi.utils.softexit import Softexit


@pytest.fixture
def plumed_field(tmp_path, monkeypatch):
    plumed = pytest.importorskip("plumed")
    try:
        probe = plumed.Plumed()
    except RuntimeError as error:
        pytest.skip(f"PLUMED kernel unavailable: {error}")
    probe.finalize()
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("ipi.engine.forcefields.softexit", Softexit())
    structure = tmp_path / "init.xyz"
    structure.write_text("2\n\nH 0 0 0\nH 1 0 0\n")
    bias = tmp_path / "plumed.dat"
    bias.write_text(
        "UNITS LENGTH=Bohr ENERGY=Ha\n"
        "d: DISTANCE ATOMS=1,2 NOPBC\n"
        "r: RESTRAINT ARG=d AT=0 KAPPA=2\n"
    )
    fields = []

    def create(mode="xyz", **options):
        ff = FFPlumed(
            init_file=SimpleNamespace(
                mode=mode, value=str(structure), units="atomic_unit"
            ),
            plumed_dat=str(bias),
            **options,
        )
        fields.append(ff)
        return ff

    yield create
    for ff in fields:
        ff.softexit()


@pytest.mark.parametrize(
    "threaded,step,compute_work", [(False, 0, True), (True, 3, False)]
)
def test_plumed_distance_restraint_matches_quadratic(
    plumed_field, threaded, step, compute_work
):
    ff = plumed_field(
        threaded=threaded,
        plumed_step=step,
        compute_work=compute_work,
        plumed_extras=["d"],
    )
    assert not ff.threaded
    assert ff.plumedrestart == (step > 0)
    positions = np.array([0.0, 0.0, 0.0, 1.0, 2.0, -1.0])
    cell = np.eye(3) * 20
    request = ForceRequest(pos=positions, cell=(cell, np.linalg.inv(cell)))
    ff.evaluate(request)
    energy, force, virial, extra = request["result"]
    assert energy == pytest.approx(6, rel=1e-6)
    np.testing.assert_allclose(force, [2, 4, -2, -2, -4, 2], rtol=1e-6)
    np.testing.assert_allclose(
        virial, -2 * np.outer(positions[3:], positions[3:]), rtol=1e-6
    )
    assert extra["raw"] == ""
    np.testing.assert_allclose(extra["d"], [np.sqrt(6)], rtol=1e-6)
    assert request._event_done.is_set()
    assert ff.mtd_update(positions, cell) == pytest.approx(0, abs=1e-12)
    assert ff.plumed_step == step + 1
    with pytest.raises(ValueError, match="Size of atom array"):
        ff.evaluate(ForceRequest(pos=np.zeros(3)))


def bias_with_physical_forces(ff):
    """Results without and with a separate harmonic pair potential U = |r|^2/2."""
    q = np.array([0.0, 0.0, 0.0, 1.0, 2.0, -1.0])
    cell = np.eye(3) * 20
    first = ForceRequest(pos=q, cell=(cell, np.linalg.inv(cell)))
    ff.evaluate(first)
    # the physical forces are arrays of the backend
    system = ab.xp.asarray(q)
    ff.system_force = SimpleNamespace(
        pot=ab.xp.asarray(3.0),
        f=ab.xp.concat([system[3:], -system[3:]]),
        vir=-ab.xp.linalg.outer(system[3:], system[3:]),
    )
    second = ForceRequest(pos=q, cell=(cell, np.linalg.inv(cell)))
    ff.evaluate(second)
    return first["result"], second["result"]


def test_plumed_returns_only_bias_energy_and_forces(plumed_field):
    alone, combined = bias_with_physical_forces(plumed_field())
    for component in range(2):
        np.testing.assert_allclose(combined[component], alone[component], atol=1e-10)


def test_plumed_returns_only_bias_virial(plumed_field):
    alone, combined = bias_with_physical_forces(plumed_field())
    np.testing.assert_allclose(combined[2], alone[2], atol=1e-10)


def test_plumed_update_recomputes_bias_after_positions_change(plumed_field):
    ff = plumed_field()
    q = np.array([0.0, 0.0, 0.0, 1.0, 0.0, 0.0])
    cell = np.eye(3) * 20
    ff.evaluate(ForceRequest(pos=q, cell=(cell, np.linalg.inv(cell))))
    # A static restraint does no update work, even if positions have moved.
    assert ff.mtd_update(2 * q, cell) == pytest.approx(0, abs=1e-12)
    np.testing.assert_array_equal(ff.lastq, 2 * q)


def test_plumed_with_torch_arrays(plumed_field, strict_torch):
    """The initial structure, the physical forces and the positions passed to
    mtd_update can be arrays of the torch backend"""
    expected = bias_with_physical_forces(plumed_field())
    ab.set_backend("torch")
    try:
        ff = plumed_field()
        for results, reference in zip(bias_with_physical_forces(ff), expected):
            for value, ref in zip(results[:3], reference[:3]):
                np.testing.assert_allclose(value, ref, atol=1e-10)
        q = ab.xp.asarray([0.0, 0.0, 0.0, 2.0, 0.0, 0.0])
        assert ff.mtd_update(q, 20 * ab.xp.eye(3)) == pytest.approx(0, abs=1e-12)
        np.testing.assert_array_equal(ff.lastq, ab.to_numpy(q))
    finally:
        ab.set_backend("numpy")


def test_plumed_rejects_unsupported_initial_structure_format(plumed_field):
    with pytest.raises(ValueError, match="Unsupported init file format"):
        plumed_field(mode="unsupported")


def test_plumed_missing_library_reports_import_error(monkeypatch):

    monkeypatch.setitem(sys.modules, "plumed", None)
    with pytest.raises(ImportError, match="Cannot find plumed"):
        FFPlumed()

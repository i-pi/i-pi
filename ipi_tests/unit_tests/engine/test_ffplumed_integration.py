"""Numerical checks of FFPlumed against a real, optional PLUMED kernel."""

from types import SimpleNamespace

import numpy as np
import pytest

from ipi.engine.forcefields import FFPlumed, ForceRequest
from ipi.utils.softexit import softexit

pytest.importorskip("plumed")

# CODATA values, independent of FFPlumed's rounded conversion constants.
BOHR_NM = 0.0529177210544
HARTREE_KJ_MOL = 2625.4996394799


@pytest.fixture
def make_forcefield(tmp_path, monkeypatch):
    """Construct native interfaces without leaving global exit callbacks."""
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(softexit, "register_function", lambda function: None)
    xyz = tmp_path / "atoms.xyz"
    xyz.write_text("2\n\nH 0 0 0\nH 3 1 0.5\n")
    instances = []

    def create(lines, extras=(), step=0):
        """Create one forcefield with an independent PLUMED input."""
        dat = tmp_path / f"plumed-{len(instances)}.dat"
        dat.write_text(lines)
        ff = FFPlumed(
            init_file=SimpleNamespace(mode="xyz", value=str(xyz), units="atomic_unit"),
            plumed_dat=str(dat),
            plumed_extras=list(extras),
            plumed_step=step,
        )
        instances.append(ff)
        return ff

    yield create
    for ff in instances:
        ff.plumed.finalize()


def evaluate(ff, positions, cell, physical):
    """Provide a harmonic physical potential when ENERGY derivatives need it."""
    if physical:
        ff.system_force = SimpleNamespace(
            pot=0.5 * np.sum(positions**2),
            f=-positions.flatten(),
            vir=-positions.T @ positions,
        )
    request = ForceRequest(
        {
            "pos": positions.flatten(),
            "cell": (cell, np.linalg.inv(cell)),
            "result": None,
        }
    )
    ff.evaluate(request)
    return request["result"]


@pytest.mark.parametrize("mode", ["distance", "distance_with_physical", "energy"])
def test_real_plumed_bias_derivatives(make_forcefield, mode):
    """Check energy conversion and Cartesian/strain derivatives of physical bias."""
    if mode == "energy":
        lines = "e: ENERGY\nb: RESTRAINT ARG=e AT=0 KAPPA=0.0001\n"
    else:
        lines = "d: DISTANCE ATOMS=1,2 NOPBC\nb: RESTRAINT ARG=d AT=0.12 KAPPA=1000\n"
    ff = make_forcefield(lines)
    positions = np.array([[0.2, -0.1, 0.3], [3.0, 1.0, 0.5]])
    cell = np.array([[15.0, 1.0, 0.5], [0.0, 14.0, 0.7], [0.0, 0.0, 13.0]])
    physical = mode != "distance"
    energy, force, virial, _ = evaluate(ff, positions, cell, physical)
    if mode == "energy":
        value = 0.5 * np.sum(positions**2) * HARTREE_KJ_MOL
        expected = 0.5 * 0.0001 * value**2 / HARTREE_KJ_MOL
    else:
        value = np.linalg.norm(positions[1] - positions[0]) * BOHR_NM
        expected = 0.5 * 1000 * (value - 0.12) ** 2 / HARTREE_KJ_MOL
    np.testing.assert_allclose(energy, expected, rtol=1e-7, atol=1e-10)
    step = 1e-5
    force_fd = np.zeros_like(positions)
    virial_fd = np.zeros((3, 3))
    for index in np.ndindex(positions.shape):
        delta = np.zeros_like(positions)
        delta[index] = step
        plus = evaluate(ff, positions + delta, cell, physical)[0]
        minus = evaluate(ff, positions - delta, cell, physical)[0]
        force_fd[index] = -(plus - minus) / (2 * step)
    for index in np.ndindex(virial.shape):
        strain = np.zeros((3, 3))
        strain[index] = step
        plus = evaluate(
            ff,
            positions @ (np.eye(3) + strain),
            (np.eye(3) + strain).T @ cell,
            physical,
        )[0]
        minus = evaluate(
            ff,
            positions @ (np.eye(3) - strain),
            (np.eye(3) - strain).T @ cell,
            physical,
        )[0]
        virial_fd[index] = -(plus - minus) / (2 * step)
    np.testing.assert_allclose(force.reshape(-1, 3), force_fd, rtol=1e-7, atol=1e-9)
    np.testing.assert_allclose(virial, virial_fd, rtol=1e-7, atol=1e-9)


def test_real_plumed_scalar_extra(make_forcefield):
    """Read a rank-zero CV and keep returned extras independent of later calls."""
    ff = make_forcefield("d: DISTANCE ATOMS=1,2 NOPBC\n", extras=["d"])
    positions = np.array([[0.0, 0.0, 0.0], [3.0, 1.0, 0.5]])
    first = evaluate(ff, positions, np.eye(3) * 15, False)[3]["d"]
    expected = np.linalg.norm(positions[1]) * BOHR_NM
    np.testing.assert_allclose(first, expected, rtol=1e-7)
    evaluate(ff, positions * 1.1, np.eye(3) * 15, False)
    np.testing.assert_allclose(first, expected, rtol=1e-7)


def test_real_opes_updates_and_state_restore(make_forcefield, tmp_path):
    """Separate force evaluations from OPES updates and restore the saved bias."""
    state = tmp_path / "opes.state"
    lines = (
        "d: DISTANCE ATOMS=1,2 NOPBC\n"
        "b: OPES_METAD ARG=d SIGMA=0.01 PACE=1 BARRIER=10 TEMP=300 "
        f"STATE_WFILE={state.name} STATE_WSTRIDE=1\n"
        "FLUSH STRIDE=1\n"
    )
    ff = make_forcefield(lines)
    positions = np.array([[0.0, 0.0, 0.0], [3.0, 1.0, 0.5]])
    cell = np.eye(3) * 15
    counters = []
    for step in range(1, 5):
        before = state.read_bytes() if state.exists() else None
        for _ in range(3):
            evaluate(ff, positions, cell, False)
        assert (state.read_bytes() if state.exists() else None) == before
        assert np.isfinite(ff.mtd_update(positions.flatten(), cell))
        assert ff.plumed_step == step
        if state.exists():
            counter_lines = [
                line
                for line in state.read_text().splitlines()
                if line.startswith("#! SET counter ")
            ]
            if counter_lines:
                assert len(counter_lines) == 1
                counters.append(int(counter_lines[0].split()[-1]))
            else:
                # Opening the output creates an empty file before the first
                # (intentionally skipped) OPES update writes any state.
                assert step == 1
                assert state.stat().st_size == 0
        positions[1, 0] += 0.1
    # PLUMED skips its first update; each subsequent single-walker update
    # adds exactly one observation, regardless of the force-evaluation count.
    assert len(counters) == 3
    np.testing.assert_array_equal(np.diff(counters), [1, 1])
    expected = evaluate(ff, positions, cell, False)
    restart_lines = (
        "d: DISTANCE ATOMS=1,2 NOPBC\n"
        "b: OPES_METAD ARG=d SIGMA=0.01 PACE=1 BARRIER=10 TEMP=300 "
        f"STATE_RFILE={state.name}\n"
    )
    restored = make_forcefield(restart_lines, step=ff.plumed_step)
    actual = evaluate(restored, positions, cell, False)
    for result, reference in zip(actual[:3], expected[:3], strict=True):
        np.testing.assert_allclose(result, reference, rtol=1e-10, atol=1e-10)

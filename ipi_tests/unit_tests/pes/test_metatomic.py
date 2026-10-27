"""Metatomic driver checks with PET-MAD-xs: batching, derivatives, UQ, rotations.

The exported model is read from ``IPI_METATOMIC_MODEL``; otherwise it is exported
with ``upet``, and the module is skipped when neither is available.
"""

import json
import os
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("metatomic.torch")
pytest.importorskip("vesin.metatomic")
import ase.io
from ase import Atoms as AseAtoms
from metatomic.torch.ase_calculator import MetatomicCalculator

from ipi.engine.atoms import Atoms
from ipi.engine.cell import Cell
from ipi.engine.forcefields import FFDirect, FFRotations
from ipi.pes.metatomic import MetatomicDriver
from ipi.utils.mathtools import random_rotation
from ipi.utils.prng import Random
from ipi.utils.units import unit_to_internal, unit_to_user

WATER_32 = Path(__file__).resolve().parents[3] / "examples/init_files/water_32.pdb"


def ev(x):
    return unit_to_internal("energy", "electronvolt", x)


def ev_ang(x):
    return unit_to_internal("force", "ev/ang", x)


@pytest.fixture(scope="module")
def upet_model(tmp_path_factory):
    path = os.environ.get("IPI_METATOMIC_MODEL")
    if path:
        return path
    upet = pytest.importorskip("upet")
    path = str(tmp_path_factory.mktemp("upet") / "pet-mad-xs-v1.6.0.pt")
    upet.save_upet(model="pet-mad", size="xs", version="1.6.0", output=path)
    return path


def _write_structure(tmp_path_factory, natoms):
    """First ``natoms`` atoms of the 32-water box, kept in the periodic cell."""
    water = ase.io.read(WATER_32)[:natoms]
    atoms = AseAtoms(
        water.get_chemical_symbols(),
        positions=water.positions,
        cell=water.cell,
        pbc=True,
    )
    path = tmp_path_factory.mktemp("structure") / f"water_{natoms}.xyz"
    ase.io.write(path, atoms)
    cell = unit_to_internal("length", "angstrom", atoms.cell.array.T)
    pos = unit_to_internal("length", "angstrom", atoms.positions)
    return str(path), atoms, cell, pos


@pytest.fixture(scope="module")
def water(tmp_path_factory):
    return _write_structure(tmp_path_factory, 12)


@pytest.fixture(scope="module")
def driver(upet_model, water):
    return MetatomicDriver(water[0], upet_model)


@pytest.mark.parametrize("device", ["cpu", "cpu;cpu"])
def test_batched_matches_single(upet_model, water, device):
    template, _, cell, pos = water
    driver = MetatomicDriver(template, upet_model, device=device)
    shift = unit_to_internal("length", "angstrom", 0.05)
    positions = [
        pos + shift * np.sin(k + np.arange(pos.size)).reshape(pos.shape)
        for k in range(3)
    ]
    cells = [cell] * len(positions)
    try:
        batch = driver.compute(cells, positions)
        assert len(batch) == len(positions)
        for c, q, result in zip(cells, positions, batch):
            single = driver.compute(c, q)
            assert result[0] == pytest.approx(single[0], abs=ev(1e-5))
            np.testing.assert_allclose(result[1], single[1], atol=ev_ang(1e-4))
            np.testing.assert_allclose(result[2], single[2], atol=ev(1e-4))
        assert batch[0][0] != pytest.approx(batch[1][0], abs=ev(1e-5))
        # repeated batches run the neighbor lists concurrently on all devices
        for _ in range(15):
            repeated = driver.compute(cells + cells[:1], positions + positions[:1])
        assert repeated[3][0] == pytest.approx(batch[0][0], abs=ev(1e-5))
    finally:
        if driver.executor is not None:
            driver.executor.shutdown()


def test_forces_and_virial_match_finite_differences(driver, water):
    _, _, cell, pos = water
    _, force, virial, _ = driver.compute(cell, pos)
    h = unit_to_internal("length", "angstrom", 5e-3)
    for i, a in [(0, 0), (1, 2), (5, 1)]:
        dq = np.zeros_like(pos)
        dq[i, a] = h
        e_plus = driver.compute(cell, pos + dq)[0]
        e_minus = driver.compute(cell, pos - dq)[0]
        assert -(e_plus - e_minus) / (2 * h) == pytest.approx(
            force[i, a], abs=ev_ang(5e-3)
        )
    eps = 1e-3
    for a in range(3):
        strain = np.eye(3)
        strain[a, a] += eps
        e_plus = driver.compute(strain @ cell, pos @ strain.T)[0]
        strain[a, a] -= 2 * eps
        e_minus = driver.compute(strain @ cell, pos @ strain.T)[0]
        assert -(e_plus - e_minus) / (2 * eps) == pytest.approx(
            virial[a, a], abs=ev(1e-2)
        )


@pytest.mark.parametrize("non_conservative", [False, True])
def test_matches_metatomic_ase_calculator(upet_model, water, non_conservative):
    template, atoms, cell, pos = water
    driver = MetatomicDriver(template, upet_model, non_conservative=non_conservative)
    energy, force, virial, _ = driver.compute(cell, pos)
    atoms = atoms.copy()
    atoms.calc = MetatomicCalculator(upet_model, non_conservative=non_conservative)
    assert energy == pytest.approx(ev(atoms.get_potential_energy()), abs=ev(1e-4))
    forces = atoms.get_forces()
    if non_conservative:
        forces -= forces.mean(axis=0)
    np.testing.assert_allclose(force, ev_ang(forces), atol=ev_ang(1e-4))
    stress_virial = -atoms.get_stress(voigt=False) * atoms.get_volume()
    np.testing.assert_allclose(virial, ev(stress_virial), atol=ev(1e-3))


def test_non_conservative_forces_close_to_conservative(upet_model, water, driver):
    template, _, cell, pos = water
    force = driver.compute(cell, pos)[1]
    nc_driver = MetatomicDriver(template, upet_model, non_conservative=True)
    nc_force = nc_driver.compute(cell, pos)[1]
    np.testing.assert_allclose(nc_force.sum(axis=0), 0.0, atol=1e-12)
    rms = lambda x: np.sqrt(np.mean(x**2))
    assert rms(nc_force - force) < 0.35 * rms(force)


def test_per_atom_uncertainty(upet_model, water):
    template, _, cell, pos = water
    driver = MetatomicDriver(template, upet_model, uncertainty_threshold=1e-6)
    extras = json.loads(driver.compute(cell, pos)[3])
    assert len(extras["energy_uncertainty"]) == len(pos)
    assert min(extras["energy_uncertainty"]) > 0.0


def test_energy_ensemble_is_centred_on_energy(upet_model, water):
    template, _, cell, pos = water
    driver = MetatomicDriver(template, upet_model, energy_ensemble=True)
    energy, _, _, raw = driver.compute(cell, pos)
    committee = json.loads(raw)["committee_pot"]
    assert len(committee) > 1
    assert np.std(committee) > 0.0
    assert np.mean(committee) == pytest.approx(energy, abs=ev(1e-3))


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="BUG (upstream LLPR): requesting energy_ensemble makes energy_uncertainty per-structure",
)
def test_per_atom_uncertainty_with_energy_ensemble(upet_model, water):
    template, _, cell, pos = water
    driver = MetatomicDriver(
        template, upet_model, uncertainty_threshold=1e-6, energy_ensemble=True
    )
    extras = json.loads(driver.compute(cell, pos)[3])
    assert len(extras["energy_uncertainty"]) == len(pos)


def test_force_virial_ensemble(upet_model, tmp_path_factory):
    template, _, cell, pos = _write_structure(tmp_path_factory, 3)
    driver = MetatomicDriver(
        template, upet_model, energy_ensemble=True, force_virial_ensemble=True
    )
    _, force, virial, raw = driver.compute(cell, pos)
    extras = json.loads(raw)
    n_members = len(extras["committee_pot"])
    committee_force = np.array(extras["committee_force"]).reshape(n_members, -1, 3)
    committee_virial = np.array(extras["committee_virial"]).reshape(n_members, 3, 3)
    assert committee_force.shape[1] == len(pos)
    np.testing.assert_allclose(committee_force.mean(axis=0), force, atol=ev_ang(1e-3))
    np.testing.assert_allclose(committee_virial.mean(axis=0), virial, atol=ev(1e-3))


def test_rotation_average_matches_manual_average(upet_model, water, driver):
    template, _, cell, pos = water
    ff = FFRotations(
        ffdirect=FFDirect(
            name="pet",
            pes="metatomic",
            pars={"template": template, "model": upet_model},
            threaded=False,
            batch_size=18,
        ),
        grid_mode="lebedev",
        grid_order=3,
        prng=Random(seed=17),
        dopbc=False,
    )
    assert len(ff._rotations) == 18
    atoms = Atoms(len(pos))
    atoms.q_[:] = pos.flatten()
    request = ff.queue(atoms, Cell(cell))
    ff.poll()
    assert request["status"] == "Done"
    energy, force, virial, extras = request["result"]

    energies, forces, weights = [], [], []
    for R, w, _ in ff._rotations:
        e, f, _, _ = driver.compute(R @ cell, pos @ R.T)
        energies.append(e)
        forces.append(f @ R)
        weights.append(w)
    weights = np.array(weights) / np.sum(weights)
    assert energy == pytest.approx(np.dot(weights, energies), abs=ev(1e-4))
    np.testing.assert_allclose(
        force.reshape(-1, 3),
        np.einsum("r,rij->ij", weights, np.array(forces)),
        atol=ev_ang(1e-4),
    )
    np.testing.assert_allclose(extras["o3grid_pots"], energies, atol=ev(1e-4))


def test_pet_mad_is_approximately_equivariant(driver, water):
    _, _, cell, pos = water
    energy, force, _, _ = driver.compute(cell, pos)
    R = random_rotation(Random(seed=3))
    rot_energy, rot_force, _, _ = driver.compute(R @ cell, pos @ R.T)
    assert rot_energy == pytest.approx(energy, abs=ev(5e-2))
    assert unit_to_user("force", "ev/ang", np.abs(rot_force @ R - force).max()) < 0.1

"""Tests for the metatomic-core driver, with the NumPy Lennard-Jones test model."""

import numpy as np
import pytest

mta = pytest.importorskip("metatomic")
pytest.importorskip("vesin")
ase_build = pytest.importorskip("ase.build")
ase_io = pytest.importorskip("ase.io")
from ase.calculators.lj import LennardJones  # noqa: E402

from ipi.pes.metatomic_core import MetatomicCoreDriver  # noqa: E402
from ipi.utils.units import unit_to_internal, unit_to_user  # noqa: E402

LJ = dict(sigma=3.4, epsilon=0.0104, cutoff=7.0, atomic_type=18)


@pytest.fixture
def argon(tmp_path):
    atoms = ase_build.bulk("Ar", "fcc", a=5.26, cubic=True) * (2, 2, 2)
    atoms.rattle(0.05, seed=0)
    ase_io.write(tmp_path / "argon.xyz", atoms)
    (tmp_path / "model.py").write_text(
        "from metatomic.testing import LennardJonesModel\n"
        f"model = LennardJonesModel(**{LJ}, length_unit='angstrom', energy_unit='eV')\n"
    )
    return atoms, tmp_path


def _ipi(atoms):
    return (
        unit_to_internal("length", "angstrom", atoms.cell.array.T),
        unit_to_internal("length", "angstrom", atoms.positions),
    )


def test_python_model_matches_ase_lennard_jones(argon):
    atoms, path = argon
    driver = MetatomicCoreDriver(
        template=str(path / "argon.xyz"),
        model=f"{path / 'model.py'}:model",
        check_consistency=True,
    )
    energy, forces, virial, extras = driver(*_ipi(atoms))

    atoms.calc = LennardJones(
        sigma=LJ["sigma"], epsilon=LJ["epsilon"], rc=LJ["cutoff"], smooth=False
    )
    reference_virial = -atoms.get_stress(voigt=False) * atoms.get_volume()
    assert np.isclose(
        unit_to_user("energy", "electronvolt", energy), atoms.get_potential_energy()
    )
    assert np.allclose(
        unit_to_user("force", "ev/ang", forces), atoms.get_forces(), atol=1e-12
    )
    assert np.allclose(
        unit_to_user("energy", "electronvolt", virial), reference_virial, atol=1e-12
    )
    assert extras == ""

    # batched evaluation goes through a single model call
    cell, pos = _ipi(atoms)
    moved = pos.copy()
    moved[0] += 0.1
    batched = driver([cell, cell], [pos, moved])
    assert np.isclose(batched[0][0], energy)
    assert np.allclose(batched[0][1], forces)
    assert not np.isclose(batched[1][0], energy)


def test_native_plugin_matches_python_model(argon):
    atoms, path = argon
    options = {f"option_{k}": v for k, v in LJ.items()}
    try:
        native = MetatomicCoreDriver(
            template=str(path / "argon.xyz"),
            model="metatomic-lj-model",
            plugins="lj",
            option_length_unit="angstrom",
            option_energy_unit="eV",
            **options,
        )
    except ValueError as e:
        pytest.skip(f"metatomic-lj-plugin is not available: {e}")
    python = MetatomicCoreDriver(
        template=str(path / "argon.xyz"), model=f"{path / 'model.py'}:model"
    )
    (e_native, f_native, *_), (e_python, f_python, *_) = (
        driver(*_ipi(atoms)) for driver in (native, python)
    )
    assert np.isclose(e_native, e_python)
    assert np.allclose(f_native, f_python, atol=1e-12)

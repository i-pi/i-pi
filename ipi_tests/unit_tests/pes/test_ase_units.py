"""ASE bridge unit conversion checked using a real Lennard-Jones calculator."""

import json
import numpy as np
import pytest
from ase import Atoms
from ase.calculators.lj import LennardJones
from ase.io import write
from ipi.pes._ase import ASEDriver
from ipi.utils.units import unit_to_internal


@pytest.fixture
def bridge(tmp_path):
    atoms = Atoms(
        "Ar2", positions=[[0, 0, 0], [1.3, 0.2, 0.1]], cell=np.eye(3) * 8, pbc=True
    )
    path = tmp_path / "template.xyz"
    write(path, atoms)
    driver = ASEDriver(str(path))
    driver.ase_calculator = LennardJones(epsilon=0.7, sigma=1.0)
    return driver, atoms


def test_ase_energy_forces_and_virial_units(bridge):
    driver, atoms = bridge
    length = unit_to_internal("length", "angstrom", 1.0)
    energy, force, virial, _ = driver(
        atoms.cell.array.T * length, atoms.positions * length
    )
    atoms.calc = LennardJones(epsilon=0.7, sigma=1.0)
    assert energy == pytest.approx(
        unit_to_internal("energy", "electronvolt", atoms.get_potential_energy())
    )
    np.testing.assert_allclose(
        force, unit_to_internal("force", "ev/ang", atoms.get_forces()), rtol=1e-12
    )
    expected = unit_to_internal(
        "energy", "electronvolt", -atoms.get_stress(voigt=False) * atoms.get_volume()
    )
    np.testing.assert_allclose(virial.reshape(3, 3), expected, atol=1e-12)


def test_ase_tensor_stress_and_extras(bridge):
    driver, atoms = bridge
    props = {
        "energy": 2.0,
        "forces": np.ones((2, 3)),
        "stress": np.diag([1.0, 2.0, 3.0]),
        "dipole": np.array([1.0, 2.0, 3.0]),
        "uncertainty": 0.2,
    }
    _, _, vir, extras = driver.post_process(props, atoms)
    np.testing.assert_allclose(
        vir,
        unit_to_internal(
            "energy", "electronvolt", -props["stress"] * atoms.get_volume()
        ),
    )
    assert json.loads(extras) == {"dipole": [1.0, 2.0, 3.0], "uncertainty": 0.2}
    props["stress"] = np.ones(4)
    with pytest.raises(ValueError, match="shape"):
        driver.post_process(props, atoms)
    props["stress"] = np.eye(3)
    props["stress"][0, 1] = 1
    with pytest.raises(ValueError, match="symmetric"):
        driver.post_process(props, atoms)


@pytest.mark.xfail(
    strict=True,
    raises=ValueError,
    reason="ASE bridge rejects its own zero-stress fallback when stress is disabled",
)
def test_ase_optional_stress_is_zero(bridge):
    driver, atoms = bridge
    driver.capabilities = ["energy", "forces"]
    _, _, vir, _ = driver.post_process(
        {"energy": 1.0, "forces": np.zeros((2, 3))}, atoms
    )
    np.testing.assert_array_equal(vir, 0)

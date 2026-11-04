"""Tests for the boundary of the generic MACE driver."""

from pathlib import Path

import numpy as np
import pytest
import torch
from ase.io import read

pytest.importorskip("mace")

from ipi.pes._mace import BatchedMACE, proper_dipole
from ipi.pes.tools import ModelResults
from mace.calculators import MACECalculator


@pytest.fixture
def mace_mp_model_path():
    """Return the model downloaded by the MACE examples workflow."""

    model_path = (
        Path(__file__).parents[3]
        / "examples"
        / "clients"
        / "mace"
        / "single"
        / "mace.model"
    )
    if not model_path.is_file():
        pytest.skip(
            "requires the MACE-MP model downloaded by "
            "examples/clients/mace/single/getmodel.sh"
        )
    return model_path


def test_batched_mace_matches_mace_calculator(mace_mp_model_path):
    """Match MACECalculator when i-PI derives forces and stress externally."""

    atoms = read(mace_mp_model_path.parent / "init.xyz")

    reference_atoms = atoms.copy()
    reference_atoms.calc = MACECalculator(
        model_paths=str(mace_mp_model_path), device="cpu"
    )
    reference = {
        "energy": reference_atoms.get_potential_energy(),
        "forces": reference_atoms.get_forces(),
        "stress": reference_atoms.get_stress(voigt=False),
    }

    result = BatchedMACE(
        model_paths=str(mace_mp_model_path), device="cpu"
    ).compute_batched([atoms])[0]

    for property_name, expected in reference.items():
        np.testing.assert_allclose(
            result[property_name], expected, rtol=1e-10, atol=1e-12
        )


def test_mace_warns_and_skips_unknown_model_outputs():
    """Unknown MACE outputs do not abort an evaluation."""

    calculator = object.__new__(BatchedMACE)
    calculator.ase_like_properties = {"energy": ()}
    results = ModelResults(calculator.ase_like_properties)

    with pytest.warns(UserWarning, match="Unknown model properties"):
        calculator._store_registered_results(
            results,
            [1],
            {
                "energy": np.asarray([1.0]),
                "unknown": np.zeros((1, 3)),
            },
        )

    assert len(results) == 1
    assert results[0] == {"energy": 1.0}


def test_plain_mace_accepts_generic_dielectric_response_instruction():
    """The base MACE driver owns generic BEC/piezo response support."""

    assert "compute_BEC" in BatchedMACE.supported_instruction_keys
    assert hasattr(BatchedMACE, "add_dielectric_response")
    assert hasattr(BatchedMACE, "compute_dmu_dR_deta")
    assert hasattr(BatchedMACE, "_response_dipole")


def test_proper_dipole_is_symmetric_strain_corrected():
    """The response dipole uses the proper-dipole strain correction."""

    dipole = np.asarray([[1.0, 2.0, 3.0]])
    strain = np.zeros((1, 3, 3))
    strain[0, 0, 1] = 0.2
    strain[0, 1, 0] = 0.4
    corrected = proper_dipole(torch.as_tensor(dipole), torch.as_tensor(strain)).numpy()
    np.testing.assert_allclose(corrected[0], [0.4, 1.7, 3.0])

"""Tests for the boundary of the generic MACE driver."""

from pathlib import Path

import numpy as np
import pytest
from ase.io import read

pytest.importorskip("mace")

from ipi.pes._mace import BatchedMACE, ase_like_properties
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


def test_plain_mace_has_no_electrical_response_implementation():
    """The base calculator remains independent of response extensions."""

    assert "BEC" not in ase_like_properties
    assert "piezoelectric" not in ase_like_properties
    assert "compute_BEC" not in BatchedMACE.supported_instruction_keys
    assert not hasattr(BatchedMACE, "add_dielectric_response")
    assert not hasattr(BatchedMACE, "compute_dmu_dR_deta")
    assert not hasattr(BatchedMACE, "_response_dipole")


def test_plain_mace_rejects_extension_instructions():
    """Extension settings are not silently ignored by plain MACE."""

    with pytest.raises(ValueError, match="Unsupported mace instruction"):
        BatchedMACE(instructions={"compute_BEC": True})

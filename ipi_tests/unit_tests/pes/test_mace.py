"""Tests for the boundary of the generic MACE driver."""

import pytest

pytest.importorskip("mace")

from ipi.pes._mace import BatchedMACE, ase_like_properties


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

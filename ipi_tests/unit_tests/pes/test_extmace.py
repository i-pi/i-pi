"""Tests for the FFDielectric additions to the MACE driver."""

import json

import pytest

pytest.importorskip("mace")
torch = pytest.importorskip("torch")

from ipi.pes._mace import BatchedMACE, MACE_driver
from ipi.pes.extmace import (
    Extended_MACE_driver,
    ExtendedMACECalculator,
    proper_dipole,
)


def test_extmace_acknowledges_applied_electric_field(monkeypatch):
    """An Efield request is acknowledged in the returned extras JSON."""
    result = (1.0, "forces", "virial", '{"dipole": [0.0, 0.0, 0.0]}')
    monkeypatch.setattr(MACE_driver, "post_process", lambda *args: result)

    driver = object.__new__(Extended_MACE_driver)
    driver.extra = {"Efield": [0.0, 0.0, 0.1]}

    _, _, _, extras = driver.post_process({}, None)

    assert json.loads(extras)["applied_fields"] == ["electric_field"]


def test_proper_dipole_normalizes_singleton_model_dimension():
    """EnergyDipoleMACE outputs with a singleton wrapper remain compatible."""

    mu = torch.tensor([[[1.0, 2.0, 3.0]], [[4.0, 5.0, 6.0]]])
    strain = torch.zeros((2, 3, 3))
    strain[:, 0, 0] = 0.25

    corrected = proper_dipole(mu, strain)

    expected = mu.squeeze(-2).clone()
    expected[:, 0] *= 0.75
    assert corrected.shape == (2, 3)
    assert torch.allclose(corrected, expected)


def test_proper_dipole_rejects_non_singleton_model_dimension():
    """A malformed per-atom-like dipole is not silently flattened."""

    mu = torch.zeros((1, 4, 3))
    strain = torch.zeros((1, 3, 3))

    with pytest.raises(ValueError, match="Only extra singleton dimensions"):
        proper_dipole(mu, strain)


def test_extmace_skips_dipole_coupling_for_zero_field(monkeypatch):
    """A transmitted but exactly zero field does not require a model dipole."""

    calculator = object.__new__(ExtendedMACECalculator)
    calculator.extras = {"Efield": [0.0, 0.0, 0.0]}
    calculator.compute_bec_response = False
    data = {"energy": torch.tensor([1.0])}

    monkeypatch.setattr(
        calculator,
        "_response_dipole",
        lambda _: pytest.fail("zero field should not request the dipole"),
    )
    monkeypatch.setattr(
        calculator,
        "get_forces_stress",
        lambda output, batch, training: output,
    )

    result = calculator.augment_output(
        data=data,
        batch={},
        training=False,
    )

    assert result is data
    assert torch.equal(result["energy"], torch.tensor([1.0]))


def test_plain_mace_has_no_electrical_response_implementation():
    """The base calculator remains usable without the extmace module."""

    assert not hasattr(BatchedMACE, "add_dielectric_response")
    assert not hasattr(BatchedMACE, "compute_dmu_dR_deta")
    assert not hasattr(BatchedMACE, "_response_dipole")

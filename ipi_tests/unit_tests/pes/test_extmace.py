"""Tests for the FFDielectric additions to the MACE driver."""

import argparse
import json
from types import SimpleNamespace

import numpy as np
import pytest

pytest.importorskip("mace")
torch = pytest.importorskip("torch")

from ipi.pes._mace import BatchedMACE, MACE_driver
from ipi.pes.extmace import (
    Extended_MACE_driver,
    ExtendedMACECalculator,
    _field_extra_from_structure,
    add_extmace_cli_arguments,
    check_electric_boundary_condition,
    coerce_epsilon_infinity,
    evaluate_extmace_structures,
    proper_dipole,
)
from ipi.utils.units import unit_to_internal, unit_to_user


def test_extmace_acknowledges_applied_electric_field(monkeypatch):
    """An electric_field request is acknowledged in the returned extras JSON."""
    result = (1.0, "forces", "virial", '{"dipole": [0.0, 0.0, 0.0]}')
    monkeypatch.setattr(MACE_driver, "post_process", lambda *args: result)

    driver = object.__new__(Extended_MACE_driver)
    driver.extra = {"electric_field": [0.0, 0.0, 0.1]}

    _, _, _, extras = driver.post_process({}, None)

    assert json.loads(extras)["applied_fields"] == ["electric_field"]


def test_extmace_acknowledges_applied_electric_displacement(monkeypatch):
    """A Dfield request is acknowledged in the returned extras JSON."""

    result = (1.0, "forces", "virial", '{"dipole": [0.0, 0.0, 0.0]}')
    monkeypatch.setattr(MACE_driver, "post_process", lambda *args: result)

    driver = object.__new__(Extended_MACE_driver)
    driver.extra = {"Dfield": [0.0, 0.0, 0.1]}

    _, _, _, extras = driver.post_process({}, None)

    assert json.loads(extras)["applied_fields"] == ["electric_displacement"]


def test_extmace_rejects_mixed_electric_boundary_conditions():
    """Fixed E and fixed D cannot be applied to the same structure."""

    with pytest.raises(ValueError, match="cannot mix"):
        check_electric_boundary_condition(
            {"electric_field": [0, 0, 0], "Dfield": [0, 0, 0]}
        )


def test_extmace_expands_voigt_epsilon_infinity():
    """JSON-compatible Voigt dielectric data is converted to Cartesian form."""

    epsilon = coerce_epsilon_infinity([1.0, 2.0, 3.0, 0.4, 0.5, 0.6])

    assert np.array_equal(
        epsilon,
        np.array([[1.0, 0.6, 0.5], [0.6, 2.0, 0.4], [0.5, 0.4, 3.0]]),
    )


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
    calculator.extras = {"electric_field": [0.0, 0.0, 0.0]}
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


def _extmace_argument_parser():
    parser = argparse.ArgumentParser()
    add_extmace_cli_arguments(parser)
    return parser


def test_extmace_cli_reads_electric_field_from_atoms_info():
    """The selected extxyz info value is converted to i-PI atomic units."""

    args = _extmace_argument_parser().parse_args(
        ["--electric-field-key", "applied_E", "--field-units", "v/ang"]
    )
    atoms = SimpleNamespace(info={"applied_E": [0.1, -0.2, 0.3]})

    extras = _field_extra_from_structure(atoms, 4, args)

    expected = unit_to_internal("electric-field", "v/ang", np.array([0.1, -0.2, 0.3]))
    np.testing.assert_allclose(extras["electric_field"], expected)


def test_extmace_cli_reads_electric_displacement_from_atoms_info():
    """The displacement-key option creates a fixed-D calculator request."""

    args = _extmace_argument_parser().parse_args(
        ["--electric-displacement-key", "applied_D"]
    )
    atoms = SimpleNamespace(info={"applied_D": [0.0, 0.0, 0.1]})

    extras = _field_extra_from_structure(atoms, 0, args)

    assert set(extras) == {"Dfield"}
    np.testing.assert_allclose(
        extras["Dfield"],
        unit_to_internal("electric-field", "v/ang", np.array([0.0, 0.0, 0.1])),
    )


def test_extmace_cli_rejects_missing_frame_field():
    """Every frame must define the selected Atoms.info key."""

    args = _extmace_argument_parser().parse_args(["--electric-field-key", "applied_E"])

    with pytest.raises(ValueError, match="structure 2.*applied_E"):
        _field_extra_from_structure(SimpleNamespace(info={}), 2, args)


def test_extmace_cli_evaluates_different_fields_separately():
    """Per-frame fields are not incorrectly shared by a batched evaluation."""

    class RecordingCalculator:
        def __init__(self):
            self.extras = {"original": True}
            self.calls = []

        def compute_batched(self, structures):
            self.calls.append((list(structures), self.extras.copy()))
            return [{"energy": float(len(self.calls))}]

    args = _extmace_argument_parser().parse_args(
        ["--electric-field-key", "E", "--field-units", "atomic_unit"]
    )
    structures = [
        SimpleNamespace(info={"E": [0.0, 0.0, 0.1]}),
        SimpleNamespace(info={"E": [0.0, 0.0, 0.2]}),
    ]
    calculator = RecordingCalculator()

    results = evaluate_extmace_structures(calculator, structures, args)

    assert [result["energy"] for result in results] == [1.0, 2.0]
    assert len(calculator.calls) == 2
    assert calculator.calls[0][1] == {"electric_field": [0.0, 0.0, 0.1]}
    assert calculator.calls[1][1] == {"electric_field": [0.0, 0.0, 0.2]}
    assert calculator.extras == {"original": True}


def test_extmace_constant_d_energy_is_differentiable(monkeypatch):
    """Client-side fixed D adds a differentiable energy functional."""

    calculator = object.__new__(ExtendedMACECalculator)
    calculator.extras = {"Dfield": [0.0, 0.0, 0.1]}
    calculator.default_epsilon_infinity = np.eye(3)
    calculator.compute_bec_response = False
    dipole = torch.tensor([[0.0, 0.0, 0.2]], requires_grad=True)
    cell = torch.eye(3).unsqueeze(0).requires_grad_()
    data = {"energy": torch.zeros(1)}

    monkeypatch.setattr(calculator, "_response_dipole", lambda _: dipole)
    monkeypatch.setattr(
        calculator, "get_forces_stress", lambda output, batch, training: output
    )

    result = calculator.augment_output(data=data, batch={"cell": cell}, training=False)
    result["energy"].sum().backward()

    assert dipole.grad is not None
    assert torch.isfinite(dipole.grad).all()


def test_extmace_constant_d_uses_mace_field_units():
    """The fixed-D gradient uses the same force convention as fixed E."""

    calculator = object.__new__(ExtendedMACECalculator)
    dfield_atomic = np.array([0.0, 0.0, 0.1])
    calculator.extras = {"Dfield": dfield_atomic.tolist()}
    calculator.default_epsilon_infinity = np.eye(3)

    dipole = torch.zeros((1, 3), dtype=torch.float64, requires_grad=True)
    cell = 10.0 * torch.eye(3, dtype=torch.float64).unsqueeze(0)
    reference = torch.zeros(1, dtype=torch.float64)
    data = {"energy": reference}

    transmitted_dfield = calculator._electric_displacement(reference)
    np.testing.assert_allclose(
        transmitted_dfield.detach().numpy(),
        unit_to_user("electric-field", "v/ang", dfield_atomic),
    )

    energy = calculator._constant_d_energy(dipole, data, {"cell": cell})
    energy.sum().backward()

    force_atomic = -dipole.grad.detach().numpy() * unit_to_internal("force", "ev/ang")
    np.testing.assert_allclose(force_atomic[0], dfield_atomic, rtol=0.0, atol=1e-14)


def test_plain_mace_has_no_electrical_response_implementation():
    """The base calculator remains usable without the extmace module."""

    assert not hasattr(BatchedMACE, "add_dielectric_response")
    assert not hasattr(BatchedMACE, "compute_dmu_dR_deta")
    assert not hasattr(BatchedMACE, "_response_dipole")

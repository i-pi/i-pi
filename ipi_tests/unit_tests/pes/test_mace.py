"""Tests for the boundary of the generic MACE driver."""

import json
import sys

import numpy as np
import pytest
from ase import Atoms
from ase.io import read, write

pytest.importorskip("mace")

from ipi.pes import _mace
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


def test_normalize_ase_like_properties_accepts_json_compatible_shapes():
    """Additional MACE outputs accept shapes that can be written in JSON."""

    normalized = BatchedMACE._normalize_ase_like_properties(
        {
            "scalar": [],
            "tensor": [3, 3],
            "per_atom": ["natoms", 3],
        }
    )

    assert normalized == {
        "scalar": (),
        "tensor": (3, 3),
        "per_atom": ("natoms", 3),
    }


@pytest.mark.parametrize(
    ("properties", "exception", "message"),
    [
        ([], TypeError, "must be a dictionary"),
        ({1: []}, TypeError, "must be strings"),
        ({"bad": "scalar"}, TypeError, "must be a list or tuple"),
        ({"bad": [-1]}, ValueError, "non-negative integers"),
        ({"bad": [3, "natoms"]}, ValueError, "first dimension"),
    ],
)
def test_normalize_ase_like_properties_rejects_invalid_shapes(
    properties, exception, message
):
    """Invalid user-defined output shapes fail before a MACE run starts."""

    with pytest.raises(exception, match=message):
        BatchedMACE._normalize_ase_like_properties(properties)


def test_mace_output_summary_describes_registered_and_unknown_outputs(capsys):
    """The one-time diagnostic identifies output registration choices."""

    calculator = object.__new__(BatchedMACE)
    calculator.ase_like_properties = {
        "energy": (),
        "forces": ("natoms", 3),
    }
    calculator.device = "cpu"

    calculator._print_output_summary(
        {
            "energy": np.asarray([1.0]),
            "forces": np.zeros((2, 3)),
            "unknown": np.zeros((1, 2)),
            "absent": None,
        },
        {
            "energy": np.asarray([1.0]),
            "forces": np.zeros((2, 3)),
            "unknown": np.zeros((1, 2)),
        },
    )

    output = capsys.readouterr().out
    assert "Produced per-atom properties (ASE arrays): forces" in output
    assert "Produced per-structure properties (ASE info): energy" in output
    assert "Produced unregistered model outputs: unknown" in output
    assert "Register an output under 'ase_like_properties'" in output


def test_mace_cli_loads_custom_ase_like_properties(tmp_path, monkeypatch):
    """The standalone CLI forwards JSON output registrations to its calculator."""

    input_path = tmp_path / "input.extxyz"
    output_path = tmp_path / "output.extxyz"
    settings_path = tmp_path / "properties.json"
    write(input_path, Atoms("H", positions=[[0.0, 0.0, 0.0]]))
    settings_path.write_text(json.dumps({"custom": [3]}))

    received = {}

    class FakeCalculator:
        def __init__(self, model_paths, device, **kwargs):
            received.update(model_paths=model_paths, device=device, kwargs=kwargs)
            self.ase_like_properties = {"energy": (), **kwargs["ase_like_properties"]}

        def compute_batched(self, structures):
            return [{"energy": 1.0, "custom": np.arange(3.0)} for _ in structures]

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "mace-cli",
            "--model",
            "model.model",
            "--input_structures",
            str(input_path),
            "--output_structures",
            str(output_path),
            "--ase_like_properties",
            str(settings_path),
        ],
    )

    _mace.run_cli(calculator_class=FakeCalculator, calculator_name="FakeCalculator")

    assert received == {
        "model_paths": "model.model",
        "device": "cpu",
        "kwargs": {"ase_like_properties": {"custom": [3]}},
    }
    result = read(output_path)
    np.testing.assert_array_equal(result.info["MACE_custom"], np.arange(3.0))

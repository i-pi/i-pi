import json

import numpy as np
import pytest

from ipi.pes.tools import Instructions, ModelResults, convert, process_input
from ipi.utils.units import unit_to_internal


def test_model_results_stores_per_structure_and_per_atom_outputs():
    results = ModelResults({"energy": (), "forces": ("natoms", 3)})
    results.store(
        [1, 2],
        {
            "energy": np.array([1.5, 2.5]),
            "forces": np.arange(9).reshape(3, 3),
        },
    )

    assert len(results) == 2
    assert results.natoms() == [1, 2]
    np.testing.assert_allclose(results[0]["energy"], 1.5)
    np.testing.assert_allclose(results[1]["forces"], np.arange(3, 9).reshape(2, 3))


def test_unknown_mace_outputs_suggest_shapes_and_ignore_entries():
    results = ModelResults({})

    with pytest.raises(ValueError) as error:
        results.store(
            [2, 3],
            {
                "atomwise": np.zeros((5, 4)),
                "structure": np.zeros((2, 6)),
                "unclassified": np.zeros((4, 2)),
            },
        )

    message = str(error.value)
    assert (
        "'atomwise': raw shape (5, 4); inferred per-atom shape [\"natoms\", 4]"
        in message
    )
    assert "'structure': raw shape (2, 6); inferred per-structure shape [6]" in message
    assert (
        "'unclassified': raw shape (4, 2); no atom/batch leading dimension" in message
    )
    assert '"ase_like_properties": {' in message
    assert '"atomwise": [' in message
    assert '"ignore": [' in message


def test_unknown_mace_output_reports_ambiguous_shape():
    results = ModelResults({})

    with pytest.raises(ValueError) as error:
        results.store([1, 1], {"ambiguous": np.zeros((2, 3))})

    message = str(error.value)
    assert "'ambiguous': raw shape (2, 3); ambiguous shape" in message
    assert '["natoms", 3] or [3]' in message
    assert "Ambiguous properties are omitted" in message
    assert '"ase_like_properties": {}' in message


def test_model_results_mean_rejects_empty_and_incompatible_results():
    with pytest.raises(ValueError, match="Cannot compute mean"):
        ModelResults.mean([])

    first = ModelResults({"energy": ()})
    first.store([1], {"energy": np.array([1.0])})
    second = ModelResults({"energy": ()})
    second.store([2], {"energy": np.array([2.0])})

    with pytest.raises(ValueError, match="same natoms and shapes"):
        ModelResults.mean([first, second])


def test_model_results_mean_averages_models():
    first = ModelResults({"energy": (), "forces": ("natoms", 1)})
    second = ModelResults({"energy": (), "forces": ("natoms", 1)})
    first.store([2], {"energy": np.array([2.0]), "forces": np.array([[1.0], [3.0]])})
    second.store([2], {"energy": np.array([4.0]), "forces": np.array([[3.0], [5.0]])})

    mean = ModelResults.mean([first, second])

    np.testing.assert_allclose(mean[0]["energy"], 3.0)
    np.testing.assert_allclose(mean[0]["forces"], [[2.0], [4.0]])


class _LengthInstructions(Instructions):
    dimensions = {"cutoff": "length"}
    units = {"length": "atomic_unit"}


def test_instructions_convert_units_from_dict_and_json_file(tmp_path):
    settings = {"cutoff": [1.0, 2.0], "cutoff_unit": "angstrom"}
    parsed = _LengthInstructions(settings)
    expected = unit_to_internal("length", "angstrom", np.array([1.0, 2.0]))

    np.testing.assert_allclose(parsed.instructions["cutoff"], expected)
    assert "cutoff_unit" not in parsed.instructions

    source = tmp_path / "instructions.json"
    source.write_text(json.dumps({"cutoff": 3.0, "cutoff_unit": None}))
    parsed_file = _LengthInstructions(str(source))

    assert parsed_file.instructions == {"cutoff": 3.0}


def test_instructions_reject_invalid_source():
    with pytest.raises(ValueError, match=r"str.*dict"):
        _LengthInstructions([])


@pytest.mark.parametrize(
    ("value", "expected"),
    [(1.5, 1.5), (2, 2), ([1, 2], np.array([1, 2]))],
)
def test_process_input(value, expected):
    result = process_input(value)
    np.testing.assert_equal(result, expected)


def test_process_input_rejects_unsupported_value():
    with pytest.raises(TypeError, match="float or a list"):
        process_input("not a number")


def test_convert_with_and_without_a_unit_family():
    assert convert(7.6) == 7.6
    assert convert(7.6, "length", "angstrom", "angstrom") == 7.6

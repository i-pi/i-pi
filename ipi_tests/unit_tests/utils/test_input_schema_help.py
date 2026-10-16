"""Generated input references must expose dimensions, defaults, and nested options."""

import xml.etree.ElementTree as ET

import numpy as np
import pytest

from ipi.utils.inputvalue import (
    Input,
    InputArray,
    InputAttribute,
    InputDictionary,
    InputValue,
)
from ipi_tests.unit_tests.utils.test_inputvalue_roundtrip import parse_value


class ExampleInput(Input):
    default_label = "EXAMPLE"
    default_help = "Example parameters"
    attribs = {
        "mode": (
            InputAttribute,
            {"dtype": str, "default": "nve", "options": ["nve", "nvt"]},
        )
    }
    fields = {
        "distance": (
            InputValue,
            {
                "dtype": float,
                "dimension": "length",
                "default": 2.0,
                "help": "Distance value",
            },
        )
    }
    dynamic = {
        "sample": (
            InputArray,
            {"dtype": float, "default": np.zeros(0), "help": "Repeated measurements"},
        )
    }


@pytest.mark.parametrize("format", ["xml", "latex", "rst"])
def test_input_reference_includes_nested_fields_defaults_and_choices(format):
    schema = ExampleInput()
    text = getattr(schema, "help_" + format)(name="example")
    for required in [
        "distance",
        "length",
        "mode",
        "nve",
        "nvt",
        "sample",
        "Repeated measurements",
        "Distance value",
    ]:
        assert required in text
    if format == "xml":
        root = ET.fromstring(text)
        assert root.tag == "example"
        assert root.find("distance/dimension").text.strip() == "length"
        assert root.find("distance/dtype").text.strip() == "float"
    elif format == "latex":
        assert text.count(r"\begin{ipifield}") == text.count(r"\end{ipifield}")
        assert r"\begin{document}" in text and r"\end{document}" in text
    else:
        assert "ATTRIBUTES" in text and "FIELDS" in text


@pytest.mark.parametrize("format", ["xml", "latex", "rst"])
def test_reference_depth_limit_omits_grandchildren(format):
    schema = ExampleInput()
    render = getattr(schema, "help_" + format)
    assert render(name="example", level=2, stop_level=1) == ""
    text = render(name="example", stop_level=0)
    assert "Repeated measurements" not in text
    assert "Distance value" not in text


@pytest.mark.parametrize(
    "dtype,label",
    [
        (bool, "boolean"),
        (int, "integer"),
        (float, "float"),
        (str, "string"),
        (dict, "dictionary"),
        (tuple, "tuple"),
    ],
)
def test_schema_type_names(dtype, label):
    assert Input().type_print(dtype) == label


@pytest.mark.parametrize("field", ["default", "dtype", "dimension"])
def test_dictionary_rejects_mismatched_metadata_lengths(field):
    kwargs = dict(
        options=["a", "b"],
        default=[1.0, 2.0],
        dtype=[float, float],
        dimension=["length", "length"],
    )
    kwargs[field] = kwargs[field][:1]
    with pytest.raises(ValueError, match="length"):
        InputDictionary(**kwargs)


def test_homogeneous_dictionary_xml_roundtrip():
    original = InputDictionary(options=["left", "right"], default=0.0, dtype=float)
    original.store({"left": 1.25, "right": -2.5})
    restored = parse_value(
        InputDictionary(options=["left", "right"], default=0.0, dtype=float),
        original.write("pair"),
    )
    assert restored.fetch() == {"left": 1.25, "right": -2.5}

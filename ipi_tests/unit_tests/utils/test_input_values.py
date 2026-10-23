"""Input values: XML round trips, validation and generated documentation."""

import xml.etree.ElementTree as ET

import numpy as np
import pytest

from ipi.inputs.forcefields import InputFFSocket
from ipi.inputs.interface import InputInterfaceSocket
from ipi.utils.inputvalue import (
    Input,
    InputArray,
    InputAttribute,
    InputDictionary,
    InputValue,
)
from ipi.utils.io.inputs.io_xml import xml_parse_string


def parse_value(value, xml):
    """Parse a single XML element into the supplied input object."""
    value.parse(xml_parse_string(xml).fields[0][1])
    return value


@pytest.mark.parametrize(
    "dtype,value",
    [(bool, True), (bool, False), (int, -7), (float, 1.23456789), (str, "hydrogen")],
)
def test_scalar_xml_roundtrip(dtype, value):
    """The text representation must preserve each supported scalar type."""
    original = InputValue(dtype=dtype)
    original.store(value)
    restored = parse_value(InputValue(dtype=dtype), original.write("value"))
    assert restored.fetch() == value


@pytest.mark.parametrize("shape", [(0,), (1,), (2, 3), (2, 2, 3)])
@pytest.mark.parametrize("units", ["atomic_unit", "angstrom"])
def test_array_xml_roundtrip(shape, units):
    """Unit conversion and flattening must not alter the restart array."""
    data = np.arange(np.prod(shape), dtype=float).reshape(shape) / 7
    original = InputArray(dtype=float, dimension="length")
    original.store(data, units=units)
    restored = parse_value(
        InputArray(dtype=float, dimension="length"), original.write("q")
    )
    result = restored.fetch()
    assert result.shape == shape
    np.testing.assert_allclose(result, data, rtol=1e-7, atol=1e-12)
    result[:] = -1
    np.testing.assert_allclose(restored.fetch(), data, rtol=1e-7, atol=1e-12)
    data[:] = -2
    assert not np.any(original.fetch() == -2)


def test_file_array_is_self_contained_after_checkpoint(tmp_path):
    """A checkpoint must embed the numbers rather than depend on the initial file."""
    path = tmp_path / "positions.txt"
    expected = np.arange(6, dtype=float).reshape(2, 3)
    np.savetxt(path, expected)
    loaded = parse_value(
        InputArray(dtype=float), f'<q mode="file" shape="(2,3)">{path}</q>'
    )
    stored = InputArray(dtype=float)
    stored.store(loaded.fetch())
    path.unlink()
    restored = parse_value(InputArray(dtype=float), stored.write("q"))
    np.testing.assert_array_equal(restored.fetch(), expected)


@pytest.mark.parametrize(
    "xml,error",
    [
        ('<q shape="(2,3)">[1,2]</q>', ValueError),
        ('<q mode="invalid">[1,2]</q>', ValueError),
        ('<q unsupported="yes">[1,2]</q>', NameError),
    ],
)
def test_invalid_array_input_is_rejected(xml, error):
    """Malformed dimensions, modes, and attributes must fail before simulation."""
    with pytest.raises(error):
        parse_value(InputArray(dtype=float), xml).fetch()


def test_invalid_boolean_is_not_silently_true():
    """Python's truthiness of nonempty strings must not accept misspelled booleans."""
    with pytest.raises(ValueError, match="bool"):
        parse_value(InputValue(dtype=bool), "<enabled>flase</enabled>")


def test_homogeneous_dictionary_xml_roundtrip():
    original = InputDictionary(options=["left", "right"], default=0.0, dtype=float)
    original.store({"left": 1.25, "right": -2.5})
    restored = parse_value(
        InputDictionary(options=["left", "right"], default=0.0, dtype=float),
        original.write("pair"),
    )
    assert restored.fetch() == {"left": 1.25, "right": -2.5}


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


@pytest.mark.parametrize("schema", [InputFFSocket, InputInterfaceSocket])
@pytest.mark.parametrize(
    "field,value",
    [
        ("port", 0),
        ("port", 65536),
        ("slots", 0),
        ("slots", 6),
        ("latency", -0.01),
        ("timeout", -1),
    ],
)
def test_transport_input_rejects_invalid_ranges(schema, field, value):
    obj = parse_value(
        schema(),
        (
            f'<socket name="test"><{field}>{value}</{field}></socket>'
            if schema is InputFFSocket
            else f"<socket><{field}>{value}</{field}></socket>"
        ),
    )
    with pytest.raises(ValueError):
        obj.check()


@pytest.mark.parametrize("transport", ["unix", "inet", "shm"])
def test_force_socket_settings_survive_xml_roundtrip(transport):
    obj = parse_value(
        InputFFSocket(),
        f'<ffsocket name="test" mode="{transport}" pbc="False"><port>31415</port><address>example</address><latency>.002</latency><timeout>7</timeout><slots>3</slots></ffsocket>',
    )
    ff = obj.fetch()
    stored = InputFFSocket()
    stored.store(ff)
    restored = parse_value(InputFFSocket(), stored.write("ffsocket")).fetch()
    assert restored.name == "test"
    assert restored.dopbc is False
    for attr in ["mode", "port", "address", "timeout", "slots"]:
        assert getattr(restored.socket, attr) == getattr(ff.socket, attr)
    assert restored.latency == pytest.approx(0.002)


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
def test_generated_input_reference(format):
    schema = ExampleInput()
    render = getattr(schema, "help_" + format)
    text = render(name="example")
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

    assert render(name="example", level=2, stop_level=1) == ""
    limited = render(name="example", stop_level=0)
    assert "Repeated measurements" not in limited
    assert "Distance value" not in limited

"""Input and restart values preserve types, shapes, units, and stored data."""

import numpy as np
import pytest

from ipi.utils.inputvalue import InputArray, InputValue
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
def test_array_xml_roundtrip_preserves_shape_units_and_independent_storage(
    shape, units
):
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

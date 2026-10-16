"""Reject invalid transport settings before starting a force service."""

import pytest
from ipi.inputs.forcefields import InputFFSocket
from ipi.inputs.interface import InputInterfaceSocket
from ipi_tests.unit_tests.utils.test_inputvalue_roundtrip import parse_value


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

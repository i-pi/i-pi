"""Tests for XML template builders in ipi.scripting.templates."""

import pytest

from ipi.scripting import forcefield_xml, motion_nvt_xml
from ipi.utils.io.inputs.io_xml import xml_parse_string


def _parse(xml_fragment):
    """Wrap a bare XML fragment in a root so xml_parse_string accepts it."""
    return xml_parse_string(f"<root>{xml_fragment}</root>")


def test_forcefield_xml_direct_happy_path():
    xml = forcefield_xml(name="foo", mode="direct", pes="dummy")
    assert "<ffdirect" in xml
    assert "name='foo'" in xml
    assert "<pes>dummy</pes>" in xml
    _parse(xml)  # must be valid XML


def test_forcefield_xml_direct_with_dict_parameters():
    xml = forcefield_xml(
        name="foo", mode="direct", pes="dummy", parameters={"a": 1, "b": "x"}
    )
    assert "<parameters>" in xml
    _parse(xml)


def test_forcefield_xml_direct_requires_pes():
    with pytest.raises(ValueError, match="pes"):
        forcefield_xml(name="foo", mode="direct")


def test_forcefield_xml_direct_rejects_unknown_pes():
    with pytest.raises(ValueError, match="Invalid value"):
        forcefield_xml(name="foo", mode="direct", pes="not_a_real_pes_name")


def test_forcefield_xml_unix_has_address_no_port():
    xml = forcefield_xml(name="s", mode="unix", address="/tmp/sock")
    assert "<ffsocket" in xml
    assert "<address>/tmp/sock</address>" in xml
    assert "<port>" not in xml
    _parse(xml)


def test_forcefield_xml_inet_requires_port():
    with pytest.raises(ValueError):
        forcefield_xml(name="s", mode="inet", address="localhost")


def test_forcefield_xml_inet_includes_port():
    xml = forcefield_xml(name="s", mode="inet", address="localhost", port=31415)
    assert "<ffsocket" in xml
    assert "<port>31415</port>" in xml
    _parse(xml)


def test_forcefield_xml_rejects_unknown_mode():
    with pytest.raises(ValueError, match="Invalid forcefield mode"):
        forcefield_xml(name="s", mode="bogus")


def test_motion_nvt_xml_defaults_to_svr():
    xml = motion_nvt_xml(timestep=0.5)
    assert 'mode="dynamics"' in xml
    assert 'mode="nvt"' in xml
    assert "mode='svr'" in xml
    assert "0.5" in xml  # timestep
    _parse(xml)


def test_motion_nvt_xml_path_integrals_uses_pile_g():
    xml = motion_nvt_xml(timestep=0.5, path_integrals=True)
    assert "mode='pile_g'" in xml
    assert "pile_lambda" in xml
    _parse(xml)

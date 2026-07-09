"""Tests for XML template builders in ipi.scripting.templates."""

import numpy as np
import pytest

from ipi.scripting import forcefield_xml, motion_nvt_xml, simulation_xml
from ipi.utils.io.inputs.io_xml import xml_parse_string


def _parse(xml_fragment):
    """Wrap a bare XML fragment in a root so xml_parse_string accepts it."""
    return xml_parse_string(f"<root>{xml_fragment}</root>")


def _atoms(symbols="H2"):
    """A minimal ASE structure in a cubic 5 Ang cell."""
    ase = pytest.importorskip("ase")
    return ase.Atoms(
        symbols,
        positions=[[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]],
        cell=5.0 * np.eye(3),
        pbc=True,
    )


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


def _basic_blocks():
    """A valid forcefield and motion block to feed simulation_xml."""
    return (
        forcefield_xml(name="myff", mode="direct", pes="dummy"),
        motion_nvt_xml(timestep=0.5),
    )


def test_simulation_xml_happy_path():
    ff, motion = _basic_blocks()
    xml = simulation_xml(_atoms(), ff, motion)

    parsed = xml_parse_string(xml)  # must be a valid, single <simulation> block
    assert parsed.fields[0][0] == "simulation"
    # verbosity / safe_stride attributes and the forcefield reference are wired in
    assert "verbosity='quiet'" in xml
    assert "safe_stride='20'" in xml
    assert "<ffdirect" in xml
    assert "forcefield='myff'" in xml
    assert "<beads" in xml and "<cell" in xml


def test_simulation_xml_without_temperature_has_no_ensemble():
    ff, motion = _basic_blocks()
    xml = simulation_xml(_atoms(), ff, motion)
    # NVE: no thermal initialization nor ensemble temperature
    assert "initialize" not in xml
    assert "ensemble" not in xml
    xml_parse_string(xml)


def test_simulation_xml_with_temperature_sets_ensemble():
    ff, motion = _basic_blocks()
    xml = simulation_xml(_atoms(), ff, motion, temperature=300.0)
    assert "<initialize" in xml
    assert "mode='thermal'" in xml
    assert "<ensemble>" in xml
    assert "<temperature units='ase'> 300.0 </temperature>" in xml
    xml_parse_string(xml)


def test_simulation_xml_overrides_output_prefix():
    ff, motion = _basic_blocks()
    xml = simulation_xml(_atoms(), ff, motion, prefix="run1")
    assert "prefix='run1'" in xml
    xml_parse_string(xml)


def test_simulation_xml_rejects_invalid_output_block_with_prefix():
    ff, motion = _basic_blocks()
    with pytest.raises(ValueError, match="valid 'output' block"):
        simulation_xml(_atoms(), ff, motion, output="<notoutput/>", prefix="x")


def test_simulation_xml_list_of_structures_sets_nbeads():
    ff, motion = _basic_blocks()
    structures = [_atoms(), _atoms(), _atoms()]
    xml = simulation_xml(structures, ff, motion, temperature=300.0)
    # path integral: nbeads propagates to the beads block and the initializer
    assert "nbeads='3'" in xml
    xml_parse_string(xml)

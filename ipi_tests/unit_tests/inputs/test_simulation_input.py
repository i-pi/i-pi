"""The attributes of <simulation>."""

import pytest

from ipi.inputs.simulation import InputSimulation
from ipi_tests.unit_tests.engine.inputs import harmonic_input


@pytest.mark.parametrize(
    "attributes,interleave",
    [
        ("", True),
        ("interleave='false'", False),
        ("threading='false'", False),
        ("interleave='true' threading='false'", True),
    ],
)
def test_threading_is_a_deprecated_name_of_interleave(
    attributes, interleave, simulation
):
    """`threading` sets `interleave` if it is not given, and the checkpoints
    are written with `interleave`."""

    xml = harmonic_input().replace('interleave="False"', attributes)
    with simulation(xml) as sim:
        assert sim.interleave == interleave
        stored = InputSimulation()
        stored.store(sim)
        checkpoint = stored.write("simulation")
    assert "threading" not in checkpoint
    assert ("interleave" in checkpoint) == (not interleave)


@pytest.mark.parametrize("override", [None, "%.16e"])
def test_floatformat_in_properties_and_xyz(override, simulation, tmp_path, monkeypatch):
    """One simulation setting controls properties, XYZ coordinates and the cell."""
    from ipi import ipi_global_settings

    monkeypatch.setitem(ipi_global_settings, "floatformat", "%16.8e")
    fmt = override or "%16.8e"
    attr = "" if override is None else f' floatformat="{override}"'
    xml = harmonic_input().replace(' floatformat="%24.16e"', attr)
    xml = xml.replace(
        "<total_steps>",
        '<output prefix="test"><properties filename="out" stride="1">[volume]</properties>'
        '<trajectory filename="pos" stride="1">positions</trajectory></output><total_steps>',
    )
    with simulation(xml) as sim:
        sim.run(1)
    rows = (tmp_path / "test.out").read_text().splitlines()
    assert (
        next(row.strip() for row in rows if not row.startswith("#"))
        == (fmt % 8000.0).strip()
    )
    xyz = (tmp_path / "test.pos_0.xyz").read_text().splitlines()
    assert xyz[1].split("CELL(abcABC):")[1].split()[:6] == [
        (fmt % x).strip() for x in [20, 20, 20, 90, 90, 90]
    ]
    assert xyz[2].split() == ["H"] + [(fmt % x).strip() for x in [0.2, -0.3, 0.1]]

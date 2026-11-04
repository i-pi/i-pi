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

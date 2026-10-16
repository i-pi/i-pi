"""Accepted exchanges preserve reduced momenta; rejected exchanges restore state."""

import copy
import xml.etree.ElementTree as ET
from types import SimpleNamespace
import numpy as np
import pytest
from ipi.engine.smotion.remd import ReplicaExchange, motion_scale
from ipi_tests.unit_tests.engine.test_simulation_invariants import harmonic_input


@pytest.mark.parametrize("accept", [False, True])
@pytest.mark.parametrize("rescale", [False, True])
def test_exchange_acceptance_or_rejection_has_consistent_state(
    simulation, accept, rescale
):
    root = ET.fromstring(harmonic_input())
    second = copy.deepcopy(root.find("system"))
    second.find("ensemble/temperature").text = "600"
    root.append(second)
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        systems = sim.syslist
        # Place the hotter replica high on the potential, making a swap unfavorable.
        systems[1].beads.q[:] = 10
        momenta = [s.beads.p.copy() for s in systems]
        temperatures = [s.ensemble.temp for s in systems]
        energies = [s.ensemble.econs for s in systems]
        attempts = iter([0.0, 0.0 if accept else 0.999])

        class Draws:
            @property
            def u(self):
                return next(attempts)

        remd = ReplicaExchange(krescale=rescale)
        remd.bind(systems, Draws(), sim.output_maker)
        # Use a moderate displacement so acceptance is nonzero but below .999.
        systems[1].beads.q[:] = 0.3
        energies = [s.ensemble.econs for s in systems]
        remd.step(0)
        if accept:
            np.testing.assert_array_equal(remd.repindex, [1, 0])
            for i, s in enumerate(systems):
                assert s.ensemble.temp == temperatures[1 - i]
                factor = (
                    np.sqrt(temperatures[1 - i] / temperatures[i]) if rescale else 1.0
                )
                np.testing.assert_allclose(s.beads.p, momenta[i] * factor)
                assert s.ensemble.econs == pytest.approx(energies[i])
        else:
            np.testing.assert_array_equal(remd.repindex, [0, 1])
            for i, s in enumerate(systems):
                assert s.ensemble.temp == temperatures[i]
                np.testing.assert_allclose(s.beads.p, momenta[i])


def test_nested_thermostat_scaling_reaches_all_auxiliary_momenta():
    states = [SimpleNamespace(s=np.ones((2, 3)) * v) for v in [1.0, 2.0, 3.0]]
    motion = SimpleNamespace(
        mlist=[SimpleNamespace(thermostat=SimpleNamespace(tlist=states[:2]))],
        barostat=SimpleNamespace(thermostat=states[2]),
    )
    motion_scale(motion, 0.5)
    for v, state in enumerate(states, 1):
        np.testing.assert_allclose(state.s, 0.5 * v)

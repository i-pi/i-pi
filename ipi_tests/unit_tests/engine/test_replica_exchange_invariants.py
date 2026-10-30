"""Accepted exchanges preserve reduced momenta; rejected exchanges restore state."""

import copy
import xml.etree.ElementTree as ET
from types import SimpleNamespace
import numpy as np
import pytest
from ipi.engine.smotion.remd import ReplicaExchange, motion_scale
from ipi_tests.unit_tests.engine.inputs import harmonic_input


@pytest.mark.parametrize("accept", [False, True])
@pytest.mark.parametrize("rescale", [False, True])
@pytest.mark.parametrize("displacement", [0.2, 0.4])
def test_exchange_acceptance_or_rejection_has_consistent_state(
    simulation, accept, rescale, displacement
):
    """Draws just below or above the analytic acceptance probability."""
    # Atomic units: kB=1. U(q)=sum(q**2)/4 and K(p)=sum(p**2/(2*m)).
    temperatures = [0.2, 0.4]
    positions = [np.zeros((1, 6)), np.full((1, 6), displacement)]
    momenta = [
        np.array([[0.1, 0.2, -0.1, 0.3, -0.2, 0.1]]),
        np.array([[0.2, 0.4, -0.2, 0.6, -0.4, 0.2]]),
    ]
    masses = np.repeat([2.0, 3.0], 3)
    potential = np.array([0.25 * np.sum(q**2) for q in positions])
    kinetic = np.array([np.sum(p**2 / (2 * masses)) for p in momenta])
    # Sugita, Kitao and Okamoto: https://arxiv.org/pdf/cond-mat/0009120
    # Eqs. (5), (8), (10)-(13): momentum rescaling cancels the kinetic
    # terms in the Boltzmann-weight ratio. Without it, retain those terms.
    exchange_energy = potential if rescale else potential + kinetic
    probability = min(
        1.0,
        np.exp(
            (1 / temperatures[0] - 1 / temperatures[1])
            * (exchange_energy[0] - exchange_energy[1])
        ),
    )
    acceptance_draw = probability - 1e-6 if accept else probability + 1e-6
    root = ET.fromstring(harmonic_input())
    second = copy.deepcopy(root.find("system"))
    for system_input, temperature, q, p in zip(
        [root.find("system"), second], temperatures, positions, momenta
    ):
        temperature_input = system_input.find("ensemble/temperature")
        temperature_input.attrib.clear()
        temperature_input.text = str(temperature)
        system_input.find("beads/q").text = str(q.ravel().tolist())
        system_input.find("beads/p").text = str(p.ravel().tolist())
    root.append(second)
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        systems = sim.syslist
        energies = [s.ensemble.econs for s in systems]
        # First draw selects the exchange attempt; second tests its probability.
        attempts = iter([0.0, acceptance_draw])

        class Draws:
            @property
            def u(self):
                return next(attempts)

        remd = ReplicaExchange(krescale=rescale)
        remd.bind(systems, Draws(), sim.output_maker)
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

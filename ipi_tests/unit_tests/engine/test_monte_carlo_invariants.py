"""Controlled Monte Carlo proposals preserve the required physical quantities."""

import xml.etree.ElementTree as ET

import numpy as np
import pytest

from ipi.inputs.motion.motion import InputMotion
from ipi_tests.unit_tests.engine.inputs import harmonic_input


class Proposal:
    """Select the two atoms once; control only the random acceptance decision."""

    def __init__(self, draw, attempts=1):
        self.indices = iter([0, 0, 1])
        self.u = draw
        self.attempts = attempts

    def poisson(self, mean):
        return self.attempts

    def integers(self, low, high):
        return next(self.indices)


@pytest.mark.parametrize("mode", ["alchemy", "atomswap"])
@pytest.mark.parametrize("nbeads", [1, 4])
@pytest.mark.parametrize("attempts", [0, 1])
def test_exchange_preserves_energy_and_species(simulation, mode, nbeads, attempts):
    root = ET.fromstring(harmonic_input(nbeads=nbeads))
    motion = root.find("system/motion")
    motion.clear()
    motion.set("mode", mode)
    ET.SubElement(motion, "fixcom").text = "False"
    options = ET.SubElement(motion, mode)
    ET.SubElement(options, "names").text = "[H, O]"
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        system = sim.syslist[0]
        beads = system.beads
        q, p, m = beads.q.copy(), beads.p.copy(), beads.m.copy()
        kinetic = np.sum(p**2 / np.repeat(m, 3)) / 2
        spring = beads.vpath * system.nm.omegan2
        potential = system.forces.pot
        system.motion.prng = Proposal(0.0, attempts)
        system.motion.step(0)
        assert np.sum(beads.p**2 / np.repeat(beads.m, 3)) / 2 == pytest.approx(kinetic)
        assert system.forces.pot == pytest.approx(potential)
        if mode == "alchemy" or not attempts or nbeads == 1:
            assert (
                beads.vpath * system.nm.omegan2 + system.motion.ealc
                == pytest.approx(spring)
            )
        if attempts and mode == "alchemy":
            np.testing.assert_array_equal(beads.names, ["O", "H"])
            np.testing.assert_array_equal(beads.m, m[::-1])
            np.testing.assert_allclose(beads.p, p * np.repeat(np.sqrt(m[::-1] / m), 3))
            np.testing.assert_array_equal(beads.q, q)
        elif attempts:
            np.testing.assert_array_equal(
                beads.q, q.reshape(nbeads, 2, 3)[:, ::-1].reshape(nbeads, 6)
            )
            np.testing.assert_array_equal(beads.p, p)
        else:
            np.testing.assert_array_equal(beads.q, q)
            np.testing.assert_array_equal(beads.p, p)
            np.testing.assert_array_equal(beads.m, m)
        # The changed species/masses and accumulated exchange energy belong in restart.

        stored = InputMotion()
        stored.store(system.motion)
        restored = stored.fetch()
        assert restored.ealc == pytest.approx(system.motion.ealc)
        np.testing.assert_array_equal(restored.names, ["H", "O"])


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="AtomSwap acceptance omits the ring-polymer spring-energy change for unequal masses",
)
def test_ring_polymer_swap_includes_spring_energy(simulation):
    root = ET.fromstring(harmonic_input(nbeads=4))
    motion = root.find("system/motion")
    motion.clear()
    motion.set("mode", "atomswap")
    ET.SubElement(motion, "fixcom").text = "False"
    ET.SubElement(ET.SubElement(motion, "atomswap"), "names").text = "[H, O]"
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        s = sim.syslist[0]
        q = np.zeros((4, 6))
        q[:, 0] = [-10, 10, -10, 10]
        s.beads._q[:] = q
        delta = (
            0.5
            * (3.0 - 2.0)
            * s.nm.omegan2
            * np.sum((q[:, 0] - np.roll(q[:, 0], 1)) ** 2)
        )
        assert np.exp(-delta / (4 * s.ensemble.temp)) < 0.999
        s.motion.prng = Proposal(0.999)
        s.motion.step(0)
        np.testing.assert_array_equal(s.beads.q, q)

"""Controlled Monte Carlo proposals preserve the required physical quantities."""

import xml.etree.ElementTree as ET

import numpy as np
import pytest

from ipi.inputs.motion.motion import InputMotion
from ipi.utils.array_backend import to_numpy, xp
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
        q, p, m = [to_numpy(x).copy() for x in (beads.q, beads.p, beads.m)]
        kinetic = np.sum(p**2 / np.repeat(m, 3)) / 2
        spring = float(system.nm.vspring)
        potential = float(system.forces.pot)
        system.motion.prng = Proposal(0.0, attempts)
        system.motion.step(0)
        q1, p1, m1 = [to_numpy(x) for x in (beads.q, beads.p, beads.m)]
        assert np.sum(p1**2 / np.repeat(m1, 3)) / 2 == pytest.approx(kinetic)
        assert float(system.forces.pot) == pytest.approx(potential)
        assert float(system.nm.vspring + system.motion.ealc) == pytest.approx(spring)
        if attempts and mode == "alchemy":
            np.testing.assert_array_equal(beads.names, ["O", "H"])
            np.testing.assert_array_equal(m1, m[::-1])
            np.testing.assert_allclose(p1, p * np.repeat(np.sqrt(m[::-1] / m), 3))
            np.testing.assert_array_equal(q1, q)
        elif attempts:
            np.testing.assert_array_equal(
                q1, q.reshape(nbeads, 2, 3)[:, ::-1].reshape(nbeads, 6)
            )
            np.testing.assert_array_equal(p1, p)
        else:
            np.testing.assert_array_equal(q1, q)
            np.testing.assert_array_equal(p1, p)
            np.testing.assert_array_equal(m1, m)
        # The changed species/masses and accumulated exchange energy belong in restart.

        stored = InputMotion()
        stored.store(system.motion)
        restored = stored.fetch()
        assert float(restored.ealc) == pytest.approx(float(system.motion.ealc))
        np.testing.assert_array_equal(restored.names, ["H", "O"])


@pytest.mark.parametrize("mode", ["alchemy", "atomswap"])
@pytest.mark.parametrize("paths", ["open_paths", "bosons"])
def test_exchange_refuses_open_paths_and_bosons(simulation, mode, paths):
    """The change of spring energy assumes closed paths of distinguishable atoms."""
    root = ET.fromstring(harmonic_input(nbeads=4))
    motion = root.find("system/motion")
    motion.clear()
    motion.set("mode", mode)
    ET.SubElement(motion, "fixcom").text = "False"
    ET.SubElement(ET.SubElement(motion, mode), "names").text = "[H, O]"
    normal_modes = ET.SubElement(root.find("system"), "normal_modes")
    ET.SubElement(normal_modes, paths).text = "[0]"
    with pytest.raises(ValueError, match="not implemented for open paths or bosons"):
        with simulation(ET.tostring(root, encoding="unicode")):
            pass


@pytest.mark.parametrize("accepted", [False, True])
def test_ring_polymer_swap_includes_spring_energy(simulation, accepted):
    """Swapping atoms with different masses changes the spring energy, which
    enters the acceptance and the conserved quantity"""
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
        s.beads.q_[:] = xp.asarray(q)
        delta = (
            0.5
            * (3.0 - 2.0)
            * float(s.nm.omegan2)
            * np.sum((q[:, 0] - np.roll(q[:, 0], 1)) ** 2)
        )
        assert np.exp(-delta / (4 * float(s.ensemble.temp))) < 0.999
        s.motion.prng = Proposal(0.0 if accepted else 0.999)
        s.motion.step(0)
        if accepted:
            # the harmonic potential of the two atoms is the same
            swapped = q.reshape(4, 2, 3)[:, ::-1].reshape(4, 6)
            np.testing.assert_array_equal(to_numpy(s.beads.q), swapped)
            assert float(s.motion.ealc) == pytest.approx(-delta)
        else:
            np.testing.assert_array_equal(to_numpy(s.beads.q), q)
            assert float(s.motion.ealc) == 0


@pytest.mark.parametrize("nbeads", [1, 4])
def test_dynamics_after_alchemical_exchange_uses_new_masses(simulation, nbeads):
    """NVE dynamics conserves the energy across an exchange of masses"""
    root = ET.fromstring(harmonic_input(nbeads=nbeads))
    system = root.find("system")
    dynamics = system.find("motion")
    system.remove(dynamics)
    motion = ET.SubElement(system, "motion", mode="multi")
    motion.append(dynamics)
    alchemy = ET.SubElement(motion, "motion", mode="alchemy")
    ET.SubElement(alchemy, "fixcom").text = "False"
    ET.SubElement(ET.SubElement(alchemy, "alchemy"), "names").text = "[H, O]"
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        s = sim.syslist[0]
        econs = []
        for step in range(10):
            # exchanges the two atoms in the first step only
            s.motion.mlist[1].prng = Proposal(0.0, int(step == 0))
            sim.run_step(step)
            econs.append(float(s.ensemble.econs))
        np.testing.assert_array_equal(s.beads.names, ["O", "H"])
        # the drift is 5e-3 if the positions are propagated with the old masses
        assert np.ptp(econs) < 1e-5 * abs(econs[0])

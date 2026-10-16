"""Higher-order and constrained-centroid dynamics preserve state through restart."""

import xml.etree.ElementTree as ET

import numpy as np
import pytest
from ipi_tests.unit_tests.engine.test_simulation_invariants import harmonic_input
from ipi_tests.unit_tests.engine.test_pressure_restart_reference import pressure_input


@pytest.mark.parametrize("mode", ["sc", "scnpt", "nvt-cc"])
@pytest.mark.parametrize("splitting", ["obabo", "baoab"])
def test_extended_dynamics_restart_matches_continuous_run(simulation, mode, splitting):
    root = ET.fromstring(
        pressure_input("isotropic", 4)
        if mode == "scnpt"
        else harmonic_input(nbeads=4, thermostat="pile_l", timestep=0.001)
    )
    dyn = root.find("system/motion/dynamics")
    dyn.set("mode", mode)
    if mode == "scnpt":
        dyn.find("barostat").set("mode", "sc-isotropic")
    dyn.set("splitting", splitting)
    xml = ET.tostring(root, encoding="unicode")

    def snapshot(sim):
        s = sim.syslist[0]
        return [
            s.beads.q.copy(),
            s.beads.p.copy(),
            s.cell.h.copy(),
            s.ensemble.econs,
            s.ensemble.time,
        ]

    with simulation(xml) as sim:
        initial_centroid = sim.syslist[0].beads.qc.copy()
        sim.run(steps=6, write_outputs=False)
        expected = snapshot(sim)
        if mode == "nvt-cc":
            np.testing.assert_allclose(
                sim.syslist[0].beads.qc, initial_centroid, atol=1e-13
            )
    with simulation(xml) as sim:
        sim.run(steps=2, write_outputs=False)
    with open("RESTART") as f:
        restart = f.read()
    with simulation(restart) as sim:
        sim.run(steps=4, write_outputs=False)
        actual = snapshot(sim)
    for a, b in zip(actual, expected):
        np.testing.assert_allclose(a, b, rtol=1e-10, atol=1e-12)

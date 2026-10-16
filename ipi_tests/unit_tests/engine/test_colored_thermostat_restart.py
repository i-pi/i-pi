"""Restart must preserve random draws and auxiliary thermostat momenta."""

import xml.etree.ElementTree as ET

import numpy as np
import pytest

from ipi_tests.unit_tests.engine.test_simulation_invariants import harmonic_input


def thermostat_input(mode, explicit_covariance=False):
    root = ET.fromstring(harmonic_input(nbeads=4, thermostat=mode, timestep=0.001))
    bath = root.find("system/motion/dynamics/thermostat")
    if "gle" in mode:
        drift = np.array([[2.0, -0.2], [0.2, 1.0]])
        covariance = np.eye(2) * 0.004
        if mode != "gle":
            drift = np.tile(drift, (4, 1, 1))
            covariance = np.tile(covariance, (4, 1, 1))
        ET.SubElement(bath, "A", shape=str(drift.shape)).text = str(
            drift.ravel().tolist()
        )
        if explicit_covariance:
            ET.SubElement(bath, "C", shape=str(covariance.shape)).text = str(
                covariance.ravel().tolist()
            )
    if mode == "multi":
        for child_mode in ["langevin", "svr"]:
            child = ET.SubElement(bath, "thermostat", mode=child_mode)
            ET.SubElement(child, "tau").text = "0.5"
    return ET.tostring(root, encoding="unicode")


@pytest.mark.parametrize(
    "mode,explicit_covariance",
    [(m, c) for m in ["gle", "nm_gle", "nm_gle_g"] for c in [False, True]]
    + [(m, False) for m in ["multi", "cl", "ffl"]],
)
def test_restart_preserves_thermostat_state(simulation, mode, explicit_covariance):
    xml = thermostat_input(mode, explicit_covariance)

    def snapshot(sim):
        system = sim.syslist[0]
        bath = system.motion.thermostat
        values = [system.beads.q.copy(), system.beads.p.copy(), bath.ethermo]
        if "gle" in mode:
            values.append(bath.s.copy())
        return values

    with simulation(xml) as sim:
        sim.run(steps=6, write_outputs=False)
        expected = snapshot(sim)
    with simulation(xml) as sim:
        sim.run(steps=3, write_outputs=False)
    with open("RESTART") as f:
        restart = f.read()
    with simulation(restart) as sim:
        sim.run(steps=3, write_outputs=False)
        actual = snapshot(sim)
    for value, reference in zip(actual, expected):
        np.testing.assert_allclose(value, reference, rtol=1e-11, atol=1e-12)

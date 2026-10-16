"""Combined motions retain ramp progress and physical time across restart."""

import xml.etree.ElementTree as ET
import numpy as np
import pytest
from ipi_tests.unit_tests.engine.test_simulation_invariants import harmonic_input


@pytest.mark.parametrize("kind,attribute", [("t", "temp"), ("p", "pext")])
@pytest.mark.parametrize("logscale", [False, True])
def test_ramp_and_dynamics_restart(simulation, kind, attribute, logscale):
    root = ET.fromstring(harmonic_input(timestep=0.01))
    system = root.find("system")
    dynamics = system.find("motion")
    system.remove(dynamics)
    multi = ET.SubElement(system, "motion", mode="multi")
    ramp = ET.SubElement(multi, "motion", mode=kind + "_ramp")
    options = ET.SubElement(ramp, kind + "_ramp")
    for key, value in {
        kind + "_start": 0.1,
        kind + "_end": 1.6,
        "total_steps": 4,
        "logscale": str(logscale),
    }.items():
        ET.SubElement(options, key).text = str(value)
    multi.append(dynamics)
    xml = ET.tostring(root, encoding="unicode")
    expected_values = (
        [0.2, 0.4, 0.8, 1.6, 1.6, 1.6]
        if logscale
        else [0.475, 0.85, 1.225, 1.6, 1.6, 1.6]
    )
    with simulation(xml) as sim:
        for expected in expected_values:
            sim.run(steps=1, write_outputs=False)
            assert getattr(sim.syslist[0].ensemble, attribute) == pytest.approx(
                expected
            )
        expected_q = sim.syslist[0].beads.q.copy()
        expected_p = sim.syslist[0].beads.p.copy()
        assert sim.syslist[0].ensemble.time == pytest.approx(0.06)
    with simulation(xml) as sim:
        sim.run(steps=2, write_outputs=False)
    with open("RESTART") as f:
        restart = f.read()
    with simulation(restart) as sim:
        assert sim.syslist[0].motion.mlist[0].current_step == 2
        sim.run(steps=4, write_outputs=False)
        assert getattr(sim.syslist[0].ensemble, attribute) == pytest.approx(1.6)
        np.testing.assert_allclose(sim.syslist[0].beads.q, expected_q, atol=1e-12)
        np.testing.assert_allclose(sim.syslist[0].beads.p, expected_p, atol=1e-12)
        assert sim.syslist[0].ensemble.time == pytest.approx(0.06)

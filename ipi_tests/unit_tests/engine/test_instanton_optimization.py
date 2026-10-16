"""Instanton optimizers recover known stationary points of quadratic actions."""

import xml.etree.ElementTree as ET

import numpy as np
import pytest
from ipi_tests.unit_tests.engine.test_mode_and_path_reference import motion_input


@pytest.mark.parametrize("optimizer", ["nichols", "NR", "lanczos", "lbfgs"])
def test_instanton_optimizer_reaches_quadratic_stationary_point(simulation, optimizer):
    count = 4 if optimizer == "lbfgs" else 1
    options = ET.Element(
        "instanton", mode="splitting" if optimizer == "lbfgs" else "rate"
    )
    ET.SubElement(options, "opt").text = optimizer
    ET.SubElement(options, "hessian_asr").text = "none"
    ET.SubElement(options, "energy_shift").text = "0"
    h = np.diag([0.5 if optimizer == "lbfgs" else -0.5, 0.7, 0.9])
    ET.SubElement(options, "hessian", shape=f"(3,{3*count})").text = str(
        np.tile(h, (1, count)).ravel().tolist()
    )
    root = motion_input("instanton", options, nbeads=count)
    root.remove(root.find("ffdirect"))
    ff = ET.SubElement(root, "ffdebye", name="harmonic", pbc="False")
    ET.SubElement(ff, "hessian", shape="(3,3)").text = str(h.ravel().tolist())
    ET.SubElement(ff, "x_reference", shape="(3)").text = "[0,0,0]"
    beads = root.find("system/beads")
    beads.set("natoms", "1")
    q = np.tile([0.02, -0.03, 0.04], (count, 1))
    q[:, 0] += np.arange(count) * 0.01
    for tag, shape, value in [
        ("q", f"({count},3)", q.ravel().tolist()),
        ("p", f"({count},3)", [0.0] * (3 * count)),
        ("m", "(1)", [1.0]),
        ("names", "(1)", ["H"]),
    ]:
        beads.find(tag).set("shape", shape)
        beads.find(tag).text = str(value)
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        s = sim.syslist[0]
        s.motion.exit_on_convergence = False
        for step in range(50):
            s.motion.step(step)
            if s.motion.finished:
                break
        np.testing.assert_allclose(s.beads.q, 0, atol=3e-4)
        np.testing.assert_allclose(s.forces.f, -s.beads.q @ h, atol=1e-12)
        assert s.motion.finished

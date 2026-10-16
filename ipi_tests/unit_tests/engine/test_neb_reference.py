"""A harmonic transverse potential has a straight minimum-energy path."""

import xml.etree.ElementTree as ET

import numpy as np
import pytest
from ipi_tests.unit_tests.engine.test_mode_and_path_reference import motion_input


@pytest.mark.parametrize("optimizer", ["fire", "damped_bfgs"])
@pytest.mark.parametrize("tangent", ["plain", "improved"])
def test_neb_relaxes_transverse_displacement_with_fixed_endpoints(
    simulation, optimizer, tangent
):
    options = ET.Element("neb_optimizer", mode=optimizer)
    ET.SubElement(options, "tangent").text = tangent
    root = motion_input("neb", options, nbeads=5)
    q = np.zeros((5, 6))
    q[:, 0] = np.linspace(-1, 1, 5)
    q[1:-1, 1] = 0.1
    root.find("system/beads/q").text = str(q.ravel().tolist())
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        s = sim.syslist[0]
        s.motion.exit_on_convergence = False
        for step in range(150):
            s.motion.step(step)
            if s.motion.finished:
                break
        np.testing.assert_array_equal(s.beads.q[[0, -1]], q[[0, -1]])
        assert np.max(np.abs(s.beads.q[1:-1, 1:])) < 1e-4
        np.testing.assert_allclose(s.forces.f, -0.5 * s.beads.q, atol=1e-12)
        assert np.isfinite(s.beads.q).all()


@pytest.mark.parametrize("optimizer", ["fire", "damped_bfgs"])
def test_climbing_image_finds_quadratic_saddle(simulation, optimizer):
    options = ET.Element("neb_optimizer", mode=optimizer)
    ET.SubElement(options, "tangent").text = "improved"
    ET.SubElement(options, "use_climb").text = "True"
    root = motion_input("neb", options, nbeads=5)
    root.remove(root.find("ffdirect"))
    ff = ET.SubElement(root, "ffdebye", name="harmonic", pbc="False")
    ET.SubElement(ff, "hessian", shape="(3,3)").text = str(
        np.diag([-1.0, 0.5, 0.5]).ravel().tolist()
    )
    ET.SubElement(ff, "x_reference", shape="(3)").text = "[0,0,0]"
    beads = root.find("system/beads")
    beads.set("natoms", "1")
    q = np.zeros((5, 3))
    q[:, 0] = [-1, -0.5, 0.05, 0.5, 1]
    q[1:-1, 1] = 0.05
    for tag, shape, value in [
        ("q", "(5,3)", q.ravel().tolist()),
        ("p", "(5,3)", [0.0] * 15),
        ("m", "(1)", [1.0]),
        ("names", "(1)", ["H"]),
    ]:
        beads.find(tag).set("shape", shape)
        beads.find(tag).text = str(value)
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        s = sim.syslist[0]
        s.motion.exit_on_convergence = False
        for step in range(350):
            s.motion.step(step)
            if s.motion.finished:
                break
        assert s.motion.cl_indx == 2
        np.testing.assert_allclose(s.beads.q[2], 0, atol=2e-4)
        np.testing.assert_array_equal(s.beads.q[[0, -1]], q[[0, -1]])
        assert s.forces.pots[2] == pytest.approx(0, abs=1e-7)

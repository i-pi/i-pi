"""Optimizers and path methods recover known minima, saddles and straight paths."""

import xml.etree.ElementTree as ET
from types import SimpleNamespace

import numpy as np
import pytest

from ipi.engine.motion.stringmep import (
    GradientMapper,
    StringGradientMapper,
    spline_derv,
    spline_resample,
)
from ipi_tests.unit_tests.engine.inputs import harmonic_input, motion_input


@pytest.mark.parametrize(
    "method", ["sd", "cg", "bfgs", "bfgstrm", "lbfgs", "damped_bfgs"]
)
@pytest.mark.parametrize("selection", [None, "fixatoms", "fixatoms_dof"])
def test_geometry_optimization_finds_harmonic_minimum(
    simulation, method, selection, request
):
    """Minimize free coordinates while preserving any fixed atoms or coordinates."""
    if method == "damped_bfgs" and selection == "fixatoms_dof":
        request.node.add_marker(
            pytest.mark.xfail(
                strict=True,
                raises=ValueError,
                reason="Damped_BFGS reshapes five free coordinates into triples",
            )
        )
    root = ET.fromstring(harmonic_input())
    motion = root.find("system/motion")
    motion.clear()
    motion.set("mode", "minimize")
    ET.SubElement(motion, "fixcom").text = "False"
    if selection:
        ET.SubElement(motion, selection, shape="(1)").text = "[0]"
    optimizer = ET.SubElement(motion, "optimizer", mode=method)
    ET.SubElement(optimizer, "exit_on_convergence").text = "False"
    fixed = {None: [], "fixatoms": [0, 1, 2], "fixatoms_dof": [0]}[selection]
    free = sorted(set(range(6)) - set(fixed))
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        system = sim.syslist[0]
        initial = system.beads.q.copy()
        sim.run(steps=150 if selection else 20, write_outputs=False)
        np.testing.assert_array_equal(system.beads.q[:, fixed], initial[:, fixed])
        np.testing.assert_allclose(
            system.beads.q[:, free], 0, atol=2e-4 if selection else 1e-5
        )
        if selection:
            assert system.forces.pot == pytest.approx(
                0.25 * np.sum(initial[:, fixed] ** 2), abs=1e-7
            )
        else:
            assert system.forces.pot < 1e-10


@pytest.mark.parametrize("optimizer", ["fire", "damped_bfgs"])
@pytest.mark.parametrize("tangent", ["plain", "improved"])
def test_neb_relaxes_transverse_displacement(simulation, optimizer, tangent):
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


@pytest.mark.parametrize("count", [4, 7, 11])
def test_string_resampling_of_straight_line(count):
    t = np.linspace(0, 1, count) ** 2
    start = np.array([0.2, -0.1, 0.3])
    direction = np.array([2.0, 1.0, -1.0])
    q = start + t[:, None] * direction
    sampled = spline_resample(q, count)
    np.testing.assert_allclose(
        sampled, start + np.linspace(0, 1, count)[:, None] * direction, atol=1e-12
    )
    np.testing.assert_allclose(
        spline_derv(q, count),
        np.tile(direction / np.linalg.norm(direction), (count, 1)),
        atol=1e-12,
    )


@pytest.mark.parametrize("projected", [False, True])
def test_path_mapper_harmonic_gradient_and_tangent_projection(simulation, projected):

    with simulation(harmonic_input(nbeads=7)) as sim:
        system = sim.syslist[0]
        q = np.zeros((7, 6))
        q[:, 0] = np.linspace(-1, 1, 7)
        q[:, 1] = 0.2
        system.beads.q[:] = q
        mapper = (StringGradientMapper if projected else GradientMapper)()
        mapper.bind(
            SimpleNamespace(
                beads=system.beads,
                cell=system.cell,
                forces=system.forces,
                fixatoms_dof=np.array([], dtype=int),
            )
        )
        energy, gradient = mapper(q[1:-1])
        assert energy == pytest.approx(0.25 * np.sum(q[1:-1] ** 2))
        expected = 0.5 * q[1:-1].copy()
        if projected:
            expected[:, 0] = 0
        np.testing.assert_allclose(gradient, expected, atol=1e-12)

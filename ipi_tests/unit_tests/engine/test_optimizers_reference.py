"""Optimizers and path methods recover known minima, saddles and straight paths."""

import xml.etree.ElementTree as ET
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from ipi.engine.motion.geop import GeopMotion
from ipi.engine.motion.stringmep import (
    GradientMapper,
    StringGradientMapper,
    spline_derv,
    spline_resample,
)
from ipi_tests.unit_tests.engine.inputs import harmonic_input, motion_input


@pytest.mark.parametrize("trust_radius", [[], [0.2]])
def test_geop_reset_restores_the_trust_radius(trust_radius):
    """Resetting a BFGSTRM optimization (as the Al6xxx KMC does at each step)
    restores the initial trust radius, in the array shared with the optimizer"""

    geop = GeopMotion(mode="bfgstrm", tr_trm=np.array(trust_radius))
    tr = geop.tr
    tr[0] = 1.5
    geop.reset()
    assert geop.tr is tr and tr[0] == (trust_radius or [0.4])[0]


def geop_input(method, selection=None, nbeads=1):
    """Optimization of the harmonic input with `method`, with the first atom or
    coordinate fixed if `selection` is fixatoms or fixatoms_dof."""
    root = ET.fromstring(harmonic_input(nbeads=nbeads))
    motion = root.find("system/motion")
    motion.clear()
    motion.set("mode", "minimize")
    ET.SubElement(motion, "fixcom").text = "False"
    if selection:
        ET.SubElement(motion, selection, shape="(1)").text = "[0]"
    optimizer = ET.SubElement(motion, "optimizer", mode=method)
    ET.SubElement(optimizer, "exit_on_convergence").text = "False"
    return ET.tostring(root, encoding="unicode")


@pytest.mark.parametrize(
    "method", ["sd", "cg", "bfgs", "bfgstrm", "lbfgs", "damped_bfgs"]
)
def test_geop_reset_restarts_a_converged_optimization(simulation, method):
    """After a reset (as the Al6xxx KMC does before each relaxation), a converged
    optimization relaxes a new structure, from the initial Hessian and with
    the same settings"""
    with simulation(geop_input(method)) as sim:
        system = sim.syslist[0]
        geop = system.motion
        for step in range(100):
            geop.step(step)
        assert geop.optimizer.converged
        if method == "lbfgs":
            settings = (geop.corrections, geop.scale)
        system.beads.q_[:] = system.beads.q + 0.3
        geop.reset()
        if method == "bfgs":
            np.testing.assert_array_equal(geop.invhessian, np.eye(6))
        if method in ["bfgstrm", "damped_bfgs"]:
            np.testing.assert_array_equal(geop.hessian, np.eye(6))
        if method == "lbfgs":
            assert (geop.corrections, geop.scale) == settings
        for step in range(100):
            geop.step(step)
        np.testing.assert_allclose(system.beads.q, 0, atol=1e-4)


@pytest.mark.parametrize("method", ["bfgs", "bfgstrm", "damped_bfgs"])
def test_geop_hessian_survives_a_restart(simulation, method):
    """The inverse Hessian of BFGS, and the Hessian of BFGSTRM and damped BFGS,
    are stored in the restart and read back"""

    def hessian(sim):
        geop = sim.syslist[0].motion
        return np.array(geop.invhessian if method == "bfgs" else geop.hessian)

    with simulation(geop_input(method)) as sim:
        sim.run(steps=3, write_outputs=False)
        stored = hessian(sim)
    # the optimizer has updated the initial (identity) Hessian
    assert not np.allclose(stored, np.eye(6))
    with simulation(Path("RESTART").read_text()) as sim:
        np.testing.assert_allclose(hessian(sim), stored)


@pytest.mark.parametrize(
    "method", ["sd", "cg", "bfgs", "bfgstrm", "lbfgs", "damped_bfgs"]
)
@pytest.mark.parametrize("selection", [None, "fixatoms", "fixatoms_dof"])
@pytest.mark.parametrize("nbeads", [1, 2])
def test_geometry_optimization_finds_harmonic_minimum(
    simulation, method, selection, nbeads
):
    """Minimize free coordinates while preserving any fixed atoms or coordinates."""
    fixed = {None: [], "fixatoms": [0, 1, 2], "fixatoms_dof": [0]}[selection]
    free = sorted(set(range(6)) - set(fixed))
    with simulation(geop_input(method, selection, nbeads)) as sim:
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


def test_conjugate_gradient_with_a_fixed_atom(simulation):
    """The conjugate directions are built from the forces on the free
    coordinates, so that the force on a fixed atom far from its minimum does not
    turn conjugate gradient into steepest descent"""
    root = ET.fromstring(harmonic_input())
    root.remove(root.find("ffdirect"))
    h = np.diag([0.5, 1.0, 2.0, 0.3, 1.1, 3.0])
    h[3, 4] = h[4, 3] = 0.2
    ff = ET.SubElement(root, "ffdebye", name="harmonic", pbc="False")
    ET.SubElement(ff, "hessian", shape="(6,6)").text = str(h.ravel().tolist())
    ET.SubElement(ff, "x_reference", shape="(6)").text = str([0.0] * 6)
    root.find("system/beads/q").text = str([2.0, 2.0, 2.0, 0.3, -0.4, 0.5])
    motion = root.find("system/motion")
    motion.clear()
    motion.set("mode", "minimize")
    ET.SubElement(motion, "fixcom").text = "False"
    ET.SubElement(motion, "fixatoms").text = "[0]"
    optimizer = ET.SubElement(motion, "optimizer", mode="cg")
    ET.SubElement(optimizer, "exit_on_convergence").text = "False"
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        s = sim.syslist[0]
        for step in range(10):
            s.motion.step(step)
        np.testing.assert_allclose(s.beads.q[0, 3:], 0, atol=1e-5)


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


SADDLE = np.diag([-1.0, 0.5, 0.5])


def saddle_input(optimizer, q):
    """Climbing-image NEB of one atom on the quadratic saddle SADDLE, from the
    path q."""
    options = ET.Element("neb_optimizer", mode=optimizer)
    ET.SubElement(options, "tangent").text = "improved"
    ET.SubElement(options, "use_climb").text = "True"
    nbeads = len(q)
    root = motion_input("neb", options, nbeads=nbeads)
    root.remove(root.find("ffdirect"))
    ff = ET.SubElement(root, "ffdebye", name="harmonic", pbc="False")
    ET.SubElement(ff, "hessian", shape="(3,3)").text = str(SADDLE.ravel().tolist())
    ET.SubElement(ff, "x_reference", shape="(3)").text = "[0,0,0]"
    beads = root.find("system/beads")
    beads.set("natoms", "1")
    for tag, shape, value in [
        ("q", f"({nbeads},3)", q.ravel().tolist()),
        ("p", f"({nbeads},3)", [0.0] * 3 * nbeads),
        ("m", "(1)", [1.0]),
        ("names", "(1)", ["H"]),
    ]:
        beads.find(tag).set("shape", shape)
        beads.find(tag).text = str(value)
    return ET.tostring(root, encoding="unicode")


@pytest.mark.parametrize("optimizer", ["fire", "damped_bfgs"])
def test_climbing_image_finds_quadratic_saddle(simulation, optimizer):
    q = np.zeros((5, 3))
    q[:, 0] = [-1, -0.5, 0.05, 0.5, 1]
    q[1:-1, 1] = 0.05
    with simulation(saddle_input(optimizer, q)) as sim:
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


def test_climbing_image_stores_its_physical_force(simulation):
    """After each climbing step, the stored force of the climbing bead is the
    physical one, not the climbing force (they differ away from the saddle)"""
    # climbs from the initial path, so that the climbing bead also has a
    # transverse force, and with FIRE, which climbs in several steps
    q = np.zeros((4, 3))
    q[:, 0] = [-1, -0.4, 0.3, 1]
    q[1:-1, 1] = 0.05
    root = ET.fromstring(saddle_input("fire", q))
    ET.SubElement(root.find("system/motion/neb_optimizer"), "stage").text = "climb"
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        s = sim.syslist[0]
        s.motion.exit_on_convergence = False
        transverse = 0
        for step in range(350):
            s.motion.step(step)
            if s.motion.stage != "climb":
                break
            i = s.motion.cl_indx
            force = np.asarray(s.forces.f[i])
            np.testing.assert_allclose(
                force, -SADDLE @ np.asarray(s.beads.q[i]), atol=1e-12
            )
            transverse += abs(force[1]) > 1e-3
        assert transverse > 1


@pytest.mark.parametrize("climb", [False, True])
def test_damped_bfgs_neb_with_a_fixed_coordinate(simulation, climb):
    """Damped BFGS keeps a fixed coordinate of an atom, relaxes the transverse
    one and, when climbing, finds the saddle"""
    q = np.zeros((4, 3))
    q[:, 0] = [-1, -0.4, 0.3, 1]
    q[1:-1, 1] = 0.05
    q[1:-1, 2] = 0.04
    root = ET.fromstring(saddle_input("damped_bfgs", q))
    root.find("system/motion/neb_optimizer/use_climb").text = str(climb)
    ET.SubElement(root.find("system/motion"), "fixatoms_dof").text = "[2]"
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        s = sim.syslist[0]
        s.motion.exit_on_convergence = False
        for step in range(300):
            s.motion.step(step)
            if s.motion.finished:
                break
        np.testing.assert_array_equal(s.beads.q[:, 2], q[:, 2])
        assert np.max(np.abs(s.beads.q[1:-1, 1])) < 1e-4
        if climb:
            assert abs(s.beads.q[s.motion.cl_indx, 0]) < 1e-4


def instanton_input(optimizer, h, q, update="powell", fixatoms=None):
    """Instanton optimization on U = q h q / 2 from the path q, with all atoms of
    unit mass."""
    count, ndof = q.shape
    natoms = ndof // 3
    options = ET.Element(
        "instanton", mode="splitting" if optimizer == "lbfgs" else "rate"
    )
    ET.SubElement(options, "opt").text = optimizer
    ET.SubElement(options, "hessian_asr").text = "none"
    ET.SubElement(options, "hessian_update").text = update
    ET.SubElement(options, "energy_shift").text = "0"
    ET.SubElement(options, "hessian", shape=f"({ndof},{ndof*count})").text = str(
        np.tile(h, (1, count)).ravel().tolist()
    )
    root = motion_input("instanton", options, nbeads=count)
    if fixatoms is not None:
        ET.SubElement(root.find("system/motion"), "fixatoms").text = str(fixatoms)
    root.remove(root.find("ffdirect"))
    ff = ET.SubElement(root, "ffdebye", name="harmonic", pbc="False")
    ET.SubElement(ff, "hessian", shape=f"({ndof},{ndof})").text = str(
        h.ravel().tolist()
    )
    ET.SubElement(ff, "x_reference", shape=f"({ndof})").text = str([0.0] * ndof)
    beads = root.find("system/beads")
    beads.set("natoms", str(natoms))
    for tag, shape, value in [
        ("q", f"({count},{ndof})", q.ravel().tolist()),
        ("p", f"({count},{ndof})", [0.0] * (ndof * count)),
        ("m", f"({natoms})", [1.0] * natoms),
        ("names", f"({natoms})", ["H"] * natoms),
    ]:
        beads.find(tag).set("shape", shape)
        beads.find(tag).text = str(value)
    return ET.tostring(root, encoding="unicode")


@pytest.mark.parametrize("optimizer", ["nichols", "NR", "lanczos", "lbfgs"])
def test_instanton_optimizer_reaches_quadratic_stationary_point(simulation, optimizer):
    count = 4 if optimizer == "lbfgs" else 1
    h = np.diag([0.5 if optimizer == "lbfgs" else -0.5, 0.7, 0.9])
    q = np.tile([0.02, -0.03, 0.04], (count, 1))
    q[:, 0] += np.arange(count) * 0.01
    with simulation(instanton_input(optimizer, h, q)) as sim:
        s = sim.syslist[0]
        s.motion.exit_on_convergence = False
        for step in range(50):
            s.motion.step(step)
            if s.motion.finished:
                break
        np.testing.assert_allclose(s.beads.q, 0, atol=3e-4)
        np.testing.assert_allclose(s.forces.f, -s.beads.q @ h, atol=1e-12)
        assert s.motion.finished


@pytest.mark.parametrize(
    "optimizer,update",
    [
        ("nichols", "powell"),
        ("nichols", "recompute"),
        ("NR", "powell"),
        ("NR", "recompute"),
        ("lanczos", "powell"),
        ("lanczos", "recompute"),
        ("lbfgs", "powell"),
    ],
)
def test_instanton_optimizer_with_a_fixed_atom(simulation, optimizer, update):
    """The free atom reaches the stationary point of the potential in which the
    other atom is held fixed"""
    count = 4 if optimizer == "lbfgs" else 1
    h = np.diag([0.5 if optimizer == "lbfgs" else -0.5, 0.7, 0.9, 1.0, 1.1, 1.2])
    h[1, 4] = h[4, 1] = 0.1
    q = np.tile([0.02, -0.03, 0.04, 0.01, 0.02, 0.03], (count, 1))
    q[:, 0] += np.arange(count) * 0.01
    with simulation(instanton_input(optimizer, h, q, update, fixatoms=[1])) as sim:
        s = sim.syslist[0]
        s.motion.exit_on_convergence = False
        for step in range(50):
            s.motion.step(step)
            if s.motion.finished:
                break
        np.testing.assert_array_equal(s.beads.q[:, 3:], q[:, 3:])
        free = [0, -h[1, 4] * q[0, 4] / h[1, 1], 0]
        np.testing.assert_allclose(
            s.beads.q[:, :3], np.tile(free, (count, 1)), atol=3e-4
        )
        assert s.motion.finished


def test_instanton_rejects_partly_fixed_atoms(simulation):
    """Instantons fix whole atoms: three fixed coordinates of two atoms are
    rejected (they used to fix both atoms)"""
    root = ET.fromstring(instanton_input("nichols", -np.eye(6), np.zeros((1, 6))))
    ET.SubElement(root.find("system/motion"), "fixatoms_dof").text = "[0, 1, 5]"
    with pytest.raises(SystemExit):
        with simulation(ET.tostring(root, encoding="unicode")):
            pass


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
        system.beads.q_[:] = q
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

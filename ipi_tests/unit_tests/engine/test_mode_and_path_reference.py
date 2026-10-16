"""Vibrational free energies and minimum-energy paths with analytic references."""

import xml.etree.ElementTree as ET
from types import SimpleNamespace

import numpy as np
import pytest

from ipi.engine.motion.vscf import IMF
from ipi.engine.motion.stringmep import spline_derv, spline_resample
from ipi_tests.unit_tests.engine.test_simulation_invariants import harmonic_input


def motion_input(mode, options, nbeads=1):
    root = ET.fromstring(harmonic_input(nbeads=nbeads))
    motion = root.find("system/motion")
    motion.clear()
    motion.set("mode", mode)
    ET.SubElement(motion, "fixcom").text = "False"
    motion.append(options)
    return root


@pytest.mark.parametrize("shift", [0.0, 0.3, -0.2])
def test_harmonic_schrodinger_spectrum_and_thermal_energy(shift):
    calc = IMF()
    calc.imm = SimpleNamespace(temp=0.2)
    x = np.linspace(-12, 12, 2001)
    basis = np.array([calc.psi(i, 1.0, 0.7, x) for i in range(20)])
    dx = x[1] - x[0]
    np.testing.assert_allclose(basis @ basis.T * dx, np.eye(20), atol=1e-10)
    free, energy, levels, vectors = calc.solve_schroedingers_equation(
        0.7, basis * np.sqrt(dx), np.full(x.size, shift), True
    )
    np.testing.assert_allclose(levels, 0.7 * (np.arange(20) + 0.5) + shift, atol=1e-10)
    np.testing.assert_allclose(vectors @ vectors.T, np.eye(20), atol=1e-12)
    assert free == pytest.approx(
        0.35 + 0.2 * np.log1p(-np.exp(-0.7 / 0.2)) + shift, abs=1e-10
    )
    assert energy == pytest.approx(0.35 / np.tanh(0.7 / 0.4) + shift, abs=1e-10)


@pytest.mark.parametrize("reuse_map", [False, True])
def test_independent_modes_recover_exact_harmonic_free_energy(simulation, reuse_map):
    options = ET.Element("normalmodes", mode="imf")
    for key, value in {
        "nint": 301,
        "nbasis": 16,
        "nevib": 20,
        "athresh": 1e-8,
        "print_1b_map": "True",
    }.items():
        ET.SubElement(options, key).text = str(value)
    masses = np.repeat([2.0, 3.0], 3)
    ET.SubElement(options, "dynmat", shape="(6,6)").text = str(
        np.diag(0.5 / masses).ravel().tolist()
    )
    root = motion_input("normalmodes", options)
    root.find("system/beads/q").text = str([0.0] * 6)
    temperature = root.find("system/ensemble/temperature")
    temperature.attrib.clear()
    temperature.text = ".2"
    xml = ET.tostring(root, encoding="unicode")
    for _ in range(2 if reuse_map else 1):
        with simulation(xml) as sim:
            sim.run(steps=6, write_outputs=False)
            calc = sim.syslist[0].motion.calc
            omega = np.sqrt(0.5 / masses)
            expected = np.sum(omega / 2 + 0.2 * np.log1p(-np.exp(-omega / 0.2)))
            assert calc.total_anhar_free_energy == pytest.approx(expected, abs=2e-7)
            assert calc.total_anhar_internal_energy == pytest.approx(
                np.sum(omega / 2 / np.tanh(omega / 0.4)), abs=2e-7
            )


@pytest.mark.parametrize("count", [4, 7, 11])
def test_string_resampling_preserves_straight_line_and_unit_tangent(count):
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
    from ipi.engine.motion.stringmep import GradientMapper, StringGradientMapper

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


@pytest.mark.parametrize("reuse_map", [False, True])
def test_vscf_pair_mapping_has_zero_harmonic_mode_coupling(simulation, reuse_map):
    options = ET.Element("normalmodes", mode="vscf")
    for key, value in {
        "nint": 51,
        "nbasis": 10,
        "nevib": 8,
        "grid": "False",
        "solve": "False",
    }.items():
        ET.SubElement(options, key).text = str(value)
    masses = np.repeat([2.0, 3.0], 3)
    ET.SubElement(options, "dynmat", shape="(6,6)").text = str(
        np.diag(0.5 / masses).ravel().tolist()
    )
    root = motion_input("normalmodes", options)
    root.find("system/beads/q").text = str([0.0] * 6)
    temperature = root.find("system/ensemble/temperature")
    temperature.attrib.clear()
    temperature.text = ".2"
    xml = ET.tostring(root, encoding="unicode")
    for _ in range(2 if reuse_map else 1):
        with simulation(xml) as sim:
            sim.run(steps=22, write_outputs=False)
            calc = sim.syslist[0].motion.calc
            # Last mapped pair: a separable quadratic has no two-mode correction.
            i, j = calc.inm, calc.jnm
            expected = np.add.outer(calc.v_indep_list[i], calc.v_indep_list[j]).ravel()
            np.testing.assert_allclose(
                np.sort(
                    np.loadtxt(sim.output_maker.prefix + "." + calc.v_coupled_filename)[
                        :, 2
                    ]
                ),
                np.sort(expected),
                atol=1e-12,
            )


def test_vscf_solver_with_exact_uncoupled_maps(simulation):
    options = ET.Element("normalmodes", mode="vscf")
    for key, value in {
        "nint": 101,
        "nbasis": 12,
        "nevib": 20,
        "grid": "True",
        "solve": "True",
    }.items():
        ET.SubElement(options, key).text = str(value)
    masses = np.repeat([2.0, 3.0], 3)
    ET.SubElement(options, "dynmat", shape="(6,6)").text = str(
        np.diag(0.5 / masses).ravel().tolist()
    )
    root = motion_input("normalmodes", options)
    root.find("system/beads/q").text = str([0.0] * 6)
    temp = root.find("system/ensemble/temperature")
    temp.attrib.clear()
    temp.text = ".2"
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        # Stored pair maps are a supported input. An uncoupled harmonic potential has exactly zero pair correction.
        for i in range(6):
            for j in range(6):
                if i != j:
                    with open(
                        f"{sim.output_maker.prefix}.vcoupled_grid.{i}.{j}.dat", "wb"
                    ) as f:
                        np.save(f, np.zeros((101, 101)))
        sim.run(steps=23, write_outputs=False)
        calc = sim.syslist[0].motion.calc
        expected = sim.syslist[0].motion.w[:, None] * (np.arange(12) + 0.5)
        np.testing.assert_allclose(calc.evals_vscf, expected, atol=1e-11)
        np.testing.assert_allclose(calc.rho_grids.sum(axis=1), 1, atol=1e-12)
        assert sim.syslist[0].motion.finished


def test_reactant_postprocessing_recovers_harmonic_frequencies(simulation):
    from ipi.utils.tools.instanton_postproc import instanton_compute, cm2au

    options = ET.Element("normalmodes", mode="imf")
    masses = np.repeat([2.0, 3.0], 3)
    ET.SubElement(options, "dynmat", shape="(6,6)").text = str(
        np.diag(0.5 / masses).ravel().tolist()
    )
    root = motion_input("normalmodes", options)
    with open("harmonic.xml", "w") as f:
        f.write(ET.tostring(root, encoding="unicode"))
    with pytest.raises(SystemExit) as stopped:
        instanton_compute("harmonic.xml", "reactant", 0.2, asr="none", nbeadsR=4)
    assert stopped.value.code == 0
    np.testing.assert_allclose(
        np.sort(np.loadtxt("eigenvalues_reactant.dat")),
        np.sort(0.5 / masses),
        atol=1e-12,
    )
    np.testing.assert_allclose(
        np.sort(np.loadtxt("freq.dat") * cm2au),
        np.sort(np.sqrt(0.5 / masses)),
        atol=1e-12,
    )


@pytest.mark.parametrize("update", ["sd", "ik", "nmik"])
def test_self_consistent_phonon_update_preserves_exact_harmonic_solution(
    simulation, update
):
    options = ET.Element("scp", mode="cl")
    for key, value in {
        "max_steps": 8,
        "max_iter": 3,
        "displace_mode": update,
        "precheck": "False",
    }.items():
        ET.SubElement(options, key).text = str(value)
    masses = np.repeat([2.0, 3.0], 3)
    ET.SubElement(options, "dynmat", shape="(6,6)").text = str(
        np.diag(0.5 / masses).ravel().tolist()
    )
    root = motion_input("scp", options)
    root.find("system/beads/q").text = str([0.0] * 6)
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        sim.run(steps=10, write_outputs=False)
        mover = sim.syslist[0].motion
        np.testing.assert_allclose(mover.beads.q, 0, atol=1e-12)
        np.testing.assert_allclose(mover.K, 0.5 * np.eye(6), atol=1e-12)
        force, error, weights = mover.phononator.weighted_force()
        np.testing.assert_allclose(force, 0, atol=1e-12)
        np.testing.assert_allclose(error, 0, atol=1e-12)
        np.testing.assert_allclose(weights, 1, atol=1e-12)
        np.testing.assert_allclose(
            mover.phononator.weighted_hessian(), 0.5 * np.eye(6), atol=1e-12
        )

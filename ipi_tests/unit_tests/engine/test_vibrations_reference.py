"""Vibrational analyses of a separable harmonic potential."""

import xml.etree.ElementTree as ET
from types import SimpleNamespace

import numpy as np
import pytest

from ipi.engine.motion.vscf import IMF, _coupling_on_grid
from ipi.utils.tools.instanton_postproc import cm2au, instanton_compute
from ipi_tests.unit_tests.engine.inputs import motion_input


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


@pytest.mark.parametrize("reuse_map,nparallel", [(False, 1), (True, 1), (False, 2)])
def test_vscf_pair_mapping_has_zero_harmonic_mode_coupling(
    simulation, reuse_map, nparallel
):
    options = ET.Element("normalmodes", mode="vscf")
    for key, value in {
        "nint": 51,
        "nbasis": 10,
        "nevib": 8,
        "grid": "False",
        "solve": "False",
        "nparallel": nparallel,
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
    """The free energies per primitive cell include the energy of the minimum
    per primitive cell once"""
    options = ET.Element("normalmodes", mode="vscf")
    for key, value in {
        "nint": 101,
        "nbasis": 12,
        "nevib": 20,
        "grid": "True",
        "solve": "True",
        "nprim": 2,
    }.items():
        ET.SubElement(options, key).text = str(value)
    masses = np.repeat([2.0, 3.0], 3)
    ET.SubElement(options, "dynmat", shape="(6,6)").text = str(
        np.diag(0.5 / masses).ravel().tolist()
    )
    root = motion_input("normalmodes", options)
    root.find("system/beads/q").text = str([0.0] * 6)
    # the energy of the minimum is 0.7
    ET.SubElement(root.find("ffdirect"), "offset").text = "-0.7"
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
        omega = np.sqrt(0.5 / masses)
        free = np.sum(omega / 2 + 0.2 * np.log1p(-np.exp(-omega / 0.2)))
        assert calc.total_imf_free_energy == pytest.approx((0.7 + free) / 2, abs=1e-8)
        assert calc.total_vscf_free_energy == pytest.approx((0.7 + free) / 2, abs=1e-8)


def test_reactant_postprocessing_recovers_harmonic_frequencies(simulation):

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


@pytest.mark.parametrize("mode", ["fd", "nmfd", "enmfd"])
def test_phonon_refinement_recovers_harmonic_force_constants(simulation, mode):
    options = ET.Element("vibrations", mode=mode)
    masses = np.repeat([2.0, 3.0], 3)
    ET.SubElement(options, "dynmat", shape="(6,6)").text = str(
        np.diag(0.8 / masses).ravel().tolist()
    )
    ET.SubElement(options, "pos_shift").text = ".001"
    ET.SubElement(options, "energy_shift").text = ".0001"
    root = motion_input("vibrations", options)
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        motion = sim.syslist[0].motion
        motion.exit_on_convergence = False
        sim.run(steps=7, write_outputs=False)
        np.testing.assert_allclose(
            motion.refdynmatrix, np.diag(0.5 / masses), atol=1e-10
        )


class SeparableMorse:
    """Driver for U = sum_k D (1 - exp(-a_k x_k))^2 + c x_i^2 x_j^2 over the
    Cartesian coordinates, with curvatures k = 2 D a_k^2 at the minimum."""

    def __init__(self, depth, curvatures, coupling=0.0, pair=(0, 1)):
        self.depth = depth
        self.a = np.sqrt(np.asarray(curvatures) / (2 * depth))
        self.coupling, self.pair = coupling, pair

    def __call__(self, cell, positions):
        x = positions.reshape(-1)
        e = np.exp(-self.a * x)
        energy = np.sum(self.depth * (1 - e) ** 2)
        force = -2 * self.depth * self.a * e * (1 - e)
        i, j = self.pair
        energy += self.coupling * x[i] ** 2 * x[j] ** 2
        force[i] -= 2 * self.coupling * x[i] * x[j] ** 2
        force[j] -= 2 * self.coupling * x[j] * x[i] ** 2
        return energy, force.reshape(positions.shape), np.zeros((3, 3)), ""


MASSES = np.repeat([2.0, 3.0], 3)
CURVATURES = np.array([0.5, 0.62, 0.74, 0.86, 0.98, 1.1])


def morse_levels(omega, depth):
    """Bound levels of a Morse oscillator of frequency omega and depth D (hbar = 1)."""
    n = np.arange(int(2 * depth / omega - 0.5) + 1)
    return omega * (n + 0.5) - (omega * (n + 0.5)) ** 2 / (4 * depth)


def thermal(levels, temperature):
    """Free and internal energy of a set of levels."""
    weights = np.exp(-(levels - levels[0]) / temperature)
    free = levels[0] - temperature * np.log(weights.sum())
    return free, np.sum(levels * weights) / weights.sum()


def coupled_morse_levels(depth, k, m, coupling, npts=160, nstates=30):
    """Levels of two Morse oscillators (curvatures k, masses m) coupled by
    c x_1^2 x_2^2, diagonalised in the product of their lowest eigenstates,
    computed with a sinc-DVR in mass-weighted coordinates."""
    states, x2 = [], []
    for ki, mi in zip(k, m):
        a = np.sqrt(ki / (2 * depth * mi))
        q = np.linspace(-3 / a / 2, 12 / a / 2, npts)
        dq = q[1] - q[0]
        n = np.arange(npts)
        diff = n[:, None] - n[None, :]
        kinetic = np.where(
            diff == 0, np.pi**2 / 3, 2.0 * (-1.0) ** diff / np.where(diff, diff, 1) ** 2
        ) / (2 * dq**2)
        levels, vectors = np.linalg.eigh(
            kinetic + np.diag(depth * (1 - np.exp(-a * q)) ** 2)
        )
        vectors = vectors[:, :nstates]
        states.append(levels[:nstates])
        x2.append(vectors.T @ np.diag(q**2 / mi) @ vectors)
    hamiltonian = np.diag(np.add.outer(*states).ravel()) + coupling * np.kron(*x2)
    return np.linalg.eigvalsh(hamiltonian)


def morse_input(mode, temperature, options):
    root = ET.Element("normalmodes", mode=mode)
    for key, value in options.items():
        ET.SubElement(root, key).text = str(value)
    ET.SubElement(root, "dynmat", shape="(6,6)").text = str(
        np.diag(CURVATURES / MASSES).ravel().tolist()
    )
    root = motion_input("normalmodes", root)
    root.find("system/beads/q").text = str([0.0] * 6)
    temp = root.find("system/ensemble/temperature")
    temp.attrib.clear()
    temp.text = str(temperature)
    return ET.tostring(root, encoding="unicode")


@pytest.mark.parametrize("temperature", [0.2, 2e-4])
def test_independent_modes_are_exact_for_separable_morse_potentials(
    simulation, temperature
):
    """IMF is exact for a potential that is a sum of Morse potentials along the
    normal modes, also at temperatures at which the Boltzmann factors underflow"""

    depth = 20.0
    options = {"nint": 301, "nbasis": 16, "nevib": 40, "athresh": 1e-9}
    with simulation(morse_input("imf", temperature, options)) as sim:
        sim.syslist[0].forces.ff["harmonic"].driver = SeparableMorse(depth, CURVATURES)
        sim.run(steps=7, write_outputs=False)
        calc = sim.syslist[0].motion.calc
        exact = np.sum(
            [
                thermal(morse_levels(omega, depth), temperature)
                for omega in np.sqrt(CURVATURES / MASSES)
            ],
            axis=0,
        )
        assert calc.total_anhar_free_energy == pytest.approx(exact[0], abs=2e-6)
        assert calc.total_anhar_internal_energy == pytest.approx(exact[1], abs=2e-6)


@pytest.mark.parametrize("alpha", [1.0, 0.5])
def test_vscf_free_energy_of_coupled_morse_modes(simulation, alpha):
    """The VSCF correction to the free energy of two coupled Morse modes is an
    upper bound to the exact one, larger by a term of second order in the
    coupling (it would be twice as large if the pair coupling were counted twice).
    Only the modes in modes.dat are treated, and the mixing of the mean field
    (alpha) does not change the result"""

    depth, coupling, temperature = 40.0, 1e-3, 0.5
    # coordinates 1 and 4 are modes 2 and 3; mode 0 (coordinate 0) is left out
    pair = (1, 4)
    options = {
        "nint": 201,
        "nbasis": 30,
        "nevib": 40,
        "athresh": 1e-9,
        "alpha": alpha,
        "grid": "True",
        "solve": "True",
    }
    # the modes to treat are read when the motion is bound
    np.savetxt("i-pi.modes.dat", [[1, 2, 3, 4, 5]], fmt="%d")
    with simulation(morse_input("vscf", temperature, options)) as sim:
        sim.syslist[0].forces.ff["harmonic"].driver = SeparableMorse(
            depth, CURVATURES, coupling, pair
        )
        sim.run(steps=200, write_outputs=False)
        assert sim.syslist[0].motion.finished
        calc = sim.syslist[0].motion.calc
        assert calc.inms == [1, 2, 3, 4, 5]
        omega = np.sqrt(CURVATURES / MASSES)
        uncoupled = [
            thermal(morse_levels(omega[k], depth), temperature)[0] for k in pair
        ]
        coupled = coupled_morse_levels(
            depth, CURVATURES[list(pair)], MASSES[list(pair)], coupling
        )
        exact = thermal(coupled, temperature)[0] - sum(uncoupled)
        # the corrections to the uncoupled free energy are compared, as the VSCF
        # solver uses a fixed basis for the one-dimensional problems
        correction = calc.total_vscf_free_energy - calc.total_imf_free_energy
        assert correction - exact > -2e-7
        assert correction - exact < 1e-2 * exact


def test_pair_coupling_of_a_polynomial_potential():
    """The pair coupling V(qi, qj) - V(qi, 0) - V(0, qj) + V(0, 0) is extracted
    exactly from a cubic polynomial sampled on a grid, indexed as [j, i]"""

    qi, qj = np.linspace(-1.2, 1.5, 10), np.linspace(-0.8, 1.1, 8)
    x, y = qi[:, None], qj[None, :]
    v = 0.2 + 1.3 * x**2 + 0.7 * y**2 + 0.4 * x * y + 0.05 * x**3 - 0.1 * x * y**2
    igrid, jgrid = np.linspace(-1.2, 1.5, 7), np.linspace(-0.8, 1.1, 5)
    expected = 0.4 * np.outer(jgrid, igrid) - 0.1 * np.outer(jgrid**2, igrid)
    np.testing.assert_allclose(
        _coupling_on_grid(qi, qj, v.ravel(), igrid, jgrid), expected, atol=1e-12
    )

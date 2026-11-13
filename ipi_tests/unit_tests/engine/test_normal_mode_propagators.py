"""Mass-scaled free normal modes must follow the specified harmonic dynamics."""

import xml.etree.ElementTree as ET

import numpy as np
import pytest
from scipy.linalg import expm
from ipi.utils.array_backend import to_numpy, xp
from ipi_tests.unit_tests.engine.inputs import harmonic_input

MODE_SETTINGS = [
    ("rpmd", []),
    ("manual", [0.02, 0.03, 0.04]),
    ("pa-cmd", [0.02]),
    ("wmax-cmd", [0.1, 0.02]),
]

# Cayley propagation is only checked with physical (RPMD) masses.
PROPAGATOR_CASES = [
    (style, freqs, propagator)
    for style, freqs in MODE_SETTINGS
    for propagator in ["exact", "cayley"]
    if propagator == "exact" or style == "rpmd"
]


# Four-bead unit-spring matrices have these eigenvalues: a cyclic chain in
# Fourier mode order, and a chain with free ends in ascending order.
CLOSED_SPRING = np.sqrt([0.0, 2.0, 4.0, 2.0])
OPEN_SPRING = np.sqrt([0.0, 2 - np.sqrt(2), 2.0, 2 + np.sqrt(2)])


def mass_factors(style, freqs, spring):
    """The mass factors of the modes of frequencies `spring`."""
    factors = np.ones(4)
    if style == "manual":
        factors[1:] = (spring[1:] / freqs) ** 2
    elif style == "pa-cmd":
        factors[1:] = (spring[1:] / freqs[0]) ** 2
    elif style == "wmax-cmd":
        # The first free internal mode must oscillate at wtarget. Adding
        # physical curvature wmax**2 must make all internal modes oscillate at
        # the same frequency. Those two conditions determine each mass factor.
        wmax, wtarget = freqs
        common_frequency = wtarget * np.sqrt(1 + (wmax / spring[1]) ** 2)
        factors[1:] = (spring[1:] ** 2 + wmax**2) / common_frequency**2
    return factors


def normal_mode_input(style, freqs, propagator):
    root = ET.fromstring(harmonic_input(nbeads=4, timestep=0.1))
    # In atomic units P*kB*T/hbar = 4*0.2.
    temperature = root.find("system/ensemble/temperature")
    temperature.attrib.clear()
    temperature.text = "0.2"
    normal = ET.SubElement(
        root.find("system"), "normal_modes", transform="matrix", propagator=propagator
    )
    ET.SubElement(normal, "frequencies", style=style).text = str(freqs)
    return ET.tostring(root, encoding="unicode")


@pytest.mark.parametrize("style,freqs,propagator", PROPAGATOR_CASES)
def test_frequencies_and_mass_factors(simulation, style, freqs, propagator):
    """Spring eigenvalues and mass factors of closed and open paths."""
    with simulation(normal_mode_input(style, freqs, propagator)) as sim:
        nm = sim.syslist[0].nm
        assert float(nm.omegan) == pytest.approx(0.8, rel=1e-12)
        for spring, factors, expected_spring in [
            (nm.omegak, nm.nm_factor, 0.8 * CLOSED_SPRING),
            (nm.o_omegak, nm.o_nm_factor, 0.8 * OPEN_SPRING),
        ]:
            np.testing.assert_allclose(
                to_numpy(spring), expected_spring, rtol=1e-12, atol=1e-14
            )
            np.testing.assert_allclose(
                to_numpy(factors),
                mass_factors(style, freqs, expected_spring),
                rtol=1e-12,
                atol=1e-14,
            )


@pytest.mark.parametrize("style,freqs", MODE_SETTINGS)
def test_mass_factors_follow_the_temperature(simulation, style, freqs):
    """The mass factors of closed and open paths are recomputed when the
    temperature changes the frequencies of the free ring polymer"""
    with simulation(normal_mode_input(style, freqs, "exact")) as sim:
        system = sim.syslist[0]
        nm = system.nm
        # computes the factors at the initial temperature
        for factors in (nm.nm_factor_, nm.o_nm_factor_):
            factors.get()
        system.ensemble.temp = 0.4
        assert float(nm.omegan) == pytest.approx(1.6, rel=1e-12)
        for factors, spring in [
            (nm.nm_factor, 1.6 * CLOSED_SPRING),
            (nm.o_nm_factor, 1.6 * OPEN_SPRING),
        ]:
            np.testing.assert_allclose(
                to_numpy(factors),
                mass_factors(style, freqs, spring),
                rtol=1e-12,
                atol=1e-14,
            )


@pytest.mark.parametrize("style,freqs,propagator", PROPAGATOR_CASES)
@pytest.mark.parametrize("open_paths", [False, True], ids=["closed", "open"])
def test_centroid_positions_and_momenta_follow_free_motion(
    simulation, style, freqs, propagator, open_paths
):
    """A free centroid drifts with its physical mass and keeps its momentum."""
    root = ET.fromstring(normal_mode_input(style, freqs, propagator))
    # Remove physical forces; the internal ring-polymer springs remain.
    root.find("ffdirect/parameters").text = "{k1: 0.0}"
    if open_paths:
        ET.SubElement(root.find("system/normal_modes"), "open_paths").text = "[0,1]"
    masses = np.repeat([2.0, 3.0], 3)
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        system = sim.syslist[0]
        initial_position = to_numpy(system.beads.q).mean(axis=0)
        initial_momentum = to_numpy(system.beads.p).mean(axis=0)
        # The input timestep is 0.1; four-bead mode zero is twice the bead mean.
        expected_position = initial_position + 0.1 * initial_momentum / masses
        sim.run(steps=1, write_outputs=False)
        np.testing.assert_allclose(
            to_numpy(system.beads.q).mean(axis=0),
            expected_position,
            rtol=1e-12,
            atol=1e-14,
        )
        np.testing.assert_allclose(
            to_numpy(system.beads.p).mean(axis=0),
            initial_momentum,
            rtol=1e-12,
            atol=1e-14,
        )
        np.testing.assert_allclose(
            to_numpy(system.nm.qnm[0]), 2 * expected_position, rtol=1e-12, atol=1e-14
        )
        np.testing.assert_allclose(
            to_numpy(system.nm.pnm[0]), 2 * initial_momentum, rtol=1e-12, atol=1e-14
        )


@pytest.mark.parametrize("paths", ["closed", "open_paths", "bosons"])
def test_spring_force_is_negative_energy_gradient(simulation, paths):
    """The spring force is minus the gradient of the spring energy for every
    bead, and follows changes of the masses"""
    root = ET.fromstring(harmonic_input(nbeads=4))
    if paths != "closed":
        normal = ET.SubElement(root.find("system"), "normal_modes")
        ET.SubElement(normal, paths).text = "[0]"
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        beads, nm = sim.syslist[0].beads, sim.syslist[0].nm
        original = np.random.default_rng(82).normal(size=(4, 6))
        beads.q_[:] = xp.asarray(original)
        force = to_numpy(nm.fspring).copy()
        numerical = np.zeros_like(original)
        step = 1e-5
        for index in np.ndindex(original.shape):
            beads.q_[index] = float(original[index] + step)
            forward = float(nm.vspring)
            beads.q_[index] = float(original[index] - step)
            backward = float(nm.vspring)
            beads.q_[index] = float(original[index])
            numerical[index] = -(forward - backward) / (2 * step)
        np.testing.assert_allclose(force, numerical, rtol=1e-7, atol=1e-13)
        np.testing.assert_allclose(force.sum(axis=0), 0, atol=1e-15)
        beads.m_[:] *= 2
        np.testing.assert_allclose(to_numpy(nm.fspring), 2 * force, rtol=1e-12)


@pytest.mark.parametrize("style,freqs,propagator", PROPAGATOR_CASES)
def test_open_and_closed_mode_propagators_match_harmonic_matrix(
    simulation, style, freqs, propagator
):
    with simulation(normal_mode_input(style, freqs, propagator)) as sim:
        nm = sim.syslist[0].nm
        dt = float(nm.dt)
        for arrays in [
            (nm.omegak, nm.nm_factor, nm.prop_pq),
            (nm.o_omegak, nm.o_nm_factor, nm.o_prop_pq),
        ]:
            omega, factors, prop = map(to_numpy, arrays)
            for k in range(1, 4):
                generator = np.array([[0.0, -omega[k] ** 2], [1 / factors[k], 0.0]])
                if propagator == "exact":
                    expected = expm(dt * generator)
                else:
                    half = dt * generator / 2
                    expected = np.linalg.solve(np.eye(2) - half, np.eye(2) + half)
                np.testing.assert_allclose(prop[k], expected, atol=1e-12)
                assert np.linalg.det(prop[k]) == pytest.approx(1.0)
        cloned = nm.copy()
        assert cloned.mode == nm.mode
        np.testing.assert_array_equal(cloned.nm_freqs, nm.nm_freqs)

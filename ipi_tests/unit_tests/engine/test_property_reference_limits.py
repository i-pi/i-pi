"""Isotope and kinetic estimators in independently solvable limits."""

import numpy as np
import pytest
from ipi_tests.unit_tests.engine.test_simulation_invariants import harmonic_input


@pytest.mark.parametrize("method", ["get_isotope_zetatd", "get_isotope_zetatd_4th"])
@pytest.mark.parametrize("alpha", [0.5, 1.0, 2.0])
@pytest.mark.parametrize("atom,indices", [("", [0, 1]), ("H", [0]), ("1", [1])])
def test_free_particle_isotope_weight_equals_spring_boltzmann_ratio(
    simulation, method, alpha, atom, indices
):
    xml = harmonic_input(nbeads=4).replace("{k1: 0.5}", "{k1: 0.0}")
    with simulation(xml) as sim:
        s = sim.syslist[0]
        q = s.beads.q.copy().reshape(4, 2, 3)
        # The discretized free-particle action is a quadratic form in adjacent displacements.
        differences = q - np.roll(q, 1, axis=0)
        temp = s.ensemble.temp
        action = (
            0.5
            * np.asarray([2.0, 3.0])
            * (4 * temp) ** 2
            * np.sum(differences**2, axis=(0, 2))
        )
        delta = (alpha - 1) * action[indices]
        weight = alpha**4.5 * np.mean(np.exp(-delta / (4 * temp)))
        expected = [delta.mean(), np.mean(delta**2), weight]
        if method.endswith("4th"):
            expected += [weight, weight]
        np.testing.assert_allclose(
            getattr(s.properties, method)(str(alpha), atom),
            expected,
            rtol=1e-12,
            atol=1e-15,
        )


@pytest.mark.parametrize("method", ["get_isotope_zetasc", "get_isotope_zetasc_4th"])
@pytest.mark.parametrize("atom", ["", "H", "1"])
def test_unchanged_isotope_has_unit_weight(simulation, method, atom):
    with simulation(harmonic_input(nbeads=4)) as sim:
        s = sim.syslist[0]
        original = s.beads.q.copy()
        actual = getattr(s.properties, method)("1", atom)
        np.testing.assert_allclose(actual, [0, 0] + [1] * (len(actual) - 2), atol=1e-13)
        np.testing.assert_array_equal(s.beads.q, original)


@pytest.mark.parametrize("alpha", [0.5, 2.0])
def test_scaled_isotope_harmonic_weight_from_variance(simulation, alpha):
    with simulation(harmonic_input(nbeads=4)) as sim:
        s = sim.syslist[0]
        q = s.beads.q.copy()
        fluctuations = q - q.mean(axis=0)
        delta = 0.25 * (1 / alpha - 1) * np.sum(fluctuations[:, :3] ** 2)
        expected = [delta, delta**2, np.exp(-delta / (4 * s.ensemble.temp))]
        np.testing.assert_allclose(
            s.properties.get_isotope_zetasc(str(alpha), "H"), expected, atol=1e-13
        )
        np.testing.assert_array_equal(s.beads.q, q)


@pytest.mark.parametrize(
    "method",
    [
        "get_isotope_zetatd",
        "get_isotope_zetatd_4th",
        "get_isotope_zetasc",
        "get_isotope_zetasc_4th",
        "get_rg",
        "get_ktens",
    ],
)
@pytest.mark.parametrize("atom", ["7", "missing"])
def test_atom_selectors_reject_absent_atoms(simulation, method, atom):
    with simulation(harmonic_input(nbeads=4)) as sim:
        prop = getattr(sim.syslist[0].properties, method)
        if method == "get_ktens" and atom == "missing":
            np.testing.assert_array_equal(prop(atom=atom), np.zeros(6))
        else:
            with pytest.raises(IndexError):
                prop(atom=atom)


@pytest.mark.parametrize("atom,indices", [("", [0, 1]), ("H", [0]), ("1", [1])])
def test_gyration_and_kinetic_tensor_of_harmonic_path(simulation, atom, indices):
    with simulation(harmonic_input(nbeads=4)) as sim:
        s = sim.syslist[0]
        q = s.beads.q.copy().reshape(4, 2, 3)
        centered = q - q.mean(axis=0)
        radii = np.sqrt(np.sum(centered**2, axis=(0, 2)) / 4)
        assert s.properties.get_rg(atom) == pytest.approx(radii[indices].mean())
        tensor = sum(
            0.5 * s.ensemble.temp * np.eye(3)
            + 0.25 * centered[:, i].T @ centered[:, i] / 4
            for i in indices
        )
        np.testing.assert_allclose(
            s.properties.get_ktens(atom),
            tensor[[0, 1, 2, 0, 0, 1], [0, 1, 2, 1, 2, 2]],
            atol=1e-13,
        )


@pytest.mark.parametrize(
    "method", ["get_yama_estimators", "get_scyama_estimators", "get_scop_estimators"]
)
def test_free_particle_energy_and_temperature_derivative(simulation, method):
    with simulation(harmonic_input(nbeads=4).replace("{k1: 0.5}", "{k1: 0.0}")) as sim:
        s = sim.syslist[0]
        t = s.ensemble.temp
        value = getattr(s.properties, method)(fd_delta=0.001)
        sign = 1 if method == "get_scop_estimators" else -1
        np.testing.assert_allclose(value, [3 * t, sign * 3 * t * t], atol=1e-13)


@pytest.mark.parametrize("atom,mass", [("H", 2), ("O", 3)])
def test_free_particle_end_to_end_distribution(simulation, atom, mass):
    with simulation(harmonic_input(nbeads=4).replace("{k1: 0.5}", "{k1: 0.0}")) as sim:
        s = sim.syslist[0]
        q = s.beads.q.copy()
        actual = s.properties.get_linlin(".2", "-.3", ".1", atom)
        assert actual == pytest.approx(
            np.exp(-float(mass) * (0.04 + 0.09 + 0.01) * s.ensemble.temp / 2)
        )
        np.testing.assert_array_equal(s.beads.q, q)
        with pytest.raises(IndexError):
            s.properties.get_linlin(atom="absent")


def test_harmonic_takahashi_imada_correction_has_force_squared_scaling(simulation):
    with simulation(harmonic_input(nbeads=4)) as sim:
        s = sim.syslist[0]
        q = s.beads.q.copy()
        temp = s.ensemble.temp
        expected = 0.25 * np.sum(q * q / s.beads.m3) / (24 * (4 * temp) ** 2 * 4)
        assert s.properties.get_ti_term() == pytest.approx(expected)
        correction = s.properties.get_ti_correction()
        exponent = -expected / temp
        np.testing.assert_allclose(
            correction, [exponent, exponent**2, np.exp(exponent)]
        )
        assert s.properties.get_ti_term("H") + s.properties.get_ti_term(
            "O"
        ) == pytest.approx(expected)
        assert s.properties.get_ti_term("absent") == 0
        with pytest.raises(IndexError):
            s.properties.get_ti_term("7")


def test_free_particle_fourth_order_correction_is_zero(simulation):
    with simulation(harmonic_input(nbeads=4).replace("{k1: 0.5}", "{k1: 0.0}")) as sim:
        np.testing.assert_array_equal(
            sim.syslist[0].properties.get_chin_correction(), [0, 0, 1]
        )

"""Path-integral estimators checked against explicit expressions.

Most tests use the quartic `reference_system` (conftest.py) and compare with
formulas evaluated in mathematical_references; the others use the harmonic
input, free-particle limits and hand-built paths.
"""

import xml.etree.ElementTree as ET

import numpy as np
import pytest

from ipi.engine.simulation import Simulation
from ipi.utils.units import Constants, unit_to_internal
from ipi_tests.unit_tests.engine.inputs import harmonic_input
from ipi_tests.unit_tests.engine.mathematical_references import (
    fourier_matrix,
    mode_mass_factors,
    permutation_reference,
    polynomial,
    scaled_derivatives,
    sc_reference,
    sc_scaled_derivatives,
    sc_weights,
    spring_action,
)

TEMPERATURE_300K = unit_to_internal("temperature", "kelvin", 300)


def harmonic_system(nbeads=4, spring=0.5, temperature=None):
    """Read-only harmonic system; `temperature` in atomic units overrides 300 K."""
    root = ET.fromstring(
        harmonic_input(nbeads=nbeads).replace("{k1: 0.5}", f"{{k1: {spring}}}")
    )
    if temperature is not None:
        node = root.find("system/ensemble/temperature")
        node.attrib.clear()
        node.text = str(temperature)
    xml = ET.tostring(root, encoding="unicode")
    return Simulation.load_from_xml(xml, read_only=True).syslist[0]


def vector(tensor):
    return tensor[[0, 1, 2, 0, 0, 1], [0, 1, 2, 1, 2, 2]]


@pytest.mark.parametrize("count", [2, 4, 6])
def test_spring_primitive_and_conserved_energies(reference_system, count):
    s, r = reference_system(count)
    action = spring_action(r.q, r.masses, r.temperature)
    force = (
        -((count * r.temperature) ** 2)
        * np.repeat(r.masses, 3)
        * (2 * r.q - np.roll(r.q, 1, axis=0) - np.roll(r.q, -1, axis=0))
    )
    assert s.properties["spring"][0] == pytest.approx(action / count)
    assert s.properties["kinetic_td"][0] == pytest.approx(
        3 * count * r.temperature - action / count
    )
    np.testing.assert_allclose(s.trajs["forces_spring"][0], force, atol=1e-13)
    energy = polynomial(r.q, r.a, r.b)[0].sum()
    kinetic = np.sum(r.p**2 / (2 * np.repeat(r.masses, 3)))
    assert s.properties["conserved"][0] == pytest.approx(
        (action + energy + kinetic) / count
    )


@pytest.mark.parametrize(
    "style,frequencies",
    [
        ("rpmd", []),
        ("manual", [0.3, 0.5, 0.7]),
        ("pa-cmd", [0.4]),
        ("wmax-cmd", [0.8, 0.3]),
    ],
)
@pytest.mark.parametrize("atom", ["", "0", "1"])
def test_kinetic_md_uses_mode_masses(reference_system, style, frequencies, atom):
    s, r = reference_system(4, style=style, frequencies=frequencies)
    modes = fourier_matrix(4) @ r.p
    factors = mode_mass_factors(4, 0.4, style, frequencies)
    mass = factors[:, None] * np.repeat(r.masses, 3)
    selected = (
        np.arange(6) if atom == "" else np.arange(3 * int(atom), 3 * int(atom) + 3)
    )
    expected = np.sum(modes[:, selected] ** 2 / (2 * mass[:, selected])) / 4
    assert s.properties.get_kinmd(atom) == pytest.approx(expected)
    for mode in range(4):
        assert s.properties.get_kinmd(atom, nm=str(mode)) == pytest.approx(
            np.sum(modes[mode, selected] ** 2 / (2 * mass[mode, selected]))
        )
    stress = np.zeros((3, 3))
    for mode in range(4):
        vectors = modes[mode].reshape(2, 3)
        for i, m in enumerate(r.masses):
            stress += np.outer(vectors[i], vectors[i]) / (m * factors[mode])
    np.testing.assert_allclose(s.properties.get_kstress_md(), stress, atol=1e-13)
    if style != "rpmd":
        with pytest.raises(ValueError, match="mass matrix"):
            s.properties.get_kinmd(bead="0")


@pytest.mark.parametrize("nbeads", [1, 4])
@pytest.mark.parametrize(
    "selection,indices", [("", [0, 1]), ("H", [0]), ("O", [1]), ("0", [0]), ("1", [1])]
)
def test_kinetic_estimators_select_atoms(nbeads, selection, indices):
    s = harmonic_system(nbeads)
    # p_H=(.1,.2,-.1), m_H=2; p_O=(.3,-.2,.1), m_O=3.
    expected_kinetic = sum([0.06 / 4, 0.14 / 6][index] for index in indices)
    kinetic, count = s.properties.get_kinmd(atom=selection, return_count=True)
    assert count == len(indices)
    assert kinetic == pytest.approx(expected_kinetic)
    assert s.properties.get_kinmd(atom=selection, bead="0") == pytest.approx(
        expected_kinetic
    )
    # Only H is spread: Var(0,.01,...,(P-1)*.01) = .01²(P²-1)/12.
    virial = 0.25 * 0.01**2 * (nbeads**2 - 1) / 12 if 0 in indices else 0
    expected_cv = len(indices) * 1.5 * Constants.kb * TEMPERATURE_300K + virial
    assert s.properties.get_kincv(selection) == pytest.approx(expected_cv, abs=1e-12)


@pytest.mark.parametrize(
    "constraint",
    [
        "none",
        "center-of-mass",
        "fixed-coordinate",
        "centroid",
        "centroid-and-coordinate",
    ],
)
@pytest.mark.parametrize("atom", ["", "0", "1"])
@pytest.mark.parametrize("mode", ["", "0", "1"])
def test_temperature_counts_constrained_degrees_of_freedom(
    reference_system, constraint, atom, mode, request
):
    if constraint == "centroid-and-coordinate" and mode != "":
        request.node.add_marker(
            pytest.mark.xfail(
                strict=True,
                raises=AssertionError,
                reason="BUG: mode-selected temperature miscounts fixed coordinates when the centroid is fixed",
            )
        )
    s, r = reference_system(4)
    mode_p = fourier_matrix(4) @ r.p
    selected_atoms = [0, 1] if atom == "" else [int(atom)]
    selected = np.concatenate([np.arange(3 * i, 3 * i + 3) for i in selected_atoms])
    removed = np.zeros((4, 6))
    s.motion.fixcom = constraint == "center-of-mass"
    s.motion.fixatoms_dof = (
        np.array([1, 4], dtype=int)
        if "coordinate" in constraint
        else np.array([], dtype=int)
    )
    s.motion.enstype = "nvt-cc" if "centroid" in constraint else "nve"
    if constraint == "center-of-mass":
        p = mode_p[0].reshape(2, 3)
        p -= r.masses[:, None] / r.masses.sum() * p.sum(axis=0)
        # Conditional Gaussian covariance removes P*T*m_i/M per component.
        removed[0] = 4 * np.repeat(r.masses / r.masses.sum(), 3)
    if "centroid" in constraint:
        mode_p[0] = 0
        removed[0] = 4
    if "coordinate" in constraint:
        mode_p[:, [1, 4]] = 0
        removed[:, [1, 4]] = 4
    s.beads.p_[:] = fourier_matrix(4).T @ mode_p
    modes = range(4) if mode == "" else [int(mode)]
    energy = sum(
        np.sum(mode_p[k, selected] ** 2 / (2 * np.repeat(r.masses, 3)[selected]))
        for k in modes
    )
    missing = sum(np.sum(removed[k, selected]) for k in modes)
    divisor = 4 if mode == "" else 1
    expected = (2 * energy / divisor + r.temperature * missing / divisor) / (
        3 * len(selected_atoms) * 4
    )
    assert s.properties.get_temp(atom, nm=mode) == pytest.approx(
        expected, rel=1e-12, abs=1e-12
    )


@pytest.mark.parametrize("alpha", [0.0, 0.3])
@pytest.mark.parametrize("key", ["kstress_md", "kstress_cv", "kstress_tdsc"])
def test_kinetic_stress_estimators(reference_system, alpha, key):
    s, r = reference_system(alpha=alpha)
    ref = sc_reference(r.q, r.masses, r.a, r.b, alpha, r.temperature)
    delta = r.q - r.q.mean(axis=0)
    pc = r.p.mean(axis=0).reshape(2, 3)
    # The CV stress estimator replaces off-diagonal centroid momentum products
    # by their Maxwell expectation, zero, and keeps the diagonal ones.
    centroid = np.diag(sum(p**2 / mass for p, mass in zip(pc, r.masses))) * r.count

    def kinetic_stress(force):
        return np.triu(
            -sum(x.reshape(2, 3).T @ f.reshape(2, 3) for x, f in zip(delta, force))
            + centroid
        )

    cv = kinetic_stress(ref["force"])
    sc = kinetic_stress(ref["force"] + ref["force1"] + ref["force2"])
    md = sum(
        p.reshape(2, 3).T @ (p / np.repeat(r.masses, 3)).reshape(2, 3) for p in r.p
    )
    expected = dict(
        kstress_md=vector(md / 120),
        kstress_cv=vector(cv / (120 * r.count)),
        kstress_tdsc=np.trace(sc) / (360 * r.count),
    )
    np.testing.assert_allclose(
        s.properties[key][0], expected[key], rtol=1e-8, atol=1e-10
    )
    if key == "kstress_tdsc":
        np.testing.assert_allclose(
            s.properties.kstress_sctd(), sc, rtol=1e-8, atol=1e-10
        )


@pytest.mark.parametrize("first,second", [(0, 1), (1, 0), (0, 0), (1, 1)])
def test_cross_atom_kinetic_tensor(reference_system, first, second):
    s, r = reference_system()
    force = polynomial(r.q, r.a, r.b)[1].reshape(r.count, 2, 3)
    delta = (r.q - r.q.mean(axis=0)).reshape(r.count, 2, 3)
    tensor = -(
        r.masses[first] * delta[:, first].T @ force[:, second]
        + r.masses[second] * force[:, first].T @ delta[:, second]
    ) / (4 * r.count * np.sqrt(r.masses[first] * r.masses[second]))
    if first == second:
        tensor += np.eye(3) * r.temperature / 2
    np.testing.assert_allclose(
        s.properties.get_kij(str(first), str(second)), vector(tensor), atol=1e-13
    )


@pytest.mark.parametrize("atom,indices", [("", [0, 1]), ("H", [0]), ("1", [1])])
def test_gyration_and_kinetic_tensors_from_path_covariance(atom, indices):
    s = harmonic_system()
    q = np.arange(24).reshape(4, 6) / 20
    s.beads.q_[:] = q
    centered = (q - q.mean(axis=0)).reshape(4, 2, 3)
    covariance = np.einsum("bai,baj->aij", centered, centered) / 4
    variances = np.diagonal(covariance, axis1=1, axis2=2)
    temp = s.ensemble.temp
    # For U = k|q|^2/2 with k = 1/2, the CV kinetic tensor is T/2 + k Cov/2.
    assert s.properties.get_rg(atom) == pytest.approx(
        np.sqrt(variances.sum(axis=1))[indices].mean()
    )
    tensor = sum(0.5 * temp * np.eye(3) + 0.25 * covariance[i] for i in indices)
    np.testing.assert_allclose(s.properties.get_ktens(atom), vector(tensor))
    np.testing.assert_allclose(s.trajs.get_rg(), np.sqrt(variances).ravel())
    np.testing.assert_allclose(
        s.trajs.get_akcv(), 0.5 * temp + 0.25 * variances.ravel()
    )
    np.testing.assert_allclose(
        s.trajs.get_akcv_od().reshape(2, 3), 0.25 * covariance[:, [0, 0, 1], [1, 2, 2]]
    )


@pytest.mark.parametrize("alpha", [0.0, 0.3, 1.0])
def test_suzuki_chin_estimators(reference_system, alpha):
    s, r = reference_system(alpha=alpha)
    ref = sc_reference(r.q, r.masses, r.a, r.b, alpha, r.temperature)
    w, c = sc_weights(r.count, alpha, r.temperature)
    correction_force = ref["force1"] + ref["force2"]
    delta = r.q - r.q.mean(axis=0)
    cv = (
        3 * r.temperature
        - np.sum(delta * ref["force"] * np.tile([1, 0], r.count // 2)[:, None])
        / r.count
    )
    assert s.properties["potential_opsc"][0] == pytest.approx(ref["energy"][::2].mean())
    assert s.properties["potential_tdsc"][0] == pytest.approx(
        np.mean(w * ref["energy"] + 2 * c * ref["square"])
    )
    assert s.properties["kinetic_prsc"][0] == pytest.approx(
        3 * r.count * r.temperature
        - spring_action(r.q, r.masses, r.temperature) / r.count
        + np.mean(c * ref["square"])
    )
    assert s.properties["kinetic_opsc"][0] == pytest.approx(cv)
    expected_td = (
        3 * r.temperature
        - np.sum(delta * (ref["force"] + correction_force)) / (2 * r.count)
        + np.mean(c * ref["square"])
    )
    assert s.properties["kinetic_tdsc"][0] == pytest.approx(
        expected_td, rel=1e-8, abs=1e-10
    )
    np.testing.assert_allclose(
        s.trajs["forces_sc"][0],
        ref["force"] + correction_force,
        rtol=1e-8,
        atol=1e-10,
    )


def test_takahashi_imada_correction_for_harmonic_forces():
    s = harmonic_system()
    q = s.beads.q.copy()
    temp = s.ensemble.temp
    # |F|^2/m summed over beads, divided by 24 (P omega_P)^2 and by P.
    expected = 0.25 * np.sum(q * q / s.beads.m3) / (24 * (4 * temp) ** 2 * 4)
    assert s.properties.get_ti_term() == pytest.approx(expected)
    exponent = -expected / temp
    np.testing.assert_allclose(
        s.properties.get_ti_correction(), [exponent, exponent**2, np.exp(exponent)]
    )
    assert s.properties.get_ti_term("H") + s.properties.get_ti_term(
        "O"
    ) == pytest.approx(expected)


@pytest.mark.parametrize("atom", [0, 1])
def test_displaced_path_estimator(reference_system, atom):
    s, r = reference_system()
    opening = np.linspace(
        -(r.count - 1) / (2 * r.count), (r.count - 1) / (2 * r.count), r.count
    )
    u = np.array([0.17, -0.23, 0.11])
    shifted = r.q.copy()
    shifted[:, 3 * atom : 3 * atom + 3] += opening[:, None] * u
    change = (polynomial(shifted, r.a, r.b)[0] - polynomial(r.q, r.a, r.b)[0]).sum()
    expected = np.exp(
        -r.masses[atom] * np.dot(u, u) * r.temperature / 2
        - change / (r.count * r.temperature)
    )
    np.testing.assert_allclose(
        [s.properties.opening(i) for i in range(r.count)], opening
    )
    assert s.properties[f"displacedpath(.17;-.23;.11;{atom})"][0] == pytest.approx(
        expected
    )


@pytest.mark.parametrize("atom,mass", [("H", 2), ("O", 3)])
def test_free_particle_end_to_end_distribution(atom, mass):
    s = harmonic_system(spring=0.0)
    q = s.beads.q.copy()
    actual = s.properties.get_linlin(".2", "-.3", ".1", atom)
    assert actual == pytest.approx(
        np.exp(-float(mass) * (0.04 + 0.09 + 0.01) * s.ensemble.temp / 2)
    )
    np.testing.assert_array_equal(s.beads.q, q)


@pytest.mark.parametrize(
    "method", ["get_yama_estimators", "get_scyama_estimators", "get_scop_estimators"]
)
def test_scaled_coordinates_free_particle_limit(method):
    s = harmonic_system(spring=0.0)
    t = s.ensemble.temp
    value = getattr(s.properties, method)(fd_delta=0.001)
    sign = 1 if method == "get_scop_estimators" else -1
    np.testing.assert_allclose(value, [3 * t, sign * 3 * t * t], atol=1e-13)


@pytest.mark.parametrize("displacement", [0.001, 0.0003, -0.0003])
def test_yamamoto_energy_and_beta_derivative(reference_system, displacement):
    s, r = reference_system()
    expected = scaled_derivatives(r.q, r.a, r.b, r.temperature)
    np.testing.assert_allclose(
        s.properties.get_yama_estimators(displacement), expected, rtol=2e-6, atol=2e-7
    )


@pytest.mark.parametrize("displacement", [1e-5, 1e-6])
def test_operator_scaled_coordinates_double_virial(reference_system, displacement):
    # Kapil et al., JCTC 15,3237 (2019), Eqs.10,14,15:
    # f_dcv = 3*F/2 - H*(q-qc)/2, using even (zero-based) beads.
    s, r = reference_system(alpha=0)
    energy, force, hessian = polynomial(r.q, r.a, r.b)
    delta = r.q - r.q.mean(axis=0)
    eop = 3 * r.temperature + np.mean(
        energy[::2] - np.sum(delta[::2] * force[::2], axis=1) / 2
    )
    double_virial_force = 1.5 * force - 0.5 * hessian * delta
    derivative = 3 * r.temperature**2 + r.temperature / r.count * np.sum(
        delta[::2] * double_virial_force[::2]
    )
    np.testing.assert_allclose(
        s.properties.get_scop_estimators(displacement),
        [eop, derivative],
        rtol=2e-6,
        atol=2e-8,
    )


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="BUG: SC scaled-coordinate derivatives omit the explicit beta dependence of the force-squared coefficient (Yamamoto Eqs.33-34)",
)
@pytest.mark.parametrize("alpha", [0.0, 0.3])
def test_suzuki_chin_scaled_coordinates(reference_system, alpha):
    s, r = reference_system(alpha=alpha)
    expected = sc_scaled_derivatives(r.q, r.masses, r.a, r.b, alpha, r.temperature)
    np.testing.assert_allclose(
        s.properties.get_scyama_estimators(0.0003), expected, rtol=3e-6, atol=2e-7
    )


@pytest.mark.parametrize("alpha", [0.5, 2.0])
@pytest.mark.parametrize(
    "atom,indices",
    [("", [0, 1]), ("H", [0]), ("1", [1])],
    ids=["all-atoms", "hydrogen", "oxygen-index"],
)
def test_free_particle_isotope_weight_is_spring_boltzmann_ratio(alpha, atom, indices):
    s = harmonic_system(spring=0.0)
    q = s.beads.q.copy().reshape(4, 2, 3)
    differences = q - np.roll(q, 1, axis=0)
    temp = s.ensemble.temp
    action = (
        0.5
        * np.array([2.0, 3.0])
        * (4 * temp) ** 2
        * np.sum(differences**2, axis=(0, 2))
    )
    delta = (alpha - 1) * action[indices]
    weight = alpha**4.5 * np.mean(np.exp(-delta / (4 * temp)))
    np.testing.assert_allclose(
        s.properties.get_isotope_zetatd(str(alpha), atom),
        [delta.mean(), np.mean(delta**2), weight],
        rtol=1e-12,
        atol=1e-15,
    )
    per_atom = (alpha - 1) * action
    per_atom[[i for i in range(2) if i not in indices]] = 0
    np.testing.assert_allclose(
        s.trajs.get_isotope_zetatd(str(alpha), atom).reshape(2, 3),
        np.column_stack([per_atom, per_atom**2, np.exp(-per_atom / (4 * temp))]),
        atol=1e-14,
    )


@pytest.mark.parametrize("alpha", [0.5, 2.0])
@pytest.mark.parametrize("atom,indices", [("", [0, 1]), ("H", [0]), ("1", [1])])
def test_scaled_isotope_weights_from_harmonic_variance(alpha, atom, indices):
    s = harmonic_system()
    q = s.beads.q.copy()
    centered = (q - q.mean(axis=0)).reshape(4, 2, 3)
    delta = 0.25 * (1 / alpha - 1) * np.sum(centered**2, axis=(0, 2)) / 4
    delta[[i for i in range(2) if i not in indices]] = 0
    expected = np.column_stack([delta, delta**2, np.exp(-delta / s.ensemble.temp)])
    np.testing.assert_allclose(
        s.trajs.get_isotope_zetasc(str(alpha), atom).reshape(2, 3), expected, atol=1e-13
    )
    np.testing.assert_array_equal(s.beads.q, q)


@pytest.mark.parametrize("mass_ratio", [0.5, 2.0])
@pytest.mark.parametrize("atom", [0, 1])
def test_fourth_order_isotope_weights(reference_system, mass_ratio, atom):
    s, r = reference_system(alpha=0)
    base = polynomial(r.q, r.a, r.b)
    omega2 = (r.count * r.temperature) ** 2

    def actions(q, masses):
        energy, force, _ = polynomial(q, r.a, r.b)
        norm = np.sum(force**2 / np.repeat(masses, 3), axis=1)
        ti = norm.sum() / (24 * omega2)
        chin = (energy[1::2].sum() - energy[::2].sum()) / 3 + norm[1::2].sum() / (
            9 * omega2
        )
        return ti, chin

    ti0, chin0 = actions(r.q, r.masses)
    exponent = -chin0 / (r.count * r.temperature)
    np.testing.assert_allclose(
        s.properties.get_chin_correction(), [exponent, exponent**2, np.exp(exponent)]
    )
    masses = r.masses.copy()
    masses[atom] *= mass_ratio
    ti1, chin1 = actions(r.q, masses)
    action = (mass_ratio - 1) * spring_action(
        r.q[:, 3 * atom : 3 * atom + 3], [r.masses[atom]], r.temperature
    )
    prefactor = mass_ratio ** (1.5 * (r.count - 1))
    expected_td = [
        action,
        action**2,
        prefactor * np.exp(-action / (r.count * r.temperature)),
        prefactor * np.exp(-(action + ti1 - ti0) / (r.count * r.temperature)),
        prefactor * np.exp(-(action + chin1 - chin0) / (r.count * r.temperature)),
    ]
    np.testing.assert_allclose(
        s.properties.get_isotope_zetatd_4th(str(mass_ratio), str(atom)),
        expected_td,
        rtol=1e-11,
    )
    shifted = r.q.copy()
    columns = slice(3 * atom, 3 * atom + 3)
    center = r.q[:, columns].mean(axis=0)
    shifted[:, columns] = center + (r.q[:, columns] - center) / np.sqrt(mass_ratio)
    ti2, chin2 = actions(shifted, masses)
    action = (polynomial(shifted, r.a, r.b)[0] - base[0]).sum()
    expected_sc = [
        action,
        action**2,
        np.exp(-action / (r.count * r.temperature)),
        np.exp(-(action + ti2 - ti0) / (r.count * r.temperature)),
        np.exp(-(action + chin2 - chin0) / (r.count * r.temperature)),
    ]
    np.testing.assert_allclose(
        s.properties.get_isotope_zetasc_4th(str(mass_ratio), str(atom)),
        expected_sc,
        rtol=1e-11,
    )


@pytest.mark.parametrize("estimator", ["get_isotope_yama", "get_isotope_thermo"])
@pytest.mark.parametrize(
    "spring",
    [
        0.5,
        pytest.param(
            -0.5,
            marks=pytest.mark.xfail(
                strict=True,
                raises=AssertionError,
                reason="Signed log accumulation gives an incorrect negative isotope kinetic estimator",
            ),
        ),
    ],
)
def test_unit_mass_ratio_recovers_mean_kinetic_estimator(estimator, spring):
    s = harmonic_system(spring=spring)
    s.beads.q_[:, 0] = [-2, -1, 1, 2]
    s.beads.q_[:, 3] = [-4, -2, 2, 4]
    q = s.beads.q.copy()
    estimate = getattr(s.properties, estimator)(alpha="1")
    # The atoms' x variances are 2.5 and 10. The mean virial term is k*12.5/4.
    expected = 1.5 * TEMPERATURE_300K + spring * 12.5 / 4
    assert estimate[0] == pytest.approx(0, abs=1e-12)
    assert np.exp(estimate[4]) == pytest.approx(2)
    assert estimate[2] == pytest.approx(expected)
    assert estimate[6] * np.exp(estimate[5] - estimate[4]) == pytest.approx(expected)
    np.testing.assert_array_equal(s.beads.q, q)


@pytest.mark.parametrize("estimator", ["get_isotope_yama", "get_isotope_thermo"])
@pytest.mark.parametrize("alpha", [0.5, 2.0])
def test_isotope_reweighting_at_changed_mass(estimator, alpha):
    # k = 1/2, T = 0.2, P = 4, m = [2, 3]; only x of each atom is spread.
    s = harmonic_system(temperature=0.2)
    positions = np.zeros((4, 6))
    positions[:, 0] = [-0.2, -0.1, 0.1, 0.2]
    positions[:, 3] = [-0.4, -0.2, 0.2, 0.4]
    s.beads.q_[:] = positions
    variance = np.array([0.025, 0.1])
    if estimator == "get_isotope_yama":
        # Scaling x by 1/sqrt(alpha) scales its variance by 1/alpha.
        kinetic = 0.3 + 0.25 * variance / alpha
        log_ratio = 0.25 * (1 / alpha - 1) * variance / 0.2
    else:
        # The spring action m (P T)^2/2 sum((x_b - x_(b-1))^2) scales with alpha;
        # the closed-path sums of squared differences are [0.22, 0.88].
        kinetic = 0.3 + 0.25 * variance
        spring_energy = 0.5 * np.array([2.0, 3.0]) * 0.8**2 * np.array([0.22, 0.88])
        log_ratio = (alpha - 1) * spring_energy / 0.8
    weights = np.exp(-log_ratio)

    estimate = getattr(s.properties, estimator)
    for atom in range(2):
        actual = estimate(alpha=str(alpha), atom=str(atom))
        assert actual[0] == pytest.approx(log_ratio[atom], rel=1e-12, abs=1e-12)
        assert actual[2] == pytest.approx(kinetic[atom], rel=1e-12, abs=1e-12)
        assert np.exp(actual[4]) == pytest.approx(weights[atom], rel=1e-12, abs=1e-12)
    expected = [
        np.mean(log_ratio),
        np.mean(log_ratio**2),
        np.mean(kinetic),
        np.mean(kinetic**2),
        np.log(np.sum(weights)),
        np.log(np.dot(weights, kinetic)),
        1.0,
    ]
    np.testing.assert_allclose(
        estimate(alpha=str(alpha)), expected, rtol=1e-12, atol=1e-12
    )
    np.testing.assert_array_equal(s.beads.q, positions)
    assert s.forces.pot == pytest.approx(0.125, abs=1e-12)


@pytest.mark.parametrize("separation", [0.15, 0.4])
def test_exchange_estimators_sum_both_permutations(separation):
    root = ET.fromstring(harmonic_input(nbeads=4))
    root.find("system/ensemble/temperature").attrib.clear()
    root.find("system/ensemble/temperature").text = "0.4"
    root.find("system/beads/m").text = "[2,2]"
    normal = ET.SubElement(root.find("system"), "normal_modes")
    ET.SubElement(normal, "bosons", id="index").text = "[0,1]"
    s = Simulation.load_from_xml(
        ET.tostring(root, encoding="unicode"), read_only=True
    ).syslist[0]
    q = np.arange(24).reshape(4, 2, 3) * 0.013
    q[:, 1] += separation
    s.beads.q_[:] = q.reshape(4, 6)
    ref = permutation_reference(q, 2, 0.4)
    for key, expected in [
        ("exchange_distinct_prob", ref["distinct"]),
        ("exchange_all_prob", ref["longest"]),
        ("fermionic_sign", ref["sign"]),
    ]:
        assert s.properties[key][0] == pytest.approx(expected, rel=1e-12)
    assert s.properties.get_kintd() == pytest.approx(ref["kinetic"])
    with pytest.raises(IndexError, match="not applicable to bosons"):
        s.properties.get_kincv()
    with pytest.raises(IndexError, match="subset"):
        s.properties.get_kintd("0")


@pytest.mark.parametrize(
    "call,error",
    [
        pytest.param(lambda s: s.properties.get_kinmd(atom="2"), IndexError, id="atom"),
        pytest.param(lambda s: s.properties.get_kinmd(bead="4"), IndexError, id="bead"),
        pytest.param(lambda s: s.properties.get_kinmd(nm="4"), IndexError, id="mode"),
        pytest.param(
            lambda s: s.properties.get_kinmd(bead="0", nm="0"),
            ValueError,
            id="bead-and-mode",
        ),
        pytest.param(
            lambda s: s.properties.get_kinmd(bead="first"), ValueError, id="bead-text"
        ),
        pytest.param(
            lambda s: s.properties.get_kinmd(nm="invalid"), ValueError, id="mode-text"
        ),
        pytest.param(
            lambda s: s.properties.get_atom_vec(s.beads.q, "", "0"),
            IndexError,
            id="atom_vec-no-atom",
        ),
        pytest.param(
            lambda s: s.properties.get_atom_vec(s.beads.q, "0", "4"),
            IndexError,
            id="atom_vec-bead",
        ),
        pytest.param(lambda s: s.properties.get_kij("2", "0"), IndexError, id="kij"),
        pytest.param(lambda s: s.properties.get_kij("0", "2"), IndexError, id="kij"),
        *[
            pytest.param(
                lambda s, method=method: getattr(s.properties, method)(atom="2"),
                IndexError,
                id=method,
            )
            for method in [
                "get_atom_ids",
                "get_kincv",
                "get_kintd",
                "get_sckinop",
                "get_sckintd",
                "get_ktens",
                "get_rg",
                "get_linlin",
                "get_isotope_yama",
                "get_isotope_thermo",
                "get_isotope_zetatd",
                "get_isotope_zetasc",
                "get_isotope_zetatd_4th",
                "get_isotope_zetasc_4th",
            ]
        ],
        *[
            pytest.param(
                lambda s, method=method: getattr(s.trajs, method)(atom="2"),
                IndexError,
                id="trajectory-" + method,
            )
            for method in ["get_isotope_zetatd", "get_isotope_zetasc"]
        ],
        *[
            pytest.param(
                lambda s, method=method: getattr(s.properties, method)(),
                Exception,
                id=method,
            )
            for method in [
                "get_exchange_distinct_prob",
                "get_exchange_longest_prob",
                "get_fermionic_sign",
            ]
        ],
    ],
)
def test_invalid_selections_are_rejected(call, error):
    with pytest.raises(error):
        call(harmonic_system())

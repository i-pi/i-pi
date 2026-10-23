"""Independent numerical references for the property/trajectory checklist.

Tests use real bound systems and a prescribed convex polynomial potential.
The test driver differentiates its explicit expression; expectations use
numpy.polynomial and the independent references in mathematical_references.
"""

import json
from types import SimpleNamespace
import xml.etree.ElementTree as ET

import numpy as np
import pytest

from ipi.engine.simulation import Simulation
from ipi.engine.motion.driven_dynamics import DrivenDynamics, ElectricField, BEC
from ipi_tests.unit_tests.engine.test_simulation_invariants import harmonic_input
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

SC_VIRIAL_BUG = pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="BUG: SC virial omits the force-product term required by the complete energy's cell-strain derivative (JCTC 2019, SI S1.31)",
)


def verifies(*entries):
    """Attach checklist identities for the external report builder."""

    def decorate(function):
        function.independent_checks = entries
        return function

    return decorate


class PrescribedPolynomial:
    """PES input U=sum(a*x^2+b*x^4); no reference helper computes its output."""

    def __init__(self, quadratic, quartic, electric=False):
        self.a = np.array(quadratic)
        self.b = np.array(quartic)
        self.electric = electric

    def __call__(self, cell, positions):
        x = positions.reshape(-1)
        energy = np.sum(self.a * x * x + self.b * x * x * x * x)
        force = -2 * self.a * x - 4 * self.b * x * x * x
        virial = np.zeros((3, 3))
        for atom in range(len(positions)):
            for i in range(3):
                for j in range(3):
                    virial[i, j] += x[3 * atom + i] * force[3 * atom + j]
        extras = {"energy": float(energy)}
        if self.electric:
            extras["dipole"] = (positions[0] - positions[1]).tolist()
        extra = json.dumps(extras)
        return energy, force.reshape(positions.shape), virial, extra


@pytest.fixture
def reference_system():
    def create(count=6, alpha=0.3, style="rpmd", frequencies=(), quartic=True):
        root = ET.fromstring(harmonic_input(nbeads=count))
        temp = root.find("system/ensemble/temperature")
        temp.attrib.clear()
        temp.text = "0.4"
        normal = ET.SubElement(root.find("system"), "normal_modes", transform="matrix")
        ET.SubElement(normal, "frequencies", style=style).text = str(list(frequencies))
        s = Simulation.load_from_xml(
            ET.tostring(root, encoding="unicode"), read_only=True
        ).syslist[0]
        a = np.array([0.35, 0.55, 0.8, 0.45, 0.65, 0.95])
        b = np.array([0.03, 0.05, 0.02, 0.04, 0.01, 0.06]) if quartic else np.zeros(6)
        driver = PrescribedPolynomial(a, b)
        s.forces.ff["harmonic"].driver = driver
        q = (
            0.35 * np.sin(np.arange(count * 6).reshape(count, 6) * 0.71)
            + np.arange(6) * 0.07
        )
        p = (
            0.23 * np.cos(np.arange(count * 6).reshape(count, 6) * 0.47)
            + np.arange(6) * 0.03
        )
        masses = np.array([2.0, 3.0])
        h = np.array([[6.0, 1.0, 0.5], [0.0, 5.0, 0.7], [0.0, 0.0, 4.0]])
        s.beads.q[:] = q
        s.beads.p[:] = p
        s.beads.m[:] = masses
        s.cell.h[:] = h
        s.forces.alpha = alpha
        return s, SimpleNamespace(
            q=q,
            p=p,
            masses=masses,
            a=a,
            b=b,
            h=h,
            temperature=0.4,
            count=count,
            alpha=alpha,
        )

    return create


def output(s, key, trajectory=False):
    return (s.trajs if trajectory else s.properties)[key][0]


@verifies(
    "properties:density",
    "properties:volume",
    "properties:cell_h",
    "properties:cell_abcABC",
)
@pytest.mark.parametrize("count", [2, 6])
def test_cell_and_density_outputs_match_lattice_geometry(reference_system, count):
    s, r = reference_system(count)
    assert output(s, "density") == pytest.approx(5 / 120)
    assert output(s, "volume") == pytest.approx(120)
    np.testing.assert_allclose(output(s, "cell_h"), [6, 5, 4, 1, 0.5, 0.7])
    lengths = np.linalg.norm(r.h, axis=0)
    angles = [
        np.degrees(np.arccos(np.dot(r.h[:, i], r.h[:, j]) / (lengths[i] * lengths[j])))
        for i, j in [(1, 2), (0, 2), (0, 1)]
    ]
    np.testing.assert_allclose(output(s, "cell_abcABC"), np.r_[lengths, angles])


@verifies(
    "properties:atom_x_path",
    "properties:atom_f_path",
    "properties:atom_v",
    "properties:atom_p",
    "properties:atom_f",
    "trajectories:velocities",
    "trajectories:v_centroid",
    "trajectories:x_centroid",
    "trajectories:p_centroid",
    "trajectories:f_centroid",
    "trajectories:x_centroid_even",
    "trajectories:x_centroid_odd",
    "trajectories:v_centroid_even",
    "trajectories:v_centroid_odd",
)
@pytest.mark.parametrize("atom", [0, 1])
def test_atomic_and_centroid_outputs_preserve_values_and_bead_order(
    reference_system, atom
):
    s, r = reference_system()
    _, force, _ = polynomial(r.q, r.a, r.b)
    columns = slice(3 * atom, 3 * atom + 3)
    velocity = r.p / np.repeat(r.masses, 3)
    np.testing.assert_allclose(
        output(s, f"atom_x_path({atom})"), r.q[:, columns].ravel()
    )
    np.testing.assert_allclose(
        output(s, f"atom_f_path({atom})"), force[:, columns].ravel()
    )
    for key, values in [("atom_v", velocity), ("atom_p", r.p), ("atom_f", force)]:
        for bead in [-1, 2]:
            expected = (
                values[:, columns].mean(axis=0) if bead < 0 else values[bead, columns]
            )
            np.testing.assert_allclose(output(s, f"{key}({atom};{bead})"), expected)
    for key, values in [
        ("velocities", velocity),
        ("v_centroid", velocity.mean(axis=0)),
        ("x_centroid", r.q.mean(axis=0)),
        ("p_centroid", r.p.mean(axis=0)),
        ("f_centroid", force.mean(axis=0)),
    ]:
        np.testing.assert_allclose(output(s, key, True), values)
    for parity, selected in [
        ("even", slice(None, None, 2)),
        ("odd", slice(1, None, 2)),
    ]:
        np.testing.assert_allclose(
            output(s, "x_centroid_" + parity, True), r.q[selected].mean(axis=0)
        )
        np.testing.assert_allclose(
            output(s, "v_centroid_" + parity, True), velocity[selected].mean(axis=0)
        )


@verifies(
    "properties:spring",
    "properties:kinetic_td",
    "properties_helpers:Properties._kinetic_td_distinguishables",
    "trajectories:forces_spring",
    "properties:conserved",
    "properties:ensemble_bias",
)
@pytest.mark.parametrize("count", [2, 4, 6])
def test_spring_primitive_conserved_and_bias_outputs_match_explicit_actions(
    reference_system, count
):
    from ipi.utils.depend import depend_value

    s, r = reference_system(count)
    action = spring_action(r.q, r.masses, r.temperature)
    force = (
        -((count * r.temperature) ** 2)
        * np.repeat(r.masses, 3)
        * (2 * r.q - np.roll(r.q, 1, axis=0) - np.roll(r.q, -1, axis=0))
    )
    assert output(s, "spring") == pytest.approx(action / count)
    assert output(s, "kinetic_td") == pytest.approx(
        3 * count * r.temperature - action / count
    )
    np.testing.assert_allclose(output(s, "forces_spring", True), force, atol=1e-13)
    bias = depend_value(name="pot", value=1.7)
    s.ensemble.bias._pot = bias
    s.ensemble._econs.add_dependency(bias)
    energy = polynomial(r.q, r.a, r.b)[0].sum()
    kinetic = np.sum(r.p**2 / (2 * np.repeat(r.masses, 3)))
    assert output(s, "ensemble_bias") == pytest.approx(1.7 / count)
    assert output(s, "conserved") == pytest.approx(
        (action + energy + kinetic + 1.7) / count
    )


@verifies("properties:kinetic_md", "properties_helpers:Properties.get_kstress_md")
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
def test_mass_scaled_kinetic_energy_and_stress_use_prescribed_mode_masses(
    reference_system, style, frequencies, atom
):
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


@verifies("properties:temperature")
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
def test_temperature_restores_only_removed_gaussian_degrees_of_freedom(
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
    s.beads.p[:] = fourier_matrix(4).T @ mode_p
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


@verifies("properties:Eenvelope", "trajectories:Eforces")
@pytest.mark.parametrize("time", [0.1, 0.7, 1.2])
def test_driven_output_wrappers_match_gaussian_envelope_and_charge_field_product(
    reference_system, time
):
    s, r = reference_system()
    charge = np.arange(r.count * 6 * 3).reshape(r.count, 6, 3) / 17
    motion = DrivenDynamics(
        timestep=0.01,
        mode="eda-nve",
        efield=ElectricField(
            amp=np.array([0.2, -0.3, 0.4]), freq=0.7, phase=0.2, peak=0.6, sigma=0.3
        ),
        bec=BEC(bec=charge),
    )
    motion.bind(
        s.ensemble, s.beads, s.nm, s.cell, s.forces, s.simul.prng, s.simul.output_maker
    )
    s.motion = s.properties.motion = motion
    s.ensemble.time = time
    envelope = np.exp(-0.5 * ((time - 0.6) / 0.3) ** 2)
    field = envelope * np.cos(0.7 * time + 0.2) * np.array([0.2, -0.3, 0.4])
    assert output(s, "Eenvelope") == pytest.approx(envelope)
    np.testing.assert_allclose(output(s, "Eforces", True), charge @ field, atol=1e-13)


@verifies("properties:kinetic_ij")
@pytest.mark.parametrize("first,second", [(0, 1), (1, 0), (0, 0), (1, 1)])
def test_cross_atom_kinetic_tensor_has_independent_mass_symmetrization(
    reference_system, first, second
):
    s, r = reference_system()
    force = polynomial(r.q, r.a, r.b)[1].reshape(r.count, 2, 3)
    delta = (r.q - r.q.mean(axis=0)).reshape(r.count, 2, 3)
    tensor = -(
        r.masses[first] * delta[:, first].T @ force[:, second]
        + r.masses[second] * force[:, first].T @ delta[:, second]
    ) / (4 * r.count * np.sqrt(r.masses[first] * r.masses[second]))
    if first == second:
        tensor += np.eye(3) * r.temperature / 2
    expected = tensor[[0, 1, 2, 0, 0, 1], [0, 1, 2, 1, 2, 2]]
    np.testing.assert_allclose(
        s.properties.get_kij(str(first), str(second)), expected, atol=1e-13
    )


@verifies(
    "properties:potential_opsc",
    "properties:potential_tdsc",
    "properties:kinetic_prsc",
    "properties:kinetic_tdsc",
    "properties:kinetic_opsc",
    "properties_helpers:Properties.kstress_sctd",
    "trajectories:forces_sc",
)
@pytest.mark.parametrize("alpha", [0.0, 0.3, 1.0])
def test_nonzero_suzuki_chin_estimators_follow_published_action(
    reference_system, alpha
):
    s, r = reference_system(alpha=alpha)
    ref = sc_reference(r.q, r.masses, r.a, r.b, alpha, r.temperature)
    w, c = sc_weights(r.count, alpha, r.temperature)
    correction_force = ref["force1"] + ref["force2"]
    delta = r.q - r.q.mean(axis=0)
    cv = (
        3 * r.temperature
        - np.sum(delta * ref["force"][::1] * np.tile([1, 0], r.count // 2)[:, None])
        / r.count
    )
    assert output(s, "potential_opsc") == pytest.approx(ref["energy"][::2].mean())
    assert output(s, "potential_tdsc") == pytest.approx(
        np.mean(w * ref["energy"] + 2 * c * ref["square"])
    )
    assert output(s, "kinetic_prsc") == pytest.approx(
        3 * r.count * r.temperature
        - spring_action(r.q, r.masses, r.temperature) / r.count
        + np.mean(c * ref["square"])
    )
    assert output(s, "kinetic_opsc") == pytest.approx(cv)
    expected_td = (
        3 * r.temperature
        - np.sum(delta * (ref["force"] + correction_force)) / (2 * r.count)
        + np.mean(c * ref["square"])
    )
    assert output(s, "kinetic_tdsc") == pytest.approx(expected_td, rel=1e-8, abs=1e-10)
    np.testing.assert_allclose(
        output(s, "forces_sc", True),
        ref["force"] + correction_force,
        rtol=1e-8,
        atol=1e-10,
    )


@verifies(
    "properties:stress_md",
    "properties:pressure_md",
    "properties:kstress_md",
    "properties:virial_md",
    "properties:stress_cv",
    "properties:pressure_cv",
    "properties:kstress_cv",
    "properties:virial_cv",
    "properties:kstress_tdsc",
)
@pytest.mark.parametrize("alpha", [0.0, 0.3])
@pytest.mark.parametrize(
    "key",
    [
        "stress_md",
        "pressure_md",
        "kstress_md",
        "virial_md",
        "stress_cv",
        "pressure_cv",
        "kstress_cv",
        "virial_cv",
        "kstress_tdsc",
    ],
)
def test_pressure_and_stress_outputs_have_independent_volume_and_bead_factors(
    reference_system, alpha, key
):
    s, r = reference_system(alpha=alpha)
    ref = sc_reference(r.q, r.masses, r.a, r.b, alpha, r.temperature)
    delta = r.q - r.q.mean(axis=0)
    pc = r.p.mean(axis=0).reshape(2, 3)
    # An unbiased CV stress estimator may replace off-diagonal centroid
    # momentum products by their Maxwell expectation, zero. Retain the
    # instantaneous squared centroid momenta only on the diagonal.
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
    virial = ref["virial"].sum(axis=0)

    def vector(tensor):
        return tensor[[0, 1, 2, 0, 0, 1], [0, 1, 2, 1, 2, 2]]

    expected = dict(
        stress_md=vector((virial + md) / 120),
        pressure_md=np.trace(virial + md) / 360,
        kstress_md=vector(md / 120),
        virial_md=vector(virial / 120),
        stress_cv=vector((virial + cv) / (120 * r.count)),
        pressure_cv=np.trace(virial + cv) / (360 * r.count),
        kstress_cv=vector(cv / (120 * r.count)),
        virial_cv=vector(virial / (120 * r.count)),
        kstress_tdsc=np.trace(sc) / (360 * r.count),
    )
    np.testing.assert_allclose(
        output(s, key), expected[key], rtol=1e-8, atol=1e-10, err_msg=key
    )
    if key == "kstress_tdsc":
        np.testing.assert_allclose(
            s.properties.kstress_sctd(), sc, rtol=1e-8, atol=1e-10
        )


@verifies("properties:displacedpath", "properties_helpers:Properties.opening")
@pytest.mark.parametrize("atom", [0, 1])
def test_displaced_path_nonzero_potential_matches_explicit_opened_path(
    reference_system, atom
):
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
    assert output(s, f"displacedpath(.17;-.23;.11;{atom})") == pytest.approx(expected)


@verifies("properties:scaledcoords")
@pytest.mark.parametrize("displacement", [0.001, 0.0003, -0.0003])
def test_yamamoto_nonharmonic_energy_and_beta_derivative(
    reference_system, displacement
):
    s, r = reference_system()
    expected = scaled_derivatives(r.q, r.a, r.b, r.temperature)
    np.testing.assert_allclose(
        s.properties.get_yama_estimators(displacement), expected, rtol=2e-6, atol=2e-7
    )


@verifies("properties:sc_op_scaledcoords")
@pytest.mark.parametrize("displacement", [1e-5, 1e-6])
def test_operator_scaled_coordinates_match_published_double_virial_expression(
    reference_system, displacement
):
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


@verifies("properties:sc_scaledcoords")
@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="BUG: SC scaled-coordinate derivatives omit the explicit beta dependence of the force-squared coefficient (Yamamoto Eqs.33-34)",
)
@pytest.mark.parametrize("alpha", [0.0, 0.3])
def test_suzuki_chin_scaled_derivatives_include_explicit_beta_dependence(
    reference_system, alpha
):
    s, r = reference_system(alpha=alpha)
    expected = sc_scaled_derivatives(r.q, r.masses, r.a, r.b, alpha, r.temperature)
    np.testing.assert_allclose(
        s.properties.get_scyama_estimators(0.0003), expected, rtol=3e-6, atol=2e-7
    )


@verifies(
    "properties:chin_weight",
    "properties:isotope_zetatd_4th",
    "properties:isotope_zetasc_4th",
)
@pytest.mark.parametrize("mass_ratio", [0.5, 2.0])
@pytest.mark.parametrize("atom", [0, 1])
def test_fourth_order_isotope_weights_include_nonzero_force_and_potential_changes(
    reference_system, mass_ratio, atom
):
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


@verifies(
    "properties:exchange_distinct_prob",
    "properties:exchange_all_prob",
    "properties:fermionic_sign",
)
@pytest.mark.parametrize("separation", [0.15, 0.4])
def test_exchange_outputs_match_both_explicit_permutations(
    reference_system, separation
):
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
    s.beads.q[:] = q.reshape(4, 6)
    ref = permutation_reference(q, 2, 0.4)
    for key, expected in [
        ("exchange_distinct_prob", ref["distinct"]),
        ("exchange_all_prob", ref["longest"]),
        ("fermionic_sign", ref["sign"]),
    ]:
        assert output(s, key) == pytest.approx(expected, rel=1e-12)
    assert s.properties.get_kintd() == pytest.approx(ref["kinetic"])

"""Force arithmetic checked with explicit Fourier projection and chain rules."""

import numpy as np
import pytest

from ipi.engine.forcefields import FFDirect
from ipi.engine.forces import Forces, ForceComponent, ScaledForceComponent, MTSForces
from ipi_tests.unit_tests.engine.test_property_mathematical_references import (
    PrescribedPolynomial,
    reference_system as _reference_system,
    output,
    verifies,
)
from ipi_tests.unit_tests.engine.mathematical_references import (
    polynomial,
    projection,
    sc_reference,
    virials,
)

reference_system = _reference_system


def combined_case(s, r, contracted=True):
    components = [
        ForceComponent(
            "first",
            weight=0.7,
            mts_weights=[0.2, 0.4],
            epsilon=-1e-5,
            interpolate_extras=["energy"],
        ),
        ForceComponent(
            "second",
            nbeads=3 if contracted else 6,
            weight=-0.3,
            mts_weights=[0.8, 0.6],
            epsilon=-1e-5,
            interpolate_extras=["energy"],
        ),
    ]
    providers = {}
    expected = []
    for index, component in enumerate(components):
        provider = FFDirect(pes="harmonic", pars={"k1": 1.0})
        a, b = r.a * (1 + index * 0.4), r.b * (1 + index * 0.7)
        provider.driver = PrescribedPolynomial(a, b)
        providers[component.ffield] = provider
        count = component.nbeads or 6
        contract, expand = projection(6, count), projection(count, 6)
        q = contract @ r.q
        energy, force, hessian = polynomial(q, a, b)
        expected.append(
            dict(
                component=component,
                q=q,
                energy=energy,
                force=force,
                curvature=hessian,
                virial=virials(q, force),
                contract=contract,
                expand=expand,
                a=a,
                b=b,
            )
        )
    result = Forces()
    result.bind(s.beads, s.cell, components, providers, [], None)
    result.alpha = r.alpha
    result.omegan2 = (6 * r.temperature) ** 2
    return result, expected


def expand(matrix, values):
    return np.tensordot(matrix, values, axes=1)


@verifies(
    "forces:ForceComponent.bind",
    "forces:ScaledForceComponent.__init__",
    "forces:MTSForces.get_forces_mts",
    "forces:MTSForces.get_vir_mts",
    "forces:MTSForces.get_virs_mts",
    "forces:Forces.bind",
    "forces:Forces.bind.make_rpc",
    "forces:Forces.virs_component",
    "forces:Forces.f_combine",
    "forces:Forces.pot_combine",
    "forces:Forces.vir_combine",
    "forces:Forces.extra_combine",
)
@pytest.mark.parametrize("contracted", [False, True])
def test_noncentroid_components_and_mts_levels_match_explicit_fourier_projection(
    reference_system, contracted
):
    s, r = reference_system()
    forces, refs = combined_case(s, r, contracted)
    for actual, ref in zip(forces.mforces, refs):
        np.testing.assert_allclose(actual.beads.q, ref["q"], atol=1e-13)
        assert actual.pot == pytest.approx(ref["energy"].sum())
        scaled = ScaledForceComponent(actual, scaling=-0.4)
        assert scaled.pot == pytest.approx(-0.4 * ref["energy"].sum())
    totals = {}
    for quantity in ["force", "energy", "virial"]:
        totals[quantity] = sum(
            x["component"].weight
            * sum(x["component"].mts_weights)
            * expand(x["expand"], x[quantity])
            for x in refs
        )
    np.testing.assert_allclose(forces.f, totals["force"], atol=1e-12)
    np.testing.assert_allclose(forces.pots, totals["energy"], atol=1e-12)
    np.testing.assert_allclose(forces.virs, totals["virial"], atol=1e-12)
    np.testing.assert_allclose(forces.extras["energy"], totals["energy"], atol=1e-12)
    for level in range(2):
        expected_force = sum(
            x["component"].weight
            * x["component"].mts_weights[level]
            * expand(x["expand"], x["force"])
            for x in refs
        )
        expected_virial = sum(
            x["component"].weight
            * x["component"].mts_weights[level]
            * expand(x["expand"], x["virial"])
            for x in refs
        )
        np.testing.assert_allclose(
            forces.mts_forces[level].f, expected_force, atol=1e-12
        )
        np.testing.assert_allclose(
            forces.mts_forces[level].virs, expected_virial, atol=1e-12
        )
        np.testing.assert_allclose(
            forces.mts_forces[level].vir, expected_virial.sum(axis=0), atol=1e-12
        )
    for index, ref in enumerate(refs):
        for weighted in [False, True]:
            for interpolate in [False, True]:
                expected = ref["virial"].copy()
                if weighted:
                    expected *= ref["component"].weight
                if interpolate:
                    expected = expand(ref["expand"], expected)
                np.testing.assert_allclose(
                    forces.virs_component(index, weighted, interpolate),
                    expected,
                    atol=1e-12,
                )


@verifies(
    "properties:pot_component",
    "properties:pot_component_raw",
    "properties:vir_component",
    "properties:vir_component_raw",
    "trajectories:forces_component",
    "trajectories:forces_component_raw",
)
@pytest.mark.parametrize(
    "key",
    [
        "pot_component",
        "pot_component_raw",
        "vir_component",
        "vir_component_raw",
        "forces_component",
        "forces_component_raw",
    ],
)
@pytest.mark.parametrize("bead", [-1, 1])
def test_component_output_wrappers_preserve_raw_and_weighted_bead_values(
    reference_system, key, bead, request
):
    if key == "pot_component_raw" and bead < 0:
        request.node.add_marker(
            pytest.mark.xfail(
                strict=True,
                raises=AssertionError,
                reason="BUG: raw component potential mean divides by the full bead count instead of the component bead count",
            )
        )
    s, r = reference_system()
    forces, refs = combined_case(s, r)
    s.properties.forces = s.forces = forces
    raw = key.endswith("_raw")
    quantity = (
        "energy"
        if key.startswith("pot")
        else "virial" if key.startswith("vir") else "force"
    )
    ref = refs[1]
    expected = ref[quantity].copy()
    if not raw:
        expected = -0.3 * expand(ref["expand"], expected)
    expected = expected.mean(axis=0) if bead < 0 else expected[bead]
    if quantity == "virial":
        expected = expected[[0, 1, 2, 0, 0, 1], [0, 1, 2, 1, 2, 2]] / 120
    np.testing.assert_allclose(
        output(s, f"{key}(1;{bead})", quantity == "force"), expected, atol=1e-12
    )


@verifies(
    "forces:Forces.forcesvirs_4th_order",
    "forces:Forces.fvir_4th_order_combine",
    "forces:Forces.get_potssc",
    "forces:Forces.get_coeffsc_part_1",
    "forces:Forces.get_coeffsc_part_2",
    "forces:Forces.get_fsc_part_1",
    "forces:Forces.get_fsc_part_2",
    "forces:Forces.get_fsc",
    "forces:Forces.get_virssc_part_1",
    "forces:Forces.get_virssc_part_2",
    "forces:Forces.get_virssc",
)
@pytest.mark.parametrize("alpha", [0.0, 0.3, 1.0])
def test_individual_sc_coefficients_and_corrections_follow_polynomial_derivatives(
    reference_system, alpha
):
    s, r = reference_system(alpha=alpha)
    force = s.forces
    force.mforces[0].epsilon = -1e-5
    ref = sc_reference(r.q, r.masses, r.a, r.b, alpha, 0.4)
    derivative = ref["derivative"].copy()
    derivative_virial = ref["derivative_virial"].copy()
    if alpha == 0:
        derivative[::2] = 0
        derivative_virial[::2] = 0
    checks = dict(
        coeffsc_part_1=ref["linear"][:, None],
        coeffsc_part_2=ref["quadratic"][:, None],
        f_4th_order=derivative,
        vir_4th_order=derivative_virial,
        potssc=ref["correction"],
        fsc_part_1=ref["force1"],
        fsc_part_2=ref["force2"],
        fsc=ref["force1"] + ref["force2"],
        virssc_part_1=ref["virial1"],
        virssc_part_2=ref["virial2"],
        virssc=ref["virial1"] + ref["virial2"],
    )
    for name, expected in checks.items():
        np.testing.assert_allclose(
            getattr(force, name), expected, rtol=1e-7, atol=2e-8, err_msg=name
        )


@verifies(
    "forces:Forces.forcesvirs_4th_order",
    "forces:Forces.fvir_4th_order_combine",
    "forces:Forces.get_fsc",
    "forces:Forces.get_virssc",
)
@pytest.mark.parametrize("contracted", [False, True])
@pytest.mark.parametrize("alpha", [0.0, 0.3])
def test_multicomponent_sc_force_and_virial_derivatives_obey_chain_rule(
    reference_system, contracted, alpha, request
):
    if contracted and alpha == 0:
        request.node.add_marker(
            pytest.mark.xfail(
                strict=True,
                raises=AssertionError,
                reason="BUG: alpha=0 SC contraction mixes positive and negative displacements before evaluating the contracted force",
            )
        )
    s, r = reference_system(alpha=alpha)
    forces, refs = combined_case(s, r, contracted)
    total_force = sum(
        x["component"].weight
        * sum(x["component"].mts_weights)
        * expand(x["expand"], x["force"])
        for x in refs
    )
    direction = total_force / np.repeat(r.masses, 3)
    derivatives, virial_derivatives = [], []
    for index, ref in enumerate(refs):
        dq = ref["contract"] @ direction
        local = 2 * ref["curvature"] * dq
        v_local = []
        for q, f, h, d in zip(ref["q"], ref["force"], ref["curvature"], dq):
            v_local.append(
                -2
                * np.triu(
                    d.reshape(2, 3).T @ f.reshape(2, 3)
                    + q.reshape(2, 3).T @ (-h * d).reshape(2, 3)
                )
            )
        df = expand(ref["expand"], local)
        dv = expand(ref["expand"], np.asarray(v_local))
        if alpha == 0:
            df[::2] = 0
            dv[::2] = 0
        actual_f, actual_v = forces.forcesvirs_4th_order(index)
        np.testing.assert_allclose(actual_f, df, rtol=1e-7, atol=2e-8)
        np.testing.assert_allclose(actual_v, dv, rtol=1e-7, atol=2e-8)
        weight = ref["component"].weight * sum(ref["component"].mts_weights)
        derivatives.append(weight * df)
        virial_derivatives.append(weight * dv)
    np.testing.assert_allclose(
        forces.f_4th_order, sum(derivatives), rtol=1e-7, atol=2e-8
    )
    np.testing.assert_allclose(
        forces.vir_4th_order, sum(virial_derivatives), rtol=1e-7, atol=2e-8
    )


@verifies(
    "forces:MTSForces.get_forces_mts",
    "forces:MTSForces.get_virs_mts",
    "forces:MTSForces.get_vir_mts",
)
@pytest.mark.xfail(
    strict=True,
    raises=(IndexError, AssertionError),
    reason="BUG: empty MTS weights should select the outer level but force indexing fails and virials are omitted",
)
@pytest.mark.parametrize("quantity", ["f", "virs", "vir"])
def test_empty_mts_weights_place_unpartitioned_component_on_outer_level(
    reference_system, quantity
):
    s, r = reference_system()
    s.forces.mforces[0].mts_weights = np.array([])
    outer = MTSForces(s.forces, 0)
    u, f, h = polynomial(r.q, r.a, r.b)
    v = virials(r.q, f)
    expected = dict(f=f, virs=v, vir=v.sum(axis=0))[quantity]
    np.testing.assert_allclose(getattr(outer, quantity), expected, atol=1e-12)

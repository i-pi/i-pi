"""Force composition: weights, contraction, MTS levels, SC corrections and requests.

Reference values come from explicit Fourier projections and from derivatives
of the polynomial and harmonic models.
"""

import time
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from threading import Barrier
from types import SimpleNamespace

import numpy as np
import pytest

import ipi.engine.forces as force_module
from ipi.engine.beads import Beads
from ipi.engine.cell import Cell
from ipi.engine.forcefields import FFDirect, FFEval
from ipi.engine.forces import Forces, ForceComponent, ScaledForceComponent
from ipi_tests.unit_tests.engine.inputs import PrescribedPolynomial, harmonic_input
from ipi_tests.unit_tests.engine.mathematical_references import (
    polynomial,
    projection,
    sc_reference,
    virials,
)


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


@pytest.mark.parametrize("contracted", [False, True])
def test_contracted_components_and_mts_levels(reference_system, contracted):
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


@pytest.mark.parametrize(
    "key,bead",
    [
        ("pot_component", -1),
        ("pot_component", 1),
        ("pot_component_raw", 1),
        ("vir_component", -1),
        ("vir_component", 1),
        ("vir_component_raw", -1),
        ("vir_component_raw", 1),
        ("forces_component", -1),
        ("forces_component", 1),
        ("forces_component_raw", -1),
        ("forces_component_raw", 1),
    ],
)
def test_component_outputs_raw_and_weighted(reference_system, key, bead):
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
        (s.trajs if quantity == "force" else s.properties)[f"{key}(1;{bead})"][0],
        expected,
        atol=1e-12,
    )


@pytest.mark.parametrize("alpha", [0.0, 0.3, 1.0])
def test_sc_coefficients_and_corrections(reference_system, alpha):
    s, r = reference_system(alpha=alpha)
    force = s.forces
    force.mforces[0].epsilon = -1e-5
    ref = sc_reference(r.q, r.masses, r.a, r.b, alpha, 0.4)
    derivative = ref["derivative"].copy()
    if alpha == 0:
        derivative[::2] = 0
    checks = dict(
        coeffsc_part_1=ref["linear"][:, None],
        coeffsc_part_2=ref["quadratic"][:, None],
        f_4th_order=derivative,
        potssc=ref["correction"],
        fsc_part_1=ref["force1"],
        fsc_part_2=ref["force2"],
        fsc=ref["force1"] + ref["force2"],
    )
    for name, expected in checks.items():
        np.testing.assert_allclose(
            getattr(force, name), expected, rtol=1e-7, atol=2e-8, err_msg=name
        )


@pytest.mark.parametrize("before,after", [(0.0, 0.3), (0.3, 0.0), (0.3, 0.5)])
def test_sc_alpha_cannot_switch_between_zero_and_nonzero(
    reference_system, before, after
):
    """The finite differences of the fourth-order correction displace only the
    odd beads for alpha = 0, and all of them otherwise: once it is computed,
    alpha can change only between non-zero values"""
    s, r = reference_system(alpha=before)
    s.forces.mforces[0].epsilon = -1e-5
    s.forces.f_4th_order
    s.forces.alpha = after
    if (before == 0) != (after == 0):
        with pytest.raises(ValueError, match="alpha"):
            s.forces.f_4th_order
    else:
        derivative = sc_reference(r.q, r.masses, r.a, r.b, after, 0.4)["derivative"]
        np.testing.assert_allclose(
            s.forces.f_4th_order, derivative, rtol=1e-7, atol=2e-8
        )


@pytest.mark.parametrize("contracted", [False, True])
@pytest.mark.parametrize("alpha", [0.0, 0.3])
def test_multicomponent_sc_force_derivatives_obey_chain_rule(
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
    derivatives = []
    for index, ref in enumerate(refs):
        dq = ref["contract"] @ direction
        local = 2 * ref["curvature"] * dq
        df = expand(ref["expand"], local)
        if alpha == 0:
            df[::2] = 0
        actual_f, _ = forces.forcesvirs_4th_order(index)
        np.testing.assert_allclose(actual_f, df, rtol=1e-7, atol=2e-8)
        weight = ref["component"].weight * sum(ref["component"].mts_weights)
        derivatives.append(weight * df)
    np.testing.assert_allclose(
        forces.f_4th_order, sum(derivatives), rtol=1e-7, atol=2e-8
    )


@pytest.mark.parametrize("weight", [0.0, 0.3, -0.2])
def test_weighted_centroid_force_and_saved_state(simulation, weight):
    root = ET.fromstring(harmonic_input(nbeads=4))
    forces = root.find("system/forces")
    ET.SubElement(
        forces, "force", forcefield="harmonic", nbeads="1", weight=str(weight)
    )
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        s = sim.syslist[0]
        q = s.beads.q.copy()
        centroid = q.mean(axis=0)
        expected_force = -0.5 * (q + weight * centroid)
        expected_pots = 0.25 * (np.sum(q * q, axis=1) + weight * np.sum(centroid**2))
        np.testing.assert_allclose(s.forces.f, expected_force)
        np.testing.assert_allclose(s.forces.pots, expected_pots)
        snapshot = s.forces.dump_state()
        for interpolate in [False, True]:
            count = 4 if interpolate else 1
            for weighted in [False, True]:
                w = weight if weighted else 1.0
                np.testing.assert_allclose(
                    s.forces.forces_component(1, weighted, interpolate),
                    np.tile(-0.5 * w * centroid, (count, 1)),
                )
                np.testing.assert_allclose(
                    s.forces.pots_component(1, weighted, interpolate),
                    np.full(count, 0.25 * w * np.sum(centroid**2)),
                )
                np.testing.assert_array_equal(
                    s.forces.virs_component(1, weighted, interpolate),
                    np.zeros((count, 3, 3)),
                )
        s.beads.q_[:] = q * 2
        np.testing.assert_allclose(s.forces.f, 2 * expected_force)
        s.beads.q_[:] = q
        s.forces.load_state(snapshot)
        np.testing.assert_allclose(s.forces.f, expected_force)
        # with zero weight, the component was not computed when the state was saved
        np.testing.assert_allclose(
            s.forces.forces_component(1, False, False), [-0.5 * centroid]
        )
        s.beads.q_[:] = q * 3
        np.testing.assert_allclose(s.forces.f, 3 * expected_force)
        np.testing.assert_allclose(s.forces.pots, 9 * expected_pots)
        with pytest.raises(ValueError, match="components"):
            s.forces.load_state(snapshot[:1])


@pytest.fixture
def force_system():
    """Real force aggregation with a small analytic harmonic provider, including virials and extras."""

    class Harmonic(FFEval):
        def evaluate(self, request):
            q = request["pos"]
            energy = 0.5 * np.dot(q, q)
            extras = {"energy": energy, "label": self.name, "raw": self.name}
            if self.extra_values is not None:
                extras = self.extra_values(q)
            request["result"] = [energy, -q, -np.outer(q, q), extras]
            request["t_finished"] = time.time()
            request["status"] = "Done"
            request._event_done.set()

    def create(components, nbeads=4, extra_values=None):
        # Match the input schema's default: one outer MTS level with unit weight.
        for component in components:
            if len(component.mts_weights) == 0:
                component.mts_weights = np.array([1.0])
        beads = Beads(1, nbeads)
        beads.m_[:] = 2
        beads.names_[:] = ["H"]
        beads.q_[:] = np.arange(1, 3 * nbeads + 1).reshape(nbeads, 3) / 10
        cell = Cell(np.eye(3) * 10)
        providers = {
            component.ffield: Harmonic(
                name=component.ffield, threaded=False, dopbc=False
            )
            for component in components
        }
        for provider in providers.values():
            provider.extra_values = extra_values
        forces = Forces()
        forces.bind(beads, cell, components, providers, [], None)
        return beads, cell, forces

    return create


@pytest.mark.parametrize("contracted", [1, 4])
@pytest.mark.parametrize("weight", [0.0, -0.5, 2.0])
def test_extras_and_virials_follow_component_weights(force_system, contracted, weight):

    beads, _, forces = force_system(
        [
            ForceComponent("first", interpolate_extras=["energy"]),
            ForceComponent(
                "second",
                nbeads=contracted,
                weight=weight,
                interpolate_extras=["energy"],
            ),
        ]
    )
    q = beads.q.copy()
    other = np.tile(q.mean(axis=0), (4, 1)) if contracted == 1 else q
    expected_energy = 0.5 * (
        np.sum(q * q, axis=1) + weight * np.sum(other * other, axis=1)
    )
    expected_virial = np.array(
        [-np.triu(np.outer(a, a) + weight * np.outer(b, b)) for a, b in zip(q, other)]
    )
    np.testing.assert_allclose(forces.extras["energy"], expected_energy)
    np.testing.assert_allclose(forces.pots, expected_energy)
    np.testing.assert_allclose(forces.f, -q - weight * other)
    np.testing.assert_allclose(forces.virs, expected_virial)
    np.testing.assert_allclose(forces.vir, expected_virial.sum(axis=0))
    # Non-numerical extras are not interpolated; full-bead raw strings concatenate.
    full_second = contracted == 4 and weight != 0
    assert forces.extras["raw"] == ["firstsecond" if full_second else "first"] * 4
    assert forces.extras["label"] == ["second" if full_second else "first"] * 4
    assert forces.extras_component(0)["label"] == ["first"] * 4
    assert forces.extras_component(0)["raw"] == ["first"] * 4
    np.testing.assert_allclose(
        forces.mforces[0].vir, np.sum([-np.triu(np.outer(a, a)) for a in q], axis=0)
    )


def test_extras_follow_weight_changes(force_system):
    """The combined extras are recomputed when the weight of a component changes"""
    beads, _, forces = force_system(
        [
            ForceComponent("first", interpolate_extras=["energy"]),
            ForceComponent("second", interpolate_extras=["energy"]),
        ]
    )
    energy = 0.5 * np.sum(beads.q**2, axis=1)
    np.testing.assert_allclose(forces.extras["energy"], 2 * energy)
    forces.mforces[1].weight = 0.5
    np.testing.assert_allclose(forces.extras["energy"], 1.5 * energy)
    forces.mforces[1].weight = 0.0
    np.testing.assert_allclose(forces.extras["energy"], energy)
    assert forces.extras["label"] == ["first"] * 4


@pytest.mark.parametrize("scaling", [0.0, -0.5, 2.0])
def test_scaled_component_updates_energy_force_and_virial(force_system, scaling):

    beads, _, forces = force_system([ForceComponent("harmonic")])
    scaled = ScaledForceComponent(forces.mforces[0], scaling=scaling)
    for displacement in [1.0, 2.0]:
        beads.q_[:] *= displacement
        q = beads.q.copy()
        np.testing.assert_allclose(scaled.f, -scaling * q)
        np.testing.assert_allclose(scaled.pots, scaling * 0.5 * np.sum(q * q, axis=1))
        np.testing.assert_allclose(
            scaled.vir, -scaling * np.sum([np.triu(np.outer(a, a)) for a in q], axis=0)
        )
        assert bool(scaled.extras) == (scaling != 0)
    scaled.scaling = 3
    np.testing.assert_allclose(scaled.f, -3 * beads.q)


@pytest.mark.parametrize("bad", ["missing", "nonnumeric", "extra-key"])
def test_inconsistent_extras_are_rejected(force_system, bad):

    def extras(q):
        if bad == "missing":
            return {}
        if bad == "nonnumeric":
            return {"energy": "not a number"}
        return {"energy": 1.0, **({"other": 2.0} if q[0] > 0.2 else {})}

    _, _, forces = force_system(
        [ForceComponent("harmonic", interpolate_extras=["energy"])], extra_values=extras
    )
    if bad == "nonnumeric":
        with pytest.raises(Exception, match="has to be numerical"):
            forces.mforces[0].extras
    else:
        with pytest.raises(KeyError):
            forces.mforces[0].extras


@pytest.mark.parametrize("zero_force", [False, True])
@pytest.mark.parametrize("alpha", [0.0, 0.3])
@pytest.mark.parametrize("epsilon", [-1e-4, 1e-4])
def test_suzuki_chin_finite_difference_modes(force_system, alpha, epsilon, zero_force):
    """Forward and centered displacements, including the zero-force fallback."""
    beads, _, forces = force_system([ForceComponent("harmonic", epsilon=epsilon)])
    forces.alpha = alpha
    forces.omegan2 = 2.0
    if zero_force:
        beads.q_[:] = 0
    q = beads.q.copy()
    expected = -2 * q / beads.m3  # -gradient of |F|^2/m for U=|q|^2/2
    if alpha == 0:
        expected[::2] = 0
    correction, _ = forces.forcesvirs_4th_order(0)
    np.testing.assert_allclose(correction, expected, rtol=1e-10, atol=1e-12)
    np.testing.assert_array_equal(beads.q, q)


@pytest.mark.parametrize("extras", [None, {"label": ["saved"] * 4}])
def test_manual_force_transfer_invalidated_by_move(force_system, extras):

    beads, _, forces = force_system([ForceComponent("harmonic")])
    q = beads.q.copy()
    virial = np.diag([1.0, 2.0, 3.0])
    forces.transfer_forces_manual(
        [q],
        [np.arange(4.0)],
        [np.full((4, 3), 7.0)],
        None if extras is None else [extras],
        vir=virial,
    )
    np.testing.assert_array_equal(forces.f, 7)
    np.testing.assert_array_equal(forces.pots, np.arange(4.0))
    np.testing.assert_array_equal(forces.vir, 4 * virial)
    if extras:
        assert forces.extras["label"] == ["saved"] * 4
    beads.q_[:] = q * 2
    np.testing.assert_allclose(forces.f, -2 * q)
    np.testing.assert_allclose(forces.pots, 2 * np.sum(q * q, axis=1))


def test_clone_and_force_transfer_preserve_independent_state(force_system):

    beads, cell, forces = force_system([ForceComponent("harmonic")])
    q = beads.q.copy()
    expected = forces.f.copy()
    clone_beads = beads.clone()
    clone = forces.clone(clone_beads, cell.clone())
    clone.transfer_forces(forces)
    np.testing.assert_array_equal(clone.f, expected)
    clone_beads.q_[:] *= 2
    np.testing.assert_allclose(clone.f, 2 * expected)
    np.testing.assert_array_equal(beads.q, q)
    np.testing.assert_array_equal(forces.f, expected)
    snapshot = forces.dump_state()
    snapshot[0][1].pop()
    with pytest.raises(ValueError, match="numbers of beads"):
        forces.load_state(snapshot)


@pytest.mark.parametrize(
    "bad",
    ["unbound", "beads", "cell", "component-beads", "unknown-provider", "mts-levels"],
)
def test_force_binding_rejects_incompatible_objects(force_system, bad):

    beads, cell, forces = force_system([ForceComponent("harmonic")])
    if bad == "unbound":
        with pytest.raises(ValueError, match="not yet been bound"):
            Forces().clone(beads, cell)
    elif bad in ["beads", "cell"]:
        with pytest.raises(ValueError, match="argument"):
            forces.clone(
                None if bad == "beads" else beads, None if bad == "cell" else cell
            )
    elif bad == "mts-levels":
        with pytest.raises(ValueError, match="mts_weights"):
            force_system(
                [
                    ForceComponent("first", mts_weights=[1]),
                    ForceComponent("second", mts_weights=[1, 0]),
                ]
            )
    else:
        component = ForceComponent(
            "missing" if bad == "unknown-provider" else "harmonic",
            nbeads=3 if bad == "component-beads" else 4,
        )
        with pytest.raises(ValueError, match="numbers of beads|not in the forcefields"):
            component.bind(beads, cell, forces.ff, None)


@pytest.mark.parametrize("invalid", ["odd-beads", "forward-contracted"])
def test_suzuki_chin_rejects_unsupported_bead_discretization(force_system, invalid):

    _, _, forces = force_system(
        [
            ForceComponent(
                "harmonic",
                nbeads=1 if invalid == "forward-contracted" else 0,
                epsilon=1e-4,
            )
        ],
        nbeads=3 if invalid == "odd-beads" else 4,
    )
    with pytest.raises(SystemExit):
        forces.forcesvirs_4th_order(0)
    if invalid == "odd-beads":
        with pytest.raises(SystemExit):
            forces.get_potssc()


def test_cancelled_bead_request_cannot_return_force_data(force_system, monkeypatch):

    _, _, forces = force_system([ForceComponent("harmonic")])
    bead_force = forces.mforces[0]._forces[0]
    bead_force.queue()
    bead_force.request["status"] = "Exit"
    messages = []
    monkeypatch.setattr(
        force_module,
        "softexit",
        SimpleNamespace(
            triggered=False,
            exiting=False,
            trigger=lambda **kw: messages.append(kw["message"]),
        ),
    )
    try:
        with pytest.raises(SystemExit):
            bead_force.get_all()
        assert len(messages) == 1
        assert "cannot return" in messages[0]
    finally:
        bead_force.ff.release(bead_force.request)


@pytest.mark.parametrize("readers", [1, 2])
def test_bead_force_waits_for_evaluation_and_releases_request(
    force_system, monkeypatch, readers
):

    _, _, forces = force_system([ForceComponent("harmonic")])
    bead_force = forces.mforces[0]._forces[0]
    provider = bead_force.ff
    # Queue first; evaluation starts only after the readers wait.
    provider.threaded = True
    bead_force.queue()
    request = bead_force.request
    waiting = Barrier(readers + 1)
    event_wait = request._event_done.wait

    def wait_for_result(timeout=None):
        waiting.wait(timeout=2)
        return event_wait(timeout)

    monkeypatch.setattr(request._event_done, "wait", wait_for_result)
    monkeypatch.setattr(
        force_module,
        "softexit",
        SimpleNamespace(triggered=False, exiting=False, trigger=lambda **kw: None),
    )
    with ThreadPoolExecutor(max_workers=readers) as executor:
        results = [executor.submit(bead_force.get_all) for _ in range(readers)]
        try:
            waiting.wait(timeout=2)
            assert all(not result.done() for result in results)
            provider.poll()
            for result in results:
                energy, force, virial, extra = result.result(timeout=2)
                assert energy == pytest.approx(0.07)
                np.testing.assert_allclose(force, [-0.1, -0.2, -0.3])
            assert not provider.requests
            assert bead_force.request is None
        finally:
            if any(not result.done() for result in results):
                request["status"] = "Exit"
                request._event_done.set()


def test_unqueued_bead_force_evaluates_and_releases_request(force_system):

    _, _, forces = force_system([ForceComponent("harmonic")])
    bead_force = forces.mforces[0]._forces[0]
    assert bead_force.request is None
    energy, force, virial, extras = bead_force.get_all()
    assert energy == pytest.approx(0.07)
    np.testing.assert_allclose(force, [-0.1, -0.2, -0.3])
    assert bead_force.request is None
    assert not bead_force.ff.requests

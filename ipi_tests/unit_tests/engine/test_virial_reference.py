"""Virials from -d(complete energy)/dstrain at fixed fractional coordinates.

Kapil et al., JCTC 15, 3237 (2019), DOI 10.1021/acs.jctc.8b01297,
Eq. 9 and SI S1.26--S1.31. The paper uses Xi=-W, opposite to i-PI.
The convex quartic model also has explicit volume dependence in its forces.
"""

from types import SimpleNamespace

import numpy as np
import pytest

from ipi.engine.barostats import BaroBZP, BaroSCBZP, BaroMTK, BaroRGB, BaroSCR
from ipi.engine.forcefields import FFDirect
from ipi.engine.forces import ForceComponent, Forces, ScaledForceComponent
from ipi.utils.prng import Random
from ipi_tests.unit_tests.engine.mathematical_references import (
    cell_strain_virial,
    sc_weights,
    volume_polynomial,
)

SC_VIRIAL_BUG = pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="BUG: SC virial omits the force-product term required by the complete energy's cell-strain derivative (JCTC 2019, SI S1.31)",
)


class VolumePolynomial:
    """Driver differentiates U directly; references use polynomial calculus."""

    def __init__(self, a, b, coupling):
        self.a, self.b, self.coupling = a, b, coupling

    def __call__(self, cell, positions):
        x = positions.reshape(-1)
        ratio = np.linalg.det(cell) / 120
        a = self.a * (1 + self.coupling * (ratio - 1))
        energy = np.sum(a * x**2 + self.b * x**4)
        force = -2 * a * x - 4 * self.b * x**3
        virial = positions.T @ force.reshape(-1, 3)
        virial -= np.eye(3) * self.coupling * ratio * np.sum(self.a * x**2)
        return energy, force.reshape(positions.shape), virial, ""


@pytest.fixture(params=[(4, False), (6, True)], ids=["single", "weighted"])
def virial_case(reference_system, request):
    """Positive weighted components keep the complete potential convex."""

    def create(alpha=0.3, coupling=0.6):
        count, multiple = request.param
        s, r = reference_system(count=count, alpha=alpha)
        specs = [(1.0, [1.0], r.a, r.b)]
        if multiple:
            specs = [
                (0.7, [0.2, 0.4], r.a, r.b),
                (0.3, [0.8, 0.6], 1.4 * r.a, 1.7 * r.b),
            ]
        providers, components = {}, []
        for index, (weight, mts, a, b) in enumerate(specs):
            name = f"polynomial{index}"
            provider = FFDirect(pes="harmonic", pars={"k1": 1.0})
            provider.driver = VolumePolynomial(a, b, coupling)
            providers[name] = provider
            components.append(
                ForceComponent(name, weight=weight, mts_weights=mts, epsilon=-1e-5)
            )
        forces = Forces()
        forces.bind(s.beads, s.cell, components, providers, [], None)
        forces.alpha = alpha
        forces.omegan2 = (count * 0.4) ** 2
        s.properties.forces = s.forces = forces
        a = sum(w * sum(mts) * a for w, mts, a, b in specs)
        b = sum(w * sum(mts) * b for w, mts, a, b in specs)
        mass = np.repeat([2.0, 3.0], 3)
        weights, coefficients = sc_weights(count, alpha, 0.4)

        def model(q, h, a=a, b=b):
            return volume_polynomial(q, h, a, b, coupling)

        def energies(q, h):
            energy, force, _, _ = model(q, h)
            square = np.sum(force**2 / mass, axis=1)
            return np.stack(
                [
                    energy,
                    (weights - 1) * energy,
                    coefficients * square,
                    square,
                    weights * energy + coefficients * square,
                ],
                axis=1,
            )

        return s, SimpleNamespace(
            **vars(r),
            specs=specs,
            coupling=coupling,
            model=model,
            energies=energies,
            weights=weights,
            coefficients=coefficients,
            mass=mass,
        )

    return create


def close(actual, expected):
    # Five-point strain differences have ~1e-11 absolute numerical error here;
    # i-PI's internal centered force displacement adds O(epsilon^2) error.
    np.testing.assert_allclose(actual, expected, rtol=2e-7, atol=2e-9)


def vector(tensor):
    return tensor[[0, 1, 2, 0, 0, 1], [0, 1, 2, 1, 2, 2]]


def kinetic(r, force, centroid=True, triangular=True):
    delta = (r.q - r.q.mean(axis=0)).reshape(-1, 3)
    result = -delta.T @ force.reshape(-1, 3)
    if centroid:
        pc = r.p.mean(axis=0).reshape(2, 3)
        result += np.diag(r.count * np.sum(pc**2 / r.masses[:, None], axis=0))
    return np.triu(result) if triangular else result


def bind_barostat(s, r, sc, bias=None):
    barostat = (BaroSCBZP if sc else BaroBZP)(
        dt=0.02, temp=r.count * 0.4, tau=2.0, pext=0.001, p=0.17
    )
    barostat.bind(
        s.beads,
        s.nm,
        s.cell,
        s.forces,
        bias=bias,
        prng=Random(seed=17),
        nmts=len(r.specs[0][1]),
    )
    return barostat


@pytest.mark.parametrize("coupling", [0.0, 0.6])
def test_ordinary_virials_are_complete_cell_derivatives(virial_case, coupling):
    s, r = virial_case(coupling=coupling)
    # Change volume and shear together, preserving an upper triangular cell.
    for transform in [
        np.eye(3),
        np.array([[1.08, 0, 0], [0.04, 0.94, 0], [0, 0.03, 1.02]]),
    ]:
        q = (r.q.reshape(-1, 3) @ transform).reshape(r.q.shape)
        h = transform.T @ r.h
        s.beads._q[:] = q
        s.cell._h[:] = h
        total = cell_strain_virial(lambda x, cell: r.energies(x, cell)[:, 0], q, h)
        close(s.forces.virs, total)
        close(s.forces.vir_combine(), total)
        close(s.forces.get_vir(), total.sum(axis=0))
        close(s.forces.pots, r.model(q, h)[0])
        close(s.forces.f, r.model(q, h)[1])
        for index, (w, mts, a, b) in enumerate(r.specs):
            raw = cell_strain_virial(
                lambda x, cell: volume_polynomial(x, cell, a, b, coupling)[0], q, h
            )
            component = s.forces.mforces[index]
            close(component.vir_gather(), raw)
            close(component.get_vir(), raw.sum(axis=0))
            close(
                ScaledForceComponent(component, scaling=-0.4).get_vir(),
                -0.4 * raw.sum(axis=0),
            )
            for bead, expected in zip(component._forces, raw):
                close(bead.get_vir(), expected)
            for weighted in [False, True]:
                for interpolate in [False, True]:
                    close(
                        s.forces.virs_component(index, weighted, interpolate),
                        raw * (w if weighted else 1),
                    )


@pytest.mark.parametrize("virial_case", [(6, True)], indirect=True)
def test_sc_complete_energy(virial_case):
    # The complete SC energy is the reference differentiated in the tests below.
    s, r = virial_case()
    for scale in [0.93, 1.0, 1.07]:
        q, h = scale * r.q, scale * r.h
        s.beads._q[:] = q
        s.cell._h[:] = h
        close(s.forces.pots + s.forces.potssc, r.energies(q, h)[:, 4])


@pytest.mark.parametrize("coupling", [0.0, 0.6])
def test_published_sc_virial_formula(virial_case, coupling):
    s, r = virial_case(coupling=coupling)
    force = r.model(r.q, r.h)[1]
    direction = force / r.mass
    step = 1e-3
    values = []
    for displacement in [-2 * step, -step, step, 2 * step]:
        s.beads._q[:] = r.q + displacement * direction
        # Later force evaluations update the dependency array in place.
        values.append(np.asarray(s.forces.virs).copy())
    s.beads._q[:] = r.q
    derivative = (values[0] - 8 * values[1] + 8 * values[2] - values[3]) / (12 * step)
    product = np.array(
        [np.triu(f.reshape(2, 3).T @ (f / r.mass).reshape(2, 3)) for f in force]
    )
    # SI S1.31 after changing Xi=-W: W_square = 2 F outer F/m - 2 D_(F/m) W.
    published = 2 * product - 2 * derivative
    expected = cell_strain_virial(r.energies, r.q, r.h)[:, 3]
    close(published, expected)


@SC_VIRIAL_BUG
@pytest.mark.parametrize("alpha", [0.0, 0.3, 1.0])
@pytest.mark.parametrize("coupling", [0.0, 0.6])
@pytest.mark.parametrize("epsilon", [-1e-5, 1e-6], ids=["centered", "forward"])
@pytest.mark.parametrize(
    "virial_case,component_index",
    [((4, False), 0), ((6, True), 0), ((6, True), 1)],
    indirect=["virial_case"],
)
def test_sc_component_virial_differentiates_force_square(
    virial_case, alpha, coupling, component_index, epsilon
):
    s, r = virial_case(alpha, coupling)
    s.forces.mforces[component_index].epsilon = epsilon
    _, _, a, b = r.specs[component_index]
    fixed = r.model(r.q, r.h)[1] / r.mass

    def cross(q, h):
        force = volume_polynomial(q, h, a, b, coupling)[1]
        return 2 * np.sum(fixed * force, axis=1)

    expected = cell_strain_virial(cross, r.q, r.h)
    if alpha == 0:
        expected[::2] = 0  # These internal entries are unused when c_even=0.
    close(s.forces.forcesvirs_4th_order(component_index)[1], expected)


@SC_VIRIAL_BUG
@pytest.mark.parametrize("virial_case", [(6, True)], indirect=True)
@pytest.mark.parametrize("quantity", ["combine", "cached"])
def test_sc_combined_raw_virial_differentiates_force_square(virial_case, quantity):
    # The component calculation above carries the detailed parameter checks.
    # Here a weighted case checks summation and the cached accessor separately.
    s, r = virial_case()
    expected = cell_strain_virial(r.energies, r.q, r.h)[:, 3]
    actual = (
        s.forces.fvir_4th_order_combine()[1]
        if quantity == "combine"
        else s.forces.vir_4th_order
    )
    close(actual, expected)


@pytest.mark.parametrize("virial_case", [(6, True)], indirect=True)
def test_sc_linear_virial_differentiates_linear_energy(virial_case):
    s, r = virial_case()
    expected = cell_strain_virial(r.energies, r.q, r.h)[:, 1]
    close(s.forces.get_virssc_part_1(), expected)


@SC_VIRIAL_BUG
@pytest.mark.parametrize("virial_case", [(6, True)], indirect=True)
@pytest.mark.parametrize("part", ["virssc_part_2", "virssc"])
def test_sc_corrected_virial_differentiates_complete_energy(virial_case, part):
    s, r = virial_case()
    strain = cell_strain_virial(r.energies, r.q, r.h)
    expected = strain[:, 2] if part == "virssc_part_2" else strain[:, 1] + strain[:, 2]
    close(getattr(s.forces, part), expected)


@pytest.mark.parametrize(
    "key",
    [
        "stress_md",
        "pressure_md",
        "virial_md",
        "stress_cv",
        "pressure_cv",
        "virial_cv",
        "vir_component",
        "vir_component_raw",
    ],
)
def test_pressure_outputs_include_explicit_cell_derivatives(virial_case, key):
    s, r = virial_case()
    _, force, _, _ = r.model(r.q, r.h)
    strain = cell_strain_virial(r.energies, r.q, r.h)
    virial = strain[:, 0].sum(axis=0)
    cv = kinetic(r, force)
    md = sum(p.reshape(2, 3).T @ (p / r.mass).reshape(2, 3) for p in r.p)
    volume = np.linalg.det(r.h)
    expected = dict(
        stress_md=vector((virial + md) / volume),
        pressure_md=np.trace(virial + md) / (3 * volume),
        virial_md=vector(virial / volume),
        stress_cv=vector((virial + cv) / (r.count * volume)),
        pressure_cv=np.trace(virial + cv) / (3 * r.count * volume),
        virial_cv=vector(virial / (r.count * volume)),
    )
    if key.startswith("vir_component"):
        w, _, a, b = r.specs[-1]
        raw = cell_strain_virial(
            lambda q, h: volume_polynomial(q, h, a, b, r.coupling)[0], r.q, r.h
        )
        raw *= 1 if key.endswith("_raw") else w
        for bead in [-1, 1]:
            close(
                s.properties[f"{key}({len(r.specs) - 1};{bead})"][0],
                vector(raw.mean(axis=0) if bead < 0 else raw[bead]) / volume,
            )
    else:
        close(s.properties[key][0], expected[key])


@SC_VIRIAL_BUG
@pytest.mark.parametrize("virial_case", [(6, True)], indirect=True)
@pytest.mark.parametrize("key", ["pressure_tdsc", "vir_tdsc"])
def test_sc_pressure_outputs_use_complete_energy_strain(virial_case, key):
    # One weighted, volume-dependent example checks each output's trace and
    # volume/bead factors; component parameter limits are tested above.
    s, r = virial_case()
    virial = cell_strain_virial(r.energies, r.q, r.h)[:, 4].sum(axis=0)
    if key == "pressure_tdsc":
        _, force, curvature, _ = r.model(r.q, r.h)
        sc_force = (
            r.weights[:, None] * force
            + 2 * r.coefficients[:, None] * curvature * force / r.mass
        )
        virial += kinetic(r, sc_force)
    expected = np.trace(virial) / (3 * r.count * np.linalg.det(r.h))
    close(s.properties[key][0], expected)


@pytest.mark.parametrize(
    "virial_case,sc",
    [
        pytest.param((4, False), False, id="single-ordinary"),
        pytest.param((6, True), False, id="weighted-ordinary"),
        pytest.param((6, True), True, id="weighted-SC-linear"),
    ],
    indirect=["virial_case"],
)
@pytest.mark.parametrize("biased", [False, True])
def test_mts_barostat_stress_uses_complete_cell_derivatives(virial_case, sc, biased):
    s, r = virial_case()
    bias = None
    bias_force, bias_virial = np.zeros_like(r.q), np.zeros((3, 3))
    if biased:
        provider = FFDirect(pes="harmonic", pars={"k1": 1.0})
        provider.driver = VolumePolynomial(0.12 * r.a, 0.2 * r.b, 0.3)
        bias = Forces()
        bias.bind(
            s.beads,
            s.cell,
            [ForceComponent("bias", mts_weights=[1.0])],
            {"bias": provider},
            [],
            None,
        )
        bias_force = volume_polynomial(r.q, r.h, 0.12 * r.a, 0.2 * r.b, 0.3)[1]
        bias_virial = cell_strain_virial(
            lambda q, h: volume_polynomial(q, h, 0.12 * r.a, 0.2 * r.b, 0.3)[0],
            r.q,
            r.h,
        ).sum(axis=0)
    barostat = bind_barostat(s, r, sc, bias)
    for level in range(len(r.specs[0][1])):
        a = sum(w * mts[level] * a for w, mts, a, b in r.specs)
        b = sum(w * mts[level] * b for w, mts, a, b in r.specs)
        force = volume_polynomial(r.q, r.h, a, b, r.coupling)[1]
        virial = cell_strain_virial(
            lambda q, h: volume_polynomial(q, h, a, b, r.coupling)[0], r.q, r.h
        )
        close(s.forces.mts_forces[level].get_virs_mts(), virial)
        close(s.forces.mts_forces[level].get_vir_mts(), virial.sum(axis=0))
        weights = r.weights if sc else np.ones(r.count)
        level_force = weights[:, None] * force + (bias_force if level == 0 else 0)
        kst = kinetic(
            r,
            level_force,
            centroid=level == len(r.specs[0][1]) - 1,
            triangular=sc,
        )
        expected = (
            kst
            + np.sum(weights[:, None, None] * virial, axis=0)
            + (bias_virial if level == 0 else 0)
        ) / np.linalg.det(r.h)
        actual_k = barostat.kstress_mts_sc(level) if sc else barostat.kstress_mts(level)
        actual_s = barostat.stress_mts_sc(level) if sc else barostat.stress_mts(level)
        close(actual_k, kst)
        close(actual_s, expected)


@pytest.mark.parametrize("virial_case", [(6, True)], indirect=True)
def test_sc_barostat_kinetic_correction(virial_case):
    s, r = virial_case()
    barostat = bind_barostat(s, r, True)
    _, force, curvature, _ = r.model(r.q, r.h)
    force2 = 2 * r.coefficients[:, None] * curvature * force / r.mass
    close(barostat.get_kstress_sc(), kinetic(r, force2, centroid=False))


@SC_VIRIAL_BUG
@pytest.mark.parametrize("virial_case", [(6, True)], indirect=True)
@pytest.mark.parametrize("quantity", ["stress", "piston_momentum"])
def test_sc_barostat_high_order_stress_and_momentum(virial_case, quantity):
    s, r = virial_case()
    barostat = bind_barostat(s, r, True)
    _, force, curvature, _ = r.model(r.q, r.h)
    force2 = 2 * r.coefficients[:, None] * curvature * force / r.mass
    virial2 = cell_strain_virial(r.energies, r.q, r.h)[:, 2].sum(axis=0)
    kst = kinetic(r, force2, centroid=False)
    stress = (kst + virial2) / np.linalg.det(r.h)
    if quantity == "stress":
        close(barostat.get_stress_sc(), stress)
    else:
        # SI S1.32--S1.36: dt/2 times V trace(stress), with unnormalized bead sum.
        expected = 0.17 + 0.02 / 2 * np.linalg.det(r.h) * np.trace(stress)
        barostat.pscstep()
        close(barostat.p, [expected])


@pytest.mark.parametrize(
    "virial_case,kind",
    [
        ((4, False), "BZP"),
        ((6, True), "BZP"),
        ((6, True), "SCBZP"),
        ((4, False), "MTK"),
        ((6, True), "MTK"),
        ((4, False), "RGB"),
        ((6, True), "RGB"),
    ],
    indirect=["virial_case"],
)
def test_barostat_momentum_receives_complete_virial_contribution(virial_case, kind):
    # At V=120 both models have identical U, forces and momenta. Their virials
    # differ only by the explicit cell derivative. All other kick terms cancel.
    # Hamilton's equation gives delta(p)=dt*trace(delta W) for an isotropic
    # logarithmic cell coordinate, or dt*delta W for the triangular cell tensor.
    momenta, derivatives = [], []
    for coupling in [0.0, 0.6]:
        s, r = virial_case(coupling=coupling)
        kwargs = dict(dt=0.02, temp=r.count * 0.4, tau=2.0)
        if kind == "RGB":
            barostat = BaroRGB(
                **kwargs,
                h0=s.cell.clone(),
                stressext=0.001 * np.eye(3),
                p=0.17 * np.eye(3),
            )
        elif kind == "MTK":
            barostat = BaroMTK(**kwargs, pext=0.001, p=0.17 * np.eye(3))
        else:
            barostat = (BaroBZP if kind == "BZP" else BaroSCBZP)(
                **kwargs, pext=0.001, p=0.17
            )
        levels = len(r.specs[0][1])
        barostat.bind(
            s.beads, s.nm, s.cell, s.forces, prng=Random(seed=17), nmts=levels
        )
        barostat._pdt[:] = np.array([0.01 * (level + 1) for level in range(levels)])
        steps, expected_steps = [], []
        cumulative = np.zeros(1 if kind in ["BZP", "SCBZP"] else (3, 3))
        for level in range(levels):
            a = sum(w * mts[level] * a for w, mts, a, b in r.specs)
            b = sum(w * mts[level] * b for w, mts, a, b in r.specs)
            weights = r.weights if kind == "SCBZP" else np.ones(r.count)
            virial = cell_strain_virial(
                lambda q, h: weights * volume_polynomial(q, h, a, b, coupling)[0],
                r.q,
                r.h,
            ).sum(axis=0)
            cumulative += (
                0.01
                * (level + 1)
                * (np.trace(virial) if kind in ["BZP", "SCBZP"] else virial)
            )
            barostat.pstep(level)
            steps.append(np.asarray(barostat.p).copy())
            expected_steps.append(cumulative.copy())
        momenta.append(np.asarray(steps))
        derivatives.append(np.asarray(expected_steps))
    close(momenta[1] - momenta[0], derivatives[1] - derivatives[0])


@pytest.mark.parametrize("coupling", [0.0, 0.6])
def test_scr_pressure_force_volume_derivative(reference_system, coupling):
    # SCR supports classical dynamics only. Bernetti & Bussi, JCP 153,
    # 114107 (2020), Eq. 7: F_lambda=2 sqrt(V)*(Pint-Pext+kBT/(2V)).
    s, r = reference_system(count=1)
    s.forces.ff["harmonic"].driver = VolumePolynomial(r.a, r.b, coupling)
    barostat = BaroSCR(dt=0.02, temp=0.4, tau=2.0, pext=0.001, compressibility=0.3)
    barostat.bind(s.beads, s.nm, s.cell, s.forces, prng=Random(seed=17))
    virial = cell_strain_virial(
        lambda q, h: volume_polynomial(q, h, r.a, r.b, coupling)[0], r.q, r.h
    ).sum(axis=0)
    kinetic_energy = np.sum(r.p**2 / (2 * np.repeat([2.0, 3.0], 3)))
    pressure = (2 * kinetic_energy + np.trace(virial)) / (3 * 120)
    close(
        barostat.get_lambda_force(),
        2 * np.sqrt(120) * (pressure - 0.001 + 0.4 / (2 * 120)),
    )


@pytest.mark.parametrize("reference", [False, True])
def test_force_position_virial_output(virial_case, tmp_path, reference):
    s, r = virial_case()
    origin = np.arange(6) * 0.02 if reference else np.zeros(6)
    step = 1e-4
    # This output is defined as F dot (q-origin)/(2P). Hold the cell fixed;
    # it intentionally excludes the separate explicit cell derivative.
    values = [
        r.model(origin + (1 + eps) * (r.q - origin), r.h)[0]
        for eps in [-2 * step, -step, step, 2 * step]
    ]
    expected = -np.sum(values[0] - 8 * values[1] + 8 * values[2] - values[3]) / (
        12 * step * 2 * r.count
    )
    path = tmp_path / "origin.dat"
    np.savetxt(path, origin)
    close(s.properties.get_fqvirial(str(path) if reference else ""), expected)

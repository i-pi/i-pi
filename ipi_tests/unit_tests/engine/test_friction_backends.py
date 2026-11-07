"""Backend equivalence and native-RNG checkpoint tests for electronic friction.

Matched noise isolates array algebra from backend-specific random generators.
Restart tests separately exercise each backend's actual random-number state.
"""

import json
import numpy as np
import pytest

from ipi.utils import array_backend as ab
from ipi.utils.io import NumpyEncoder
from ipi.engine.outputs import CheckpointOutput
from ipi.scripting import InteractiveSimulation
from ipi_tests.unit_tests.engine.test_friction_restart import make_simulation
from ipi_tests.unit_tests.engine.test_friction_external_meanfield import (
    create,
    load,
    advance,
)


@pytest.fixture(params=["numpy", "torch", "jax"])
def backend(request):
    if request.param != "numpy":
        pytest.importorskip(request.param)
    if request.param == "torch":
        request.getfixturevalue("strict_torch")
    ab.set_backend(request.param)
    yield request.param
    ab.set_backend("numpy")


def snapshot(sim):
    sy = sim.syslist[0]
    f = sy.motion.friction
    if f.variable_friction and f.bath.s is None:
        f.bath._get_non_markovian_jacobian()
    return [
        np.array(ab.to_numpy(a), copy=True)
        for a in (sy.beads.q, sy.beads.p, f.bath.s, f.ediss, f.energy_mf)
    ]


def matched_run(factory, steps=15):
    sim = factory()
    snapshot(sim)
    bath = sim.syslist[0].motion.friction.bath
    bath.s = ab.xp.asarray(
        np.linspace(0.01, 0.1, np.prod(bath.s.shape)).reshape(bath.s.shape)
    )
    rng = np.random.default_rng(374)
    sim.prng.gvec = lambda shape: ab.xp.asarray(rng.standard_normal(shape))
    states = []
    for _ in range(steps):
        advance(sim)
        states.append(snapshot(sim))
    return states


@pytest.mark.parametrize("variable", [False, True])
@pytest.mark.parametrize("splitting", ["fbabf", "symmetric-gle"])
def test_gle_backend_equivalence(tmp_path, backend, variable, splitting):
    def factory():
        return make_simulation(tmp_path, variable, splitting)

    actual = matched_run(factory)
    ab.set_backend("numpy")
    reference = matched_run(factory)
    for a, r in zip(actual, reference):
        for av, rv in zip(a, r):
            np.testing.assert_allclose(av, rv, rtol=2e-11, atol=2e-13)


@pytest.mark.parametrize("fixed", ["", "3,4,5"])
def test_ring_driver_backend_equivalence(tmp_path, backend, fixed):
    actual = matched_run(lambda: create(tmp_path, fixed=fixed))
    ab.set_backend("numpy")
    reference = matched_run(lambda: create(tmp_path, fixed=fixed))
    for a, r in zip(actual, reference):
        for av, rv in zip(a, r):
            np.testing.assert_allclose(av, rv, rtol=2e-11, atol=2e-13)


@pytest.mark.parametrize("kind", ["static", "variable", "ring"])
@pytest.mark.parametrize("split_step", [0, 5])
def test_native_rng_restart(tmp_path, backend, kind, split_step):
    sim = (
        create(tmp_path, fixed="3,4,5")
        if kind == "ring"
        else make_simulation(tmp_path, kind == "variable", "symmetric-gle")
    )
    for _ in range(split_step):
        advance(sim)
    path = tmp_path / "checkpoint.xml"
    output = CheckpointOutput(filename=str(path), stride=1)
    output.bind(sim)
    output.write(store=False)
    rng = json.dumps(sim.prng.state, cls=NumpyEncoder)
    reference = []
    for _ in range(6):
        advance(sim)
        reference.append(snapshot(sim))
    final_rng = json.dumps(sim.prng.state, cls=NumpyEncoder)
    resumed = (
        load(path.read_text())
        if kind == "ring"
        else InteractiveSimulation(path.read_text())
    )
    assert json.dumps(resumed.prng.state, cls=NumpyEncoder) == rng
    for ref in reference:
        advance(resumed)
        for a, r in zip(snapshot(resumed), ref):
            np.testing.assert_allclose(a, r, rtol=2e-12, atol=2e-14)
    assert json.dumps(resumed.prng.state, cls=NumpyEncoder) == final_rng


def make_markovian(tmp_path, variable):
    # Reuse the nonlinear, conservative custom driver from the GLE test.
    from ipi_tests.unit_tests.engine.test_friction_restart import DRIVER

    path = tmp_path / "markov_driver.py"
    path.write_text(DRIVER)
    xml = f"""<simulation verbosity="quiet" floatformat="%.17e">
    <prng><seed>932</seed></prng><total_steps>1000</total_steps>
    <ffdirect name="test" pbc="false"><pes>custom</pes><pes_path>{path}</pes_path></ffdirect>
    <output/><system><beads natoms="1" nbeads="1">
    <q shape="(1,3)">[.2,.3,.4]</q><p shape="(1,3)">[.4,-.1,.3]</p>
    <m shape="(1)">[2]</m><names shape="(1)">[H]</names></beads>
    <cell shape="(3,3)">[10,0,0,0,10,0,0,0,10]</cell>
    <ensemble><temperature>.01</temperature></ensemble>
    <forces><force forcefield="test"/></forces>
    <motion mode="dynamics"><fixcom>False</fixcom><dynamics mode="nve-f" splitting="fbabf">
    <timestep>.02</timestep><friction><bath_mode>markovian</bath_mode><debug_mf_mode>off</debug_mf_mode>
    <variable_friction>{variable}</variable_friction></friction></dynamics></motion>
    </system></simulation>"""
    return InteractiveSimulation(xml)


@pytest.mark.parametrize("variable", [False, True])
def test_markovian_backend_equivalence(tmp_path, backend, variable):
    actual = matched_run(lambda: make_markovian(tmp_path, variable))
    ab.set_backend("numpy")
    reference = matched_run(lambda: make_markovian(tmp_path, variable))
    for a, r in zip(actual, reference):
        for av, rv in zip(a, r):
            np.testing.assert_allclose(av, rv, rtol=2e-11, atol=2e-13)


def test_fixed_substrate_energy_convergence(tmp_path, backend):
    errors = []
    for dt in [0.02, 0.01]:
        sim = create(tmp_path / str(dt), dt=dt, fixed="3,4,5")
        sy = sim.syslist[0]
        f = sy.motion.friction
        f.bath._get_non_markovian_jacobian()
        f.bath.s = ab.xp.full(f.bath.s.shape, 0.1)
        f.bath.os_step = lambda _: None

        def energy():
            return float(
                sy.forces.pot
                + sy.nm.kin
                + sy.nm.vspring
                + f.energy_mf
                + 0.5 * ab.xp.sum(f.bath.s**2)
            )

        initial = energy()
        error = 0.0
        for _ in range(round(1 / dt)):
            advance(sim)
            error = max(error, abs(energy() - initial))
        errors.append(error)
    assert 3.5 < errors[0] / errors[1] < 4.5, errors


@pytest.mark.parametrize("variable", [False, True])
def test_doublewell_backend_equivalence(tmp_path, backend, variable):
    def factory():
        return make_simulation(tmp_path, variable, "symmetric-gle", doublewell=True)

    actual = matched_run(factory, steps=50)
    ab.set_backend("numpy")
    reference = matched_run(factory, steps=50)
    for a, r in zip(actual, reference):
        for av, rv in zip(a, r):
            np.testing.assert_allclose(av, rv, rtol=2e-11, atol=2e-13)


@pytest.mark.parametrize("splitting", ["ofbabfo", "bafofab"])
def test_nvt_f_backend_equivalence(tmp_path, backend, splitting):
    def factory():
        return make_simulation(tmp_path, True, splitting, mode="nvt-f")

    actual = matched_run(factory)
    ab.set_backend("numpy")
    reference = matched_run(factory)
    for a, r in zip(actual, reference):
        for av, rv in zip(a, r):
            np.testing.assert_allclose(av, rv, rtol=2e-11, atol=2e-13)


@pytest.mark.parametrize("drift", [0.0, 1e-200, 0.9, 1.0])
def test_scalar_ou_heat_remains_finite(tmp_path, backend, drift):
    """Strong damping must not divide by an underflowed decay factor."""
    from ipi.engine.friction import _apply_mass_scaled_ou

    p = ab.xp.asarray([2.0, -3.0])
    sm = ab.xp.asarray([2.0, 3.0])
    noise = ab.xp.asarray([0.4, -0.2])
    scale = (1.0 - drift**2) ** 0.5
    actual, heat = _apply_mass_scaled_ou(p, sm, drift, scale, noise)
    expected = drift * np.array([1.0, -1.0]) + scale * np.array([0.4, -0.2])
    np.testing.assert_allclose(ab.to_numpy(actual), expected * [2.0, 3.0])
    np.testing.assert_allclose(float(heat), 0.5 * (2.0 - expected @ expected))

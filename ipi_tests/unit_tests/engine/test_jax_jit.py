"""Compilation of the MD steps with jax.jit (`--jit`)."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

import numpy as np
import pytest

import ipi.utils.array_backend as ab
from ipi_tests.unit_tests.engine.inputs import (
    DUMMY,
    HARMONIC,
    NO_FORCES,
    PROPERTIES,
    TOLERANCE,
    case_input,
    write_input,
)

pytest.importorskip("jax")

# the options of the simulations with the JAX backend
JAX = dict(backend="jax")


@pytest.fixture(autouse=True)
def jax_random(device_random):
    """Draws the random numbers with JAX also with --ipi-backend jax, as the
    numpy draws used then would be fixed in the compiled steps."""


def jit_step(sim):
    """The object compiling the steps of the first system, if any."""

    from ipi.engine.jit import compiled_step

    return compiled_step(sim.syslist[0])


def assert_same_state(a, b):
    for name in PROPERTIES:
        np.testing.assert_allclose(b.properties(name), a.properties(name), **TOLERANCE)
    a, b = a.syslist[0], b.syslist[0]
    for x, y in [(a.beads.q, b.beads.q), (a.beads.p, b.beads.p), (a.cell.h, b.cell.h)]:
        np.testing.assert_allclose(ab.to_numpy(y), ab.to_numpy(x), **TOLERANCE)


@pytest.mark.parametrize(
    "case,forcefield",
    [
        ("nve", DUMMY),
        ("nvt_langevin", NO_FORCES),
        ("nvt_svr", NO_FORCES),
        ("pimd_nvt", DUMMY),
        ("pimd_pile_g", DUMMY),
        ("nvt_gle", DUMMY),
        ("pimd_nm_gle", DUMMY),
        ("pimd_nm_gle_g", DUMMY),
        ("npt_gle_barostat", NO_FORCES),
        ("npt_isotropic", DUMMY),
        ("npt_flexible", NO_FORCES),
        ("pimd_npt", DUMMY),
        ("mts", DUMMY),
        ("nst", DUMMY),
    ],
)
def test_compiled_steps_match_eager_steps(case, forcefield, simulation):
    """The compiled steps give the same trajectory and outputs as eager steps
    (to rounding errors), and draw the same random numbers"""

    xml = case_input(case, forcefield)
    with simulation(xml, **JAX) as eager, simulation(xml, **JAX, jit=True) as compiled:
        for step in range(8):
            eager.run_step(step)
            compiled.run_step(step)
            assert_same_state(eager, compiled)
        assert jit_step(compiled).compiled is not None
        assert compiled.prng.get_state()["jax"] == eager.prng.get_state()["jax"]


def untraceable(self, cell, pos):
    """Zero forces computed with numpy, also from JAX arrays."""

    pos = np.asarray(pos)
    return 0.0, np.zeros_like(pos), np.zeros((3, 3)), ""


@pytest.mark.parametrize("driver", ["numpy", "untraceable"])
def test_forces_of_untraceable_drivers_are_eager(driver, simulation, monkeypatch):
    """Steps that need the forces of a driver working on numpy arrays, or on
    JAX arrays without being traceable by jax.jit, are not compiled, and run
    eagerly"""

    from ipi.pes.dummy import Zero_driver

    if driver == "numpy":
        xml = write_input("nvt_langevin", HARMONIC)
    else:
        monkeypatch.setattr(Zero_driver, "compute_structure", untraceable)
        monkeypatch.setattr(Zero_driver, "traceable", False)
        xml = write_input("nvt_langevin", DUMMY)
    with simulation(xml, **JAX) as eager, simulation(xml, **JAX, jit=True) as compiled:
        for step in range(5):
            eager.run_step(step)
            compiled.run_step(step)
        assert jit_step(compiled).compiled is None
        assert "forcefield 'ff'" in jit_step(compiled).reason
        assert_same_state(eager, compiled)


def test_steps_with_rotations_are_eager(simulation):
    """The steps that need the forces of FFRotations, whose random key cannot
    be traced, run eagerly"""

    xml = case_input("rotations", DUMMY)
    with simulation(xml, **JAX) as eager, simulation(xml, **JAX, jit=True) as compiled:
        for step in range(5):
            eager.run_step(step)
            compiled.run_step(step)
        assert jit_step(compiled).compiled is None
        assert "forcefield 'ff'" in jit_step(compiled).reason
        assert_same_state(eager, compiled)


def test_temperature_changes_reach_the_compiled_step(simulation):
    """The coefficients of the thermostat are arrays, recomputed when the
    temperature changes and passed to the compiled step"""

    xml = write_input("pimd_nvt", NO_FORCES)
    with simulation(xml, **JAX) as eager, simulation(xml, **JAX, jit=True) as compiled:
        for step in range(8):
            if step == 5:
                for sim in (eager, compiled):
                    sim.syslist[0].ensemble.temp *= 2.0
            eager.run_step(step)
            compiled.run_step(step)
        assert_same_state(eager, compiled)


@pytest.mark.parametrize("step", [1, 2, 3])
@pytest.mark.parametrize("case", ["nve", "pimd_nvt"])
def test_state_written_during_the_eager_steps(case, step, simulation):
    """Momenta and positions written between the eager steps, which change the
    objects that a step reads before writing them, do not stop the compilation"""

    xml = write_input(case, DUMMY)
    with simulation(xml, **JAX) as eager, simulation(xml, **JAX, jit=True) as compiled:
        for sim in (eager, compiled):
            for i in range(8):
                if i == step:
                    sim.thermalize_momenta()
                    sim.syslist[0].beads.q = sim.syslist[0].beads.q * 1.001
                sim.run_step(i)
        assert jit_step(compiled).compiled is not None
        assert_same_state(eager, compiled)


def test_changed_input_types_stop_the_compiled_step(select_backend):
    """An input of the compiled step that changes shape or type stops it"""

    import jax.numpy as jnp

    from ipi.utils.depend import depend_array, depend_value
    from ipi.utils.jaxjit import JaxReplay
    from ipi.utils.jit import JitError

    select_backend("jax")
    scale = depend_value(name="scale", value=jnp.asarray(2.0))
    y = depend_array(name="y", value=jnp.zeros(3))
    replay = JaxReplay(lambda: y.set(y.get() + scale.get()))
    for _ in range(4):
        replay()
    assert replay.compiled is not None
    scale.set(jnp.asarray(2, dtype=jnp.int32))
    with pytest.raises(JitError, match="shape or type"):
        replay()


def test_failed_compilation_is_final(simulation, monkeypatch):
    """A step whose tracing fails is not called again, as the objects written
    by the traced call hold tracers: the following calls raise a JitError, and
    the key of the random number generator is restored"""

    import jax

    from ipi.engine.thermostats import ThermoLangevin
    from ipi.utils import depend
    from ipi.utils.jaxjit import _TraceRecorder
    from ipi.utils.jit import JitError

    step = ThermoLangevin.step

    def failing(self):
        if isinstance(depend._recorder, _TraceRecorder):
            raise RuntimeError("tracing failed")
        step(self)

    monkeypatch.setattr(ThermoLangevin, "step", failing)
    with simulation(write_input("nvt_langevin", NO_FORCES), **JAX, jit=True) as sim:
        with pytest.raises(JitError, match="cannot be compiled"):
            for i in range(8):
                sim.run_step(i)
        assert not isinstance(sim.syslist[0].prng.jax_key, jax.core.Tracer)
        for _ in range(2):
            with pytest.raises(JitError, match="has failed"):
                sim.run_step(i + 1)


@pytest.mark.parametrize("kind", ["array", "value"])
def test_host_writes_cannot_be_compiled(kind, select_backend):
    from ipi.utils.depend import depend_array, depend_value
    from ipi.utils.jaxjit import JaxReplay
    from ipi.utils.jit import JitError

    select_backend("jax")
    if kind == "array":
        dep = depend_array(name="host", value=np.zeros(3), host=True)
        replay = JaxReplay(lambda: dep.__setitem__(0, 1.0))
    else:
        dep = depend_value(name="host", value=1.0)
        replay = JaxReplay(lambda: dep.set(2.0))
    with pytest.raises(JitError, match="modifies the host|modifies host"):
        for _ in range(8):
            replay()
    assert replay.compiled is None


def test_points_that_must_run_eagerly_stop_the_tracing(select_backend):
    """A point of the function that must run eagerly in the traced call, but
    not in the observed ones, stops the compilation with its reason"""

    import jax.numpy as jnp

    from ipi.utils import jit
    from ipi.utils.depend import depend_array
    from ipi.utils.jaxjit import JaxReplay
    from ipi.utils.jit import JitError

    class Replay(JaxReplay):
        def split(self, obj):
            if self.inputs is not None:
                self.cannot(f"the point '{obj}' must run eagerly")

    select_backend("jax")
    y = depend_array(name="y", value=jnp.zeros(3))

    def func():
        jit.split("p")
        y.set(y.get() + 1.0)

    replay = Replay(func)
    with pytest.raises(JitError, match="the point 'p' must run eagerly"):
        for _ in range(8):
            replay()
    assert replay.compiled is None


def test_fixed_values_are_compared_as_arrays(select_backend):
    """A value fixed in the compiled step that holds arrays can be replaced by
    an equal one, but not by a different one"""

    import jax.numpy as jnp

    from ipi.utils.depend import depend_array, depend_value
    from ipi.utils.jaxjit import JaxReplay
    from ipi.utils.jit import JitError

    select_backend("jax")
    fixed = depend_value(name="fixed", value=(np.ones(2), np.ones(2)))
    y = depend_array(name="y", value=jnp.zeros(2))
    replay = JaxReplay(lambda: y.set(y.get() + fixed.get()[1]))
    for _ in range(4):
        replay()
    assert replay.compiled is not None
    fixed.set((np.ones(2), np.ones(2)))
    replay()
    fixed.set((np.ones(2), 2 * np.ones(2)))
    with pytest.raises(JitError, match="cannot change"):
        replay()


def test_steps_that_read_other_objects_at_each_call_are_eager(select_backend):
    """A function that reads different depend objects at every call, e.g. as
    the state is written differently between the calls, runs eagerly"""

    import jax.numpy as jnp

    from ipi.utils.depend import depend_array
    from ipi.utils.jaxjit import JaxReplay

    select_backend("jax")
    a = depend_array(name="a", value=jnp.ones(2))
    b = depend_array(name="b", value=jnp.ones(2))
    y = depend_array(name="y", value=jnp.zeros(2))
    calls = []

    def func():
        calls.append(None)
        y.set(y.get() + (a if len(calls) % 2 else b).get())

    replay = JaxReplay(func)
    for _ in range(20):
        replay()
    assert replay.compiled is None
    assert "change at every step" in replay.reason
    assert float(y.get()[0]) == 20.0


def test_read_only_simulations_are_not_compiled(tmp_path, monkeypatch, select_backend):
    """Simulations bound to be read only are not stepped, so their steps are
    not replaced"""

    from ipi.engine.simulation import Simulation

    monkeypatch.chdir(tmp_path)
    select_backend("jax", jit=True)
    xml = write_input("nvt_langevin", NO_FORCES)
    sim = Simulation.load_from_xml(xml, read_only=True)
    assert jit_step(sim) is None

"""Replay of the MD steps as CUDA graphs (`--jit` with the torch backend)."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

import contextlib
import json
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest

import ipi.utils.array_backend as ab
from ipi_tests.unit_tests.engine.test_backend_equivalence import (
    HARMONIC,
    INPUTS,
    ISOTROPIC,
    NO_FORCES,
    mts_input,
    write_input,
)

torch = pytest.importorskip("torch")
cuda = pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA not available")

DUMMY = ("<ffdirect name='ff'><pes>dummy</pes></ffdirect>", "<force forcefield='ff'/>")
DUMMY_BATCHED = (
    "<ffdirect name='ff'><pes>dummy</pes><batch_size>4</batch_size></ffdirect>",
    "<force forcefield='ff'/>",
)
PROPERTIES = ["kinetic_md", "potential", "conserved", "temperature"]
# replays match eager steps to rounding errors, as the functions compiled for
# the replays (see ipi.utils.jit.compiled) may sum in a different order
TOLERANCE = dict(rtol=1e-12, atol=1e-14)


def simulation(xml, jit):
    from ipi.scripting import InteractiveSimulation

    return InteractiveSimulation(xml, backend="torch", device="cuda", jit=jit)


def graph_step(sim):
    """The object replaying the steps of the first system, if any."""

    from ipi.engine.jit import GraphStep

    step = getattr(getattr(sim.syslist[0].motion, "integrator", None), "step", None)
    return step if isinstance(step, GraphStep) else None


@cuda
@pytest.mark.parametrize(
    "case,forcefield,nsegments",
    [
        ("nve", NO_FORCES, 1),
        ("nvt_langevin", HARMONIC, 2),
        ("nvt_svr", HARMONIC, 2),
        ("pimd_nvt", DUMMY, 5),
        ("pimd_pile_g", HARMONIC, 2),
        ("pimd_nvt", DUMMY_BATCHED, 2),
        ("mts", None, 14),
        ("npt_isotropic", HARMONIC, 2),
        ("npt_flexible", HARMONIC, 2),
        ("pimd_npt", HARMONIC, 2),
        ("nst", None, 2),
    ],
)
def test_replay_matches_eager_steps(
    case, forcefield, nsegments, tmp_path, monkeypatch, device_random
):
    """The replayed steps give the same trajectory and outputs as eager steps
    (to rounding errors), with one CUDA graph per segment between eager force
    evaluations (of a force component, or only of the driver for direct
    forcefields)"""

    monkeypatch.chdir(tmp_path)
    xml = INPUTS[case](HARMONIC) if case in INPUTS else write_input(case, forcefield)
    try:
        eager, graph = simulation(xml, False), simulation(xml, True)
        for step in range(8):
            eager.run_step(step)
            graph.run_step(step)
            for name in PROPERTIES:
                np.testing.assert_allclose(
                    graph.properties(name), eager.properties(name), **TOLERANCE
                )
            a, b = eager.syslist[0], graph.syslist[0]
            torch.testing.assert_close(b.beads.q, a.beads.q, **TOLERANCE)
            torch.testing.assert_close(b.beads.p, a.beads.p, **TOLERANCE)
            torch.testing.assert_close(b.cell.h, a.cell.h, **TOLERANCE)
        assert len(graph_step(graph).segments) == nsegments
    finally:
        ab.set_backend("numpy")


@cuda
def test_momenta_written_between_replays(tmp_path, monkeypatch):
    """Momenta written through the normal modes between replayed steps (as by
    InteractiveSimulation.thermalize_momenta) are used by the next replays"""

    monkeypatch.chdir(tmp_path)
    xml = write_input("nve", HARMONIC)
    try:
        eager, graph = simulation(xml, False), simulation(xml, True)
        for sim in (eager, graph):
            for step in range(6):
                sim.run_step(step)
            sim.syslist[0].nm.pnm_[:] = 0.0
            for step in range(6, 9):
                sim.run_step(step)
        a, b = eager.syslist[0].beads, graph.syslist[0].beads
        torch.testing.assert_close(b.q, a.q, **TOLERANCE)
        torch.testing.assert_close(b.p, a.p, **TOLERANCE)
    finally:
        ab.set_backend("numpy")


@cuda
@pytest.mark.parametrize("step", [1, 2])
@pytest.mark.parametrize("case", ["pimd_nvt", "npt_isotropic"])
def test_state_written_during_the_eager_steps(
    case, step, tmp_path, monkeypatch, device_random
):
    """Momenta and positions written between the eager steps, which change the
    objects that the recorded step reads and where it is split, do not change
    the replayed steps"""

    monkeypatch.chdir(tmp_path)
    xml = write_input(case, HARMONIC)
    try:
        eager, graph = simulation(xml, False), simulation(xml, True)
        for sim in (eager, graph):
            for i in range(8):
                if i == step:
                    sim.thermalize_momenta()
                    sim.syslist[0].beads.q = sim.syslist[0].beads.q * 1.001
                sim.run_step(i)
        assert graph_step(graph).segments is not None
        a, b = eager.syslist[0], graph.syslist[0]
        torch.testing.assert_close(b.beads.q, a.beads.q, **TOLERANCE)
        torch.testing.assert_close(b.beads.p, a.beads.p, **TOLERANCE)
        torch.testing.assert_close(b.cell.h, a.cell.h, **TOLERANCE)
    finally:
        ab.set_backend("numpy")


@cuda
def test_forces_converted_in_part_before_the_recording(
    tmp_path, monkeypatch, device_random
):
    """The results of a forcefield computed before the recorded step, but only
    in part converted to device arrays (here the potential), are converted
    eagerly when the step is recorded"""

    monkeypatch.chdir(tmp_path)
    xml = write_input("nvt_langevin", HARMONIC)
    try:
        eager, graph = simulation(xml, False), simulation(xml, True)
        for sim in (eager, graph):
            for i in range(6):
                if i == 2:
                    sim.syslist[0].beads.q = sim.syslist[0].beads.q * 1.001
                    sim.properties("potential")
                sim.run_step(i)
        assert graph_step(graph).segments is not None
        a, b = eager.syslist[0].beads, graph.syslist[0].beads
        torch.testing.assert_close(b.q, a.q, **TOLERANCE)
        torch.testing.assert_close(b.p, a.p, **TOLERANCE)
    finally:
        ab.set_backend("numpy")


@cuda
def test_recording_in_a_new_thread(tmp_path, monkeypatch, device_random):
    """A step is recorded in a thread that has run eager steps, as some
    libraries (e.g. cuBLAS) initialize per-thread state that cannot be
    recorded: steps run in threads (threading='true') then still replay"""

    from concurrent.futures import ThreadPoolExecutor

    monkeypatch.chdir(tmp_path)
    xml = write_input("nvt_langevin", NO_FORCES)
    try:
        eager, graph = simulation(xml, False), simulation(xml, True)
        for step in range(2):
            eager.run_step(step)
            graph.run_step(step)
        with ThreadPoolExecutor(max_workers=1) as executor:
            for step in range(2, 6):
                eager.run_step(step)
                executor.submit(graph.run_step, step).result()
        assert graph_step(graph).segments is not None
        a, b = eager.syslist[0].beads, graph.syslist[0].beads
        torch.testing.assert_close(b.q, a.q, **TOLERANCE)
        torch.testing.assert_close(b.p, a.p, **TOLERANCE)
    finally:
        ab.set_backend("numpy")


@cuda
def test_unsupported_steps_are_eager(tmp_path, monkeypatch, device_random):
    monkeypatch.chdir(tmp_path)
    try:
        sim = simulation(write_input("nvt_gle", HARMONIC), True)
        assert graph_step(sim) is None
        sim.run_step(0)
    finally:
        ab.set_backend("numpy")


@cuda
@pytest.mark.parametrize(
    "case,owner,name",
    [("nvt_langevin", "motion", "dt"), ("npt_isotropic", "ensemble", "pext")],
)
def test_changed_parameters_stop_the_replay(
    case, owner, name, tmp_path, monkeypatch, device_random
):
    """The time step, and the pressure of the barostats, are fixed in the graphs"""

    from ipi.utils.jit import JitError

    monkeypatch.chdir(tmp_path)
    try:
        sim = simulation(write_input(case, HARMONIC), True)
        for step in range(4):
            sim.run_step(step)
        obj = getattr(sim.syslist[0], owner)
        setattr(obj, name, getattr(obj, name) * 1.1)
        with pytest.raises(JitError, match="cannot change"):
            sim.run_step(4)
    finally:
        ab.set_backend("numpy")


@cuda
@pytest.mark.parametrize("segment", [1, 2])
def test_failed_recording_stops_the_replay(
    segment, tmp_path, monkeypatch, device_random
):
    """An error raised while recording the first or the last segment (both
    apply the thermostat) stops the steps, also at the following calls"""

    from ipi.engine.thermostats import ThermoLangevin
    from ipi.utils import depend
    from ipi.utils.jit import JitError

    step, calls = ThermoLangevin.step, []

    def failing(self):
        if depend._recorder is not None:
            calls.append(None)
            if len(calls) == segment:
                raise RuntimeError("recording failed")
        step(self)

    monkeypatch.setattr(ThermoLangevin, "step", failing)
    monkeypatch.chdir(tmp_path)
    try:
        sim = simulation(write_input("nvt_langevin", HARMONIC), True)
        for i in range(2):
            sim.run_step(i)
        with pytest.raises(JitError, match="cannot be replayed"):
            sim.run_step(2)
        for i in range(3, 5):
            with pytest.raises(JitError, match="has failed"):
                sim.run_step(i)
    finally:
        ab.set_backend("numpy")


@cuda
def test_driver_extras_are_updated(tmp_path, monkeypatch, device_random):
    """The extras returned by the driver reach the forces at each replayed step"""

    from ipi.pes.dummy import Zero_driver

    calls = []

    def compute_structure(self, cell, pos):
        calls.append(None)
        return pos.sum() * 0.0, pos * 0.0, cell * 0.0, json.dumps({"n": len(calls)})

    monkeypatch.setattr(Zero_driver, "compute_structure", compute_structure)
    monkeypatch.chdir(tmp_path)
    try:
        sim = simulation(write_input("nve", DUMMY), True)
        for step in range(6):
            sim.run_step(step)
            assert sim.syslist[0].forces.extras["n"] == [len(calls)]
    finally:
        ab.set_backend("numpy")


@cuda
@pytest.mark.parametrize("case", ["nvt_langevin", "pimd_nvt"])
def test_temperature_changes_are_replayed(case, tmp_path, monkeypatch, device_random):
    """The coefficients of the Langevin thermostats are scalars on the device,
    which the replays read: they follow a change of the temperature"""

    monkeypatch.chdir(tmp_path)
    xml = write_input(case, HARMONIC)
    try:
        eager, graph = simulation(xml, False), simulation(xml, True)
        for step in range(8):
            if step == 4:
                eager.syslist[0].ensemble.temp *= 1.5
                graph.syslist[0].ensemble.temp *= 1.5
            eager.run_step(step)
            graph.run_step(step)
        assert graph_step(graph).segments is not None
        for name in PROPERTIES:
            np.testing.assert_allclose(
                graph.properties(name), eager.properties(name), **TOLERANCE
            )
        a, b = eager.syslist[0].beads, graph.syslist[0].beads
        torch.testing.assert_close(b.p, a.p, **TOLERANCE)
    finally:
        ab.set_backend("numpy")


@cuda
@pytest.mark.parametrize("forcefield", [NO_FORCES, DUMMY_BATCHED])
def test_replay_is_sync_free(forcefield, tmp_path, monkeypatch, device_random):
    """A replayed step does not wait for the device, also with the forces of a
    driver that does not wait for it"""

    monkeypatch.chdir(tmp_path)
    try:
        sim = simulation(write_input("pimd_nvt", forcefield), True)
        for step in range(4):
            sim.run_step(step)
        torch.cuda.set_sync_debug_mode("error")
        for step in range(4, 8):
            sim.run_step(step)
    finally:
        torch.cuda.set_sync_debug_mode("default")
        ab.set_backend("numpy")


# The tests below run without a GPU: the CUDA streams and graphs are replaced
# by objects that run the recorded step eagerly and replay nothing, so that a
# replay runs only the Python code around the graphs (the staging and the
# restoring of the depend objects, and the eager calls).


class FakeStream:
    def wait_stream(self, stream):
        pass


class FakeGraph:
    def __init__(self):
        self.capturing = False

    def register_generator_state(self, generator):
        pass

    def capture_begin(self, pool=None, capture_error_mode=None):
        self.capturing = True

    def capture_end(self):
        if not self.capturing:
            raise RuntimeError("not capturing")
        self.capturing = False

    def replay(self):
        pass


@pytest.fixture
def fake_cuda(monkeypatch):
    """Replays the steps with the torch backend on the CPU, with fake CUDA
    streams and graphs, and CPU tensors in place of CUDA tensors."""

    from ipi.utils import cudagraphs, jit

    monkeypatch.setattr(torch.cuda, "Stream", FakeStream)
    monkeypatch.setattr(torch.cuda, "current_stream", FakeStream)
    monkeypatch.setattr(torch.cuda, "set_stream", lambda stream: None)
    monkeypatch.setattr(torch.cuda, "stream", lambda stream: contextlib.nullcontext())
    monkeypatch.setattr(torch.cuda, "CUDAGraph", FakeGraph)
    monkeypatch.setattr(cudagraphs, "_is_cuda_tensor", torch.is_tensor)
    monkeypatch.setattr(jit, "_compiled", lambda func: func)
    ab.set_backend("torch", device="cpu")
    monkeypatch.setitem(ab._current, "jit", True)
    yield
    ab.set_backend("numpy")


def fake_simulation(xml):
    from ipi.scripting import InteractiveSimulation

    return InteractiveSimulation(xml)


@pytest.mark.parametrize(
    "case,forcefield,nsegments",
    [
        ("nve", NO_FORCES, 1),
        ("nvt_langevin", HARMONIC, 2),
        ("nvt_svr", HARMONIC, 2),
        ("pimd_nvt", DUMMY, 5),
        ("pimd_pile_g", HARMONIC, 2),
        ("pimd_nvt", DUMMY_BATCHED, 2),
        ("mts", None, 14),
    ],
)
def test_recording_is_split_at_force_evaluations(
    case, forcefield, nsegments, tmp_path, monkeypatch, fake_cuda
):
    """A step is recorded in segments that end at the eager force evaluations
    (of a force component, or only of the driver for direct forcefields)"""

    monkeypatch.chdir(tmp_path)
    xml = mts_input(HARMONIC) if case == "mts" else write_input(case, forcefield)
    sim = fake_simulation(xml)
    for step in range(4):
        sim.run_step(step)
    assert len(graph_step(sim).segments) == nsegments


def test_replays_call_the_driver(tmp_path, monkeypatch, fake_cuda):
    """At each replay the driver is called, its tensor results are copied into
    those read by the graphs, and its extras reach the forces"""

    from ipi.pes.dummy import Zero_driver

    calls = []

    def compute_structure(self, cell, pos):
        calls.append(None)
        pot = pos.sum() * 0.0 + len(calls)
        return pot, pos * 0.0, cell * 0.0, json.dumps({"n": len(calls)})

    monkeypatch.setattr(Zero_driver, "compute_structure", compute_structure)
    monkeypatch.chdir(tmp_path)
    sim = fake_simulation(write_input("nve", DUMMY))
    bead = sim.syslist[0].forces.mforces[0]._forces[0]
    for step in range(6):
        sim.run_step(step)
        assert sim.syslist[0].forces.extras["n"] == [len(calls)]
        assert float(bead.pot) == len(calls)
    assert graph_step(sim).segments is not None


def test_changed_parameters_stop_the_fake_replay(tmp_path, monkeypatch, fake_cuda):
    from ipi.utils.jit import JitError

    monkeypatch.chdir(tmp_path)
    sim = fake_simulation(write_input("nvt_langevin", HARMONIC))
    for step in range(4):
        sim.run_step(step)
    sim.syslist[0].motion.dt *= 1.1
    with pytest.raises(JitError, match="cannot change"):
        sim.run_step(4)


def test_replaced_generator_stops_the_replay(tmp_path, monkeypatch, fake_cuda):
    from ipi.utils.jit import JitError

    monkeypatch.chdir(tmp_path)
    sim = fake_simulation(write_input("nvt_langevin", NO_FORCES))
    for step in range(3):
        sim.run_step(step)
    sim.syslist[0].prng._torch_rng = torch.Generator()
    with pytest.raises(JitError, match="generator"):
        sim.run_step(3)


def test_function_is_recorded_after_eager_calls_in_the_thread(fake_cuda):
    """A function runs eagerly twice, and once in each new thread, before it is
    recorded and then replayed"""

    from ipi.utils.cudagraphs import GraphReplay

    calls = []
    replay = GraphReplay(lambda: calls.append(None))
    replay()
    replay()
    with ThreadPoolExecutor(max_workers=1) as executor:
        executor.submit(replay).result()
        assert replay.segments is None
        executor.submit(replay).result()
    assert len(replay.segments) == 1 and len(calls) == 4
    replay()
    assert len(calls) == 4


def test_inputs_of_a_segment(fake_cuda):
    """The depend objects that a segment reads before modifying them are the
    inputs that are brought up to date before its replays"""

    from ipi.utils.cudagraphs import GraphReplay
    from ipi.utils.depend import depend_value

    read = depend_value(name="read", value=torch.tensor(1.0))
    written = depend_value(name="written", value=torch.tensor(1.0))

    def func():
        written.set(read.get() * 2)
        written.get()

    replay = GraphReplay(func)
    for _ in range(3):
        replay()
    assert [dep._name for dep, _, _ in replay.segments[0].inputs] == ["read"]


def test_updated_objects_are_inputs_of_a_segment(fake_cuda):
    """The depend objects that a segment updates in place, or of which it
    assigns a part (also by writing a view), are inputs that are brought up to
    date before its replays"""

    from ipi.utils.cudagraphs import GraphReplay
    from ipi.utils.depend import depend_array

    added = depend_array(name="added", value=torch.zeros(3))
    assigned = depend_array(name="assigned", value=torch.zeros(3))
    replaced = depend_array(name="replaced", value=torch.zeros(3))
    viewed = depend_array(name="viewed", value=torch.zeros(3))
    view = viewed[1:]

    def func():
        added.__iadd__(1.0)
        assigned[0] = 1.0
        replaced[...] = 1.0
        view.set(torch.ones(2))

    replay = GraphReplay(func)
    for _ in range(3):
        replay()
    inputs = sorted(dep._name for dep, _, _ in replay.segments[0].inputs)
    assert inputs == ["added", "assigned", "viewed"]


@pytest.mark.parametrize("kind", ["array", "value"])
def test_host_writes_cannot_be_recorded(kind, fake_cuda):
    from ipi.utils.cudagraphs import GraphReplay
    from ipi.utils.jit import JitError
    from ipi.utils.depend import depend_array, depend_value

    if kind == "array":
        dep = depend_array(name="host", value=np.zeros(3), host=True)
        replay = GraphReplay(lambda: dep.__setitem__(0, 1.0))
    else:
        dep = depend_value(name="host", value=1.0)
        replay = GraphReplay(lambda: dep.set(2.0))
    replay()
    replay()
    with pytest.raises(JitError, match="host"):
        replay()
    assert replay.segments is None


@pytest.mark.parametrize("kind", ["capture", "eager"])
def test_failed_recording_is_final(kind, fake_cuda):
    """A function whose recording fails, during a capture or in an eager call
    between two segments, is not called again, as the kernels of the segment
    being recorded have not run: the following calls raise a JitError"""

    from ipi.utils import depend, jit
    from ipi.utils.cudagraphs import GraphReplay
    from ipi.utils.jit import JitError

    calls = []

    def fail():
        if len(calls) > jit.WARMUP_CALLS:
            raise RuntimeError("recording failed")
        return []

    def func():
        calls.append(None)
        if kind == "eager":
            jit.eager(fail)
        elif depend._recorder is not None:
            fail()

    replay = GraphReplay(func)
    replay()
    replay()
    with pytest.raises(Exception, match="recording failed|cannot be replayed"):
        replay()
    for _ in range(2):
        with pytest.raises(JitError, match="has failed"):
            replay()
    assert replay.segments is None and len(calls) == 3


@pytest.mark.parametrize("end_fails", [False, True])
@pytest.mark.parametrize("kind", ["other", "jit"])
def test_errors_during_a_capture(kind, end_fails, fake_cuda, monkeypatch):
    """An error raised while a segment is captured ends the capture, also when
    ending it fails, and is raised as a JitError; JitErrors are raised as they
    are"""

    from ipi.utils import depend
    from ipi.utils.cudagraphs import GraphReplay
    from ipi.utils.jit import JitError

    def failing_end(graph):
        raise RuntimeError("the capture cannot be ended")

    if end_fails:
        monkeypatch.setattr(FakeGraph, "capture_end", failing_end)
    error = JitError("raised in the step") if kind == "jit" else ValueError()

    def func():
        if depend._recorder is not None:
            raise error

    replay = GraphReplay(func)
    replay()
    replay()
    message = "raised in the step" if kind == "jit" else "cannot be replayed"
    with pytest.raises(JitError, match=message):
        replay()
    assert replay.segments is None


def test_errors_of_eager_calls_are_raised_as_they_are(fake_cuda):
    """An error raised by a function called eagerly between two recorded
    segments (e.g. by the driver of a forcefield) is raised as it is"""

    from ipi.utils import jit
    from ipi.utils.cudagraphs import GraphReplay

    calls = []

    def driver():
        calls.append(None)
        if len(calls) > jit.WARMUP_CALLS:
            raise ValueError("driver failed")
        return [torch.zeros(1)]

    replay = GraphReplay(lambda: jit.eager(driver))
    replay()
    replay()
    with pytest.raises(ValueError, match="driver failed"):
        replay()
    assert replay.segments is None


def test_staging_updates_the_inputs_of_a_graph(fake_cuda):
    """Before a replay, scalars replaced by new objects are copied into those
    read by the graph"""

    from ipi.utils.cudagraphs import _stage
    from ipi.utils.jit import JitError
    from ipi.utils.depend import depend_array, depend_value

    recorded = torch.tensor(1.0, dtype=torch.float64)
    scalar = depend_value(name="scalar", value=recorded)
    scalar.set(torch.tensor(2.0, dtype=torch.float64))
    _stage(scalar, recorded)
    scalar.set(3.0)
    _stage(scalar, recorded)
    assert float(recorded) == 3.0
    scalar.set(torch.zeros(2, dtype=torch.float64))
    with pytest.raises(JitError, match="changed shape"):
        _stage(scalar, recorded)

    array = depend_array(name="array", value=torch.zeros(3, dtype=torch.float64))
    _stage(array, array._data)
    with pytest.raises(JitError, match="has been replaced"):
        _stage(array, torch.zeros(3, dtype=torch.float64))


def test_values_fixed_in_a_graph_are_checked(fake_cuda):
    """Before a replay, the host values fixed in the graphs must be unchanged"""

    from ipi.utils.jit import JitError, _check
    from ipi.utils.depend import depend_array, depend_value

    host = depend_value(name="host", value=[1.0, 2.0])
    _check(host, [1.0, 2.0], None)
    for value in [[1.0, 3.0], [1.0, 2.0, 3.0]]:
        host.set(value)
        with pytest.raises(JitError, match="has changed"):
            _check(host, [1.0, 2.0], None)
    single = depend_value(name="single", value=[2.0, 2.0])
    with pytest.raises(JitError, match="has changed"):
        _check(single, [2.0], None)
    held = np.zeros(2)
    host.set(held)
    copy = held.copy()
    held[0] = 1.0
    with pytest.raises(JitError, match="has changed"):
        _check(host, held, copy)

    names = depend_array(name="names", value=np.zeros(3), host=True)
    copy = names._data.copy()
    names[0] = 1.0
    with pytest.raises(JitError, match="has changed"):
        _check(names, names._data, copy)


def test_fixed_values_are_checked_before_the_replay(fake_cuda):
    """A change of a value fixed in a later segment stops the replay before
    any segment runs, rather than after a part of the step"""

    from ipi.utils import jit
    from ipi.utils.cudagraphs import GraphReplay
    from ipi.utils.depend import depend_value
    from ipi.utils.jit import JitError

    calls = []
    dt = depend_value(name="dt", value=1.0)
    y = depend_value(name="y", value=torch.zeros(1))

    def func():
        jit.eager(lambda: calls.append(None) or [])
        y.set(y.get() + dt.get())

    replay = GraphReplay(func)
    for _ in range(4):
        replay()
    assert len(replay.segments) == 2
    dt.set(2.0)
    with pytest.raises(JitError, match="cannot change"):
        replay()
    assert len(calls) == 4


@pytest.mark.parametrize("capability,is_compiled", [((6, 1), False), ((9, 0), True)])
def test_functions_are_compiled_on_recent_gpus(capability, is_compiled, monkeypatch):
    from ipi.utils.jit import _compiled

    monkeypatch.setattr(torch.cuda, "get_device_capability", lambda d: capability)
    monkeypatch.setattr(torch, "compile", lambda func, **kwargs: lambda x: "compiled")

    def func(x):
        return x

    assert (_compiled(func)(1) == "compiled") == is_compiled


@pytest.mark.parametrize("stage", ["wrapping", "first call"])
def test_functions_that_fail_to_compile_are_called_as_they_are(stage, monkeypatch):
    """Compilation errors, raised by torch.compile when it wraps the function
    or at the first call of the compiled function, make the function run as it
    is from then on"""

    from ipi.utils.jit import _compiled

    attempts = []

    def failing(*args):
        attempts.append(None)
        raise RuntimeError("compilation failed")

    def compile(func, **kwargs):
        if stage == "wrapping":
            failing()
        return failing

    monkeypatch.setattr(torch.cuda, "get_device_capability", lambda d: (9, 0))
    monkeypatch.setattr(torch, "compile", compile)

    def func(x):
        return 2 * x

    call = _compiled(func)
    assert call(1) == 2 and call(2) == 4
    assert len(attempts) == 1


def test_compilation_errors_during_a_recording_are_raised(monkeypatch):
    """A function first compiled while a step is recorded as CUDA graphs, whose
    capture the error has invalidated, is not called as it is"""

    from ipi.utils import depend
    from ipi.utils.jit import _compiled

    def failing(*args):
        raise RuntimeError("compilation failed")

    monkeypatch.setattr(torch.cuda, "get_device_capability", lambda d: (9, 0))
    monkeypatch.setattr(torch, "compile", lambda func, **kwargs: failing)
    monkeypatch.setattr(depend, "_recorder", object())

    def func(x):
        return 2 * x

    with pytest.raises(RuntimeError, match="compilation failed"):
        _compiled(func)(1)


def test_only_cuda_tensors_are_device_values():
    from ipi.utils.cudagraphs import _is_cuda_tensor

    assert not _is_cuda_tensor(torch.zeros(1)) and not _is_cuda_tensor(1.0)


def test_values_that_cannot_be_compared_are_different():
    """Values that cannot be converted to arrays, or whose comparison raises,
    count as changed"""

    from ipi.utils.jit import _equal

    class Unconvertible:
        def __array__(self, dtype=None, copy=None):
            raise TypeError("cannot be converted")

    class Incomparable:
        def __eq__(self, other):
            raise TypeError("cannot be compared")

    assert _equal([1.0, 2.0], [1.0, 2.0]) and not _equal([1.0], [1.0, 1.0])
    assert not _equal(Unconvertible(), 1.0)
    assert not _equal(Incomparable(), Incomparable())


def test_eager_call_results_are_copied_into_the_recorded_ones(fake_cuda):
    from ipi.utils.cudagraphs import _Call
    from ipi.utils.jit import JitError

    results = iter(
        [
            [[torch.tensor(1.0), "first"]],
            [[torch.tensor(2.0), "second"]],
            [[torch.zeros(2), "third"]],
            [[torch.tensor(1.0)], []],
        ]
    )
    call = _Call(lambda: next(results), ())
    call()
    recorded = call.outputs
    scalar = recorded[0][0]
    call()
    assert recorded[0][0] is scalar and float(scalar) == 2.0
    assert recorded[0][1] == "second"
    with pytest.raises(JitError, match="changed shape"):
        call()
    with pytest.raises(JitError, match="changed structure"):
        call()


@pytest.mark.parametrize(
    "case,change,reason",
    [
        ("nve", None, None),
        ("nvt_svr", None, None),
        ("pimd_pile_g", None, None),
        ("nvt_gle", None, "the thermostat 'ThermoGLE'"),
        ("npt_flexible", None, None),
        (
            "npt_isotropic",
            (ISOTROPIC, ISOTROPIC.replace("'langevin'", "'svr'")),
            "the thermostat 'ThermoSVR' of the barostat",
        ),
        (
            "npt_isotropic",
            (
                ISOTROPIC,
                "<barostat mode='stochastic-rescaling'><tau units='femtosecond'>200"
                "</tau><compressibility units='bar^-1'>4.5e-5</compressibility>"
                "</barostat>",
            ),
            "the barostat 'BaroSCR'",
        ),
        ("nvt_langevin", ("mode='dynamics'", "mode='dummy'"), "the motion class"),
        ("pimd_nvt", ("mode='nvt'", "mode='nvt-cc'"), "the ensemble 'nvt-cc'"),
        (
            "nvt_langevin",
            (
                "<motion mode='dynamics'>",
                "<motion mode='dynamics'><fixatoms>[0]</fixatoms>",
            ),
            "fixed atoms",
        ),
        (
            "pimd_nvt",
            (
                "<motion mode='dynamics'>",
                "<normal_modes><open_paths>[0]</open_paths></normal_modes>"
                "<motion mode='dynamics'>",
            ),
            "open paths or bosons",
        ),
    ],
)
def test_unsupported_steps(case, change, reason, tmp_path, monkeypatch):
    from ipi.engine.jit import replay_steps, unsupported
    from ipi.scripting import InteractiveSimulation

    monkeypatch.chdir(tmp_path)
    xml = write_input(case, NO_FORCES)
    if change is not None:
        xml = xml.replace(*change)
    sim = InteractiveSimulation(xml)
    if reason is None:
        assert unsupported(sim, sim.syslist[0]) is None
    else:
        assert reason in unsupported(sim, sim.syslist[0])
        replay_steps(sim)
        assert graph_step(sim) is None


def two_systems(xml, threading, smotion=""):
    """Duplicates the system of an input, with an optional motion of the systems."""

    start, end = xml.index("<system>"), xml.index("</system>") + len("</system>")
    system = xml[start:end]
    systems = system.replace("<system>", "<system prefix='a'>") + system.replace(
        "<system>", "<system prefix='b'>"
    )
    xml = xml[:start] + systems + smotion + xml[end:]
    return xml.replace("threading='false'", f"threading='{threading}'")


@pytest.mark.parametrize(
    "threading,smotion,reason",
    [
        ("true", "", "stepping systems in threads"),
        (
            "false",
            "<smotion mode='remd'><remd><stride>4</stride></remd></smotion>",
            "a motion of multiple systems",
        ),
    ],
)
def test_unsupported_multiple_systems(
    threading, smotion, reason, tmp_path, monkeypatch
):
    from ipi.engine.jit import unsupported
    from ipi.scripting import InteractiveSimulation

    monkeypatch.chdir(tmp_path)
    sim = InteractiveSimulation(
        two_systems(write_input("nvt_langevin", NO_FORCES), threading, smotion)
    )
    assert unsupported(sim, sim.syslist[0]) == reason

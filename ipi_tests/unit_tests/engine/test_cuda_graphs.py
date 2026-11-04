"""Replay of the MD steps as CUDA graphs (`--jit` with the torch backend)."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

import contextlib
import json
from types import SimpleNamespace
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest

import ipi.utils.array_backend as ab
from ipi_tests.unit_tests.engine.inputs import (
    DUMMY,
    DUMMY_BATCHED,
    HARMONIC,
    NO_FORCES,
    PROPERTIES,
    TOLERANCE,
    case_input,
    write_input,
)

torch = pytest.importorskip("torch")
cuda = pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA not available")

# the options of the simulations replayed on the GPU
CUDA = dict(backend="torch", device="cuda")


def recording():
    """Whether a function is being recorded, rather than run (and observed)
    eagerly."""

    from ipi.utils import depend

    return depend._recorder is not None and not depend._recorder.observes


def graph_step(sim):
    """The object replaying the steps of the first system, if any."""

    from ipi.engine.jit import GraphStep

    step = getattr(getattr(sim.syslist[0].motion, "integrator", None), "step", None)
    return step if isinstance(step, GraphStep) else None


@pytest.fixture
def eager_driver(monkeypatch):
    """Makes the `dummy` driver not traceable, so that it is called eagerly
    between the CUDA graphs."""

    from ipi.pes.dummy import Zero_driver

    monkeypatch.setattr(Zero_driver, "traceable", False)


@cuda
@pytest.mark.parametrize(
    "case,forcefield,nsegments",
    [
        ("nve", NO_FORCES, 1),
        ("nvt_langevin", HARMONIC, 2),
        ("nvt_svr", HARMONIC, 2),
        ("pimd_nvt", DUMMY, 1),
        ("pimd_pile_g", HARMONIC, 2),
        ("pimd_nvt", DUMMY_BATCHED, 1),
        ("mts", HARMONIC, 2),
        ("npt_isotropic", HARMONIC, 2),
        ("npt_flexible", HARMONIC, 2),
        ("pimd_npt", HARMONIC, 2),
        ("nst", HARMONIC, 2),
    ],
)
def test_replay_matches_eager_steps(
    case, forcefield, nsegments, simulation, device_random
):
    """The replayed steps give the same trajectory and outputs as eager steps
    (to rounding errors), with one CUDA graph per segment between eager force
    evaluations (of a force component), and a single one with a traceable
    driver"""

    assert_replay_matches_eager_steps(
        case_input(case, forcefield), nsegments, simulation
    )


EAGER_DRIVER_CASES = [
    ("pimd_nvt", DUMMY, 5),
    ("pimd_nvt", DUMMY_BATCHED, 2),
    ("mts", HARMONIC, 14),
]


@cuda
@pytest.mark.parametrize("case,forcefield,nsegments", EAGER_DRIVER_CASES)
def test_replay_with_eager_driver_calls(
    case, forcefield, nsegments, simulation, device_random, eager_driver
):
    """A driver that is not traceable is called eagerly between the CUDA graphs,
    which end at each call"""

    xml = case_input(case, forcefield)
    assert_replay_matches_eager_steps(xml, nsegments, simulation)


def assert_replay_matches_eager_steps(xml, nsegments, simulation):
    with simulation(xml, **CUDA) as eager, simulation(xml, **CUDA, jit=True) as graph:
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


@cuda
def test_momenta_written_between_replays(simulation):
    """Momenta written through the normal modes between replayed steps (as by
    InteractiveSimulation.thermalize_momenta) are used by the next replays"""

    xml = write_input("nve", HARMONIC)
    with simulation(xml, **CUDA) as eager, simulation(xml, **CUDA, jit=True) as graph:
        for sim in (eager, graph):
            for step in range(6):
                sim.run_step(step)
            sim.syslist[0].nm.pnm_[:] = 0.0
            for step in range(6, 9):
                sim.run_step(step)
        a, b = eager.syslist[0].beads, graph.syslist[0].beads
        torch.testing.assert_close(b.q, a.q, **TOLERANCE)
        torch.testing.assert_close(b.p, a.p, **TOLERANCE)


@cuda
@pytest.mark.parametrize("step", [1, 2])
@pytest.mark.parametrize("case", ["pimd_nvt", "npt_isotropic"])
def test_state_written_during_the_eager_steps(case, step, simulation, device_random):
    """Momenta and positions written between the eager steps, which change the
    objects that the recorded step reads and where it is split, do not change
    the replayed steps"""

    xml = write_input(case, HARMONIC)
    with simulation(xml, **CUDA) as eager, simulation(xml, **CUDA, jit=True) as graph:
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


@cuda
def test_forces_converted_in_part_before_the_recording(simulation, device_random):
    """The results of a forcefield computed before the recorded step, but only
    in part converted to device arrays (here the potential), are converted
    eagerly when the step is recorded"""

    xml = write_input("nvt_langevin", HARMONIC)
    with simulation(xml, **CUDA) as eager, simulation(xml, **CUDA, jit=True) as graph:
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


@cuda
def test_recording_in_a_new_thread(simulation, device_random):
    """A step is recorded in a thread that has run eager steps, as some
    libraries (e.g. cuBLAS) initialize per-thread state that cannot be
    recorded: steps run in threads (threading='true') then still replay"""

    from concurrent.futures import ThreadPoolExecutor

    xml = write_input("nvt_langevin", NO_FORCES)
    with simulation(xml, **CUDA) as eager, simulation(xml, **CUDA, jit=True) as graph:
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


@cuda
def test_unsupported_steps_are_eager(simulation, device_random):
    with simulation(write_input("nvt_gle", HARMONIC), **CUDA, jit=True) as sim:
        assert graph_step(sim) is None
        sim.run_step(0)


@cuda
@pytest.mark.parametrize("segment", [1, 2])
def test_failed_recording_stops_the_replay(
    segment, simulation, monkeypatch, device_random
):
    """An error raised while recording the first or the last segment (both
    apply the thermostat) stops the steps, also at the following calls"""

    from ipi.engine.thermostats import ThermoLangevin
    from ipi.utils.jit import JitError

    step, calls = ThermoLangevin.step, []

    def failing(self):
        if recording():
            calls.append(None)
            if len(calls) == segment:
                raise RuntimeError("recording failed")
        step(self)

    monkeypatch.setattr(ThermoLangevin, "step", failing)
    with simulation(write_input("nvt_langevin", HARMONIC), **CUDA, jit=True) as sim:
        for i in range(2):
            sim.run_step(i)
        with pytest.raises(JitError, match="cannot be replayed"):
            sim.run_step(2)
        for i in range(3, 5):
            with pytest.raises(JitError, match="has failed"):
                sim.run_step(i)


@cuda
def test_driver_extras_are_updated(
    simulation, monkeypatch, device_random, eager_driver
):
    """The extras returned by an eager driver reach the forces at each replayed
    step"""

    from ipi.pes.dummy import Zero_driver

    calls = []

    def compute_structure(self, cell, pos):
        calls.append(None)
        return pos.sum() * 0.0, pos * 0.0, cell * 0.0, json.dumps({"n": len(calls)})

    monkeypatch.setattr(Zero_driver, "compute_structure", compute_structure)
    with simulation(write_input("nve", DUMMY), **CUDA, jit=True) as sim:
        for step in range(6):
            sim.run_step(step)
            assert sim.syslist[0].forces.extras["n"] == [len(calls)]


@cuda
@pytest.mark.parametrize("case", ["nvt_langevin", "pimd_nvt"])
def test_temperature_changes_are_replayed(case, simulation, device_random):
    """The coefficients of the Langevin thermostats are scalars on the device,
    which the replays read: they follow a change of the temperature"""

    xml = write_input(case, HARMONIC)
    with simulation(xml, **CUDA) as eager, simulation(xml, **CUDA, jit=True) as graph:
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


@cuda
@pytest.mark.parametrize("forcefield", [NO_FORCES, DUMMY_BATCHED])
def test_replay_is_sync_free(forcefield, simulation, device_random):
    """A replayed step does not wait for the device, also with the forces of a
    driver that does not wait for it"""

    xml = write_input("pimd_nvt", forcefield)
    with simulation(xml, **CUDA, jit=True) as sim:
        try:
            for step in range(4):
                sim.run_step(step)
            torch.cuda.set_sync_debug_mode("error")
            for step in range(4, 8):
                sim.run_step(step)
        finally:
            torch.cuda.set_sync_debug_mode("default")


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
def fake_cuda(simulation, select_backend, monkeypatch):
    """Replays the steps with the torch backend on the CPU, with fake CUDA
    streams and graphs, and CPU tensors in place of CUDA tensors: the
    simulations are created with `simulation(xml)`."""

    from ipi.utils import cudagraphs, jit

    monkeypatch.setattr(torch.cuda, "Stream", FakeStream)
    monkeypatch.setattr(torch.cuda, "current_stream", FakeStream)
    monkeypatch.setattr(torch.cuda, "set_stream", lambda stream: None)
    monkeypatch.setattr(torch.cuda, "stream", lambda stream: contextlib.nullcontext())
    monkeypatch.setattr(torch.cuda, "CUDAGraph", FakeGraph)
    monkeypatch.setattr(cudagraphs, "_is_cuda_tensor", torch.is_tensor)
    monkeypatch.setattr(jit, "_compiled", lambda func: func)
    select_backend("torch", device="cpu")
    monkeypatch.setitem(ab._current, "jit", True)


@pytest.mark.parametrize(
    "case,forcefield,nsegments",
    [
        ("nve", NO_FORCES, 1),
        ("nvt_langevin", HARMONIC, 2),
        ("nvt_svr", HARMONIC, 2),
        ("pimd_nvt", DUMMY, 1),
        ("pimd_pile_g", HARMONIC, 2),
        ("pimd_nvt", DUMMY_BATCHED, 1),
        ("mts", HARMONIC, 2),
    ],
)
def test_recording_is_split_at_force_evaluations(
    case, forcefield, nsegments, simulation, fake_cuda
):
    """A step is recorded in segments that end at the eager force evaluations
    (of a force component), and not at the calls of a traceable driver"""

    with simulation(case_input(case, forcefield)) as sim:
        for step in range(4):
            sim.run_step(step)
        assert len(graph_step(sim).segments) == nsegments


@pytest.mark.parametrize("case,forcefield,nsegments", EAGER_DRIVER_CASES)
def test_recording_is_split_at_eager_driver_calls(
    case, forcefield, nsegments, simulation, fake_cuda, eager_driver
):
    """The segments end at the calls of a driver that is not traceable"""

    with simulation(case_input(case, forcefield)) as sim:
        for step in range(4):
            sim.run_step(step)
        assert len(graph_step(sim).segments) == nsegments


def test_traceable_drivers_are_recorded(simulation, monkeypatch, fake_cuda):
    """A traceable driver is called by the eager steps and the recorded one,
    and not by the replays, which keep the extras of the recorded call"""

    from ipi.pes.dummy import Zero_driver

    calls = []

    def compute_structure(self, cell, pos):
        calls.append(None)
        return pos.sum() * 0.0, pos * 0.0, cell * 0.0, json.dumps({"n": len(calls)})

    monkeypatch.setattr(Zero_driver, "compute_structure", compute_structure)
    with simulation(write_input("nve", DUMMY)) as sim:
        for step in range(3):
            sim.run_step(step)
        recorded = len(calls)
        for step in range(3, 6):
            sim.run_step(step)
        assert len(graph_step(sim).segments) == 1
        assert len(calls) == recorded
        assert sim.syslist[0].forces.extras["n"] == [recorded]


def test_replays_call_the_driver(simulation, monkeypatch, fake_cuda, eager_driver):
    """At each replay an eager driver is called, its tensor results are copied
    into those read by the graphs, and its extras reach the forces"""

    from ipi.pes.dummy import Zero_driver

    calls = []

    def compute_structure(self, cell, pos):
        calls.append(None)
        pot = pos.sum() * 0.0 + len(calls)
        return pot, pos * 0.0, cell * 0.0, json.dumps({"n": len(calls)})

    monkeypatch.setattr(Zero_driver, "compute_structure", compute_structure)
    with simulation(write_input("nve", DUMMY)) as sim:
        bead = sim.syslist[0].forces.mforces[0]._forces[0]
        for step in range(6):
            sim.run_step(step)
            assert sim.syslist[0].forces.extras["n"] == [len(calls)]
            assert float(bead.pot) == len(calls)
        assert graph_step(sim).segments is not None


def test_changed_parameters_stop_the_fake_replay(simulation, fake_cuda):
    from ipi.utils.jit import JitError

    with simulation(write_input("nvt_langevin", HARMONIC)) as sim:
        for step in range(4):
            sim.run_step(step)
        sim.syslist[0].motion.dt *= 1.1
        with pytest.raises(JitError, match="cannot change"):
            sim.run_step(4)


def test_replaced_generator_stops_the_replay(simulation, fake_cuda):
    from ipi.utils.jit import JitError

    with simulation(write_input("nvt_langevin", NO_FORCES)) as sim:
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

    from ipi.utils import jit
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
        elif recording():
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

    from ipi.utils.cudagraphs import GraphReplay
    from ipi.utils.jit import JitError

    def failing_end(graph):
        raise RuntimeError("the capture cannot be ended")

    if end_fails:
        monkeypatch.setattr(FakeGraph, "capture_end", failing_end)
    error = JitError("raised in the step") if kind == "jit" else ValueError()

    def func():
        if recording():
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
    data = held.tobytes()
    held[0] = 1.0
    with pytest.raises(JitError, match="has changed"):
        _check(host, held, data)
    # a replaced host array must have the same shape and type
    host.set(np.zeros(2, dtype=np.int64))
    with pytest.raises(JitError, match="has changed"):
        _check(host, held, np.zeros(2).tobytes())

    names = depend_array(name="names", value=np.zeros(3), host=True)
    data = names._data.tobytes()
    names[0] = 1.0
    with pytest.raises(JitError, match="has changed"):
        _check(names, names._data, data)


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
    monkeypatch.setattr(depend, "_recorder", SimpleNamespace(observes=False))

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
    recorded = call.results
    scalar = recorded[0][0]
    call()
    assert recorded[0][0] is scalar and float(scalar) == 2.0
    assert recorded[0][1] == "second"
    with pytest.raises(JitError, match="changed shape"):
        call()
    with pytest.raises(JitError, match="changed structure"):
        call()

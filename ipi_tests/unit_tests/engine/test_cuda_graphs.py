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
    ISOTROPIC,
    NO_FORCES,
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


def mts_input():
    """PIMD with two forcefields applied at two levels of multiple time stepping."""

    xml = write_input("pimd_nvt", HARMONIC)
    return (
        xml.replace(
            HARMONIC[0],
            HARMONIC[0] + "<ffdirect name='ff2'><pes>dummy</pes></ffdirect>",
        )
        .replace(
            "<force forcefield='ff'/>",
            "<force forcefield='ff'><mts_weights>[1,0]</mts_weights></force>"
            "<force forcefield='ff2'><mts_weights>[0,1]</mts_weights></force>",
        )
        .replace(
            "<timestep units='femtosecond'>0.5</timestep>",
            "<timestep units='femtosecond'>0.5</timestep><nmts>[1,3]</nmts>",
        )
    )


def nst_input():
    """Classical dynamics at constant stress, with the anisotropic barostat."""

    box = (8 / 0.05) ** (1.0 / 3.0)  # the box of write_input, in angstrom
    return (
        write_input("npt_isotropic", HARMONIC)
        .replace(
            "<pressure units='bar'>1</pressure>",
            "<stress units='bar'>[1, 0, 0, 0, 1, 0, 0, 0, 1]</stress>",
        )
        .replace("mode='npt'", "mode='nst'")
        .replace(
            "<barostat mode='isotropic'>",
            "<barostat mode='anisotropic'><h0 units='angstrom'>"
            f"[{box}, 0, 0, 0, {box}, 0, 0, 0, {box}]</h0>",
        )
    )


INPUTS = {"mts": mts_input, "nst": nst_input}


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
    xml = INPUTS[case]() if case in INPUTS else write_input(case, forcefield)
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
    monkeypatch.setattr(torch.cuda, "graph_pool_handle", lambda: None)
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
    xml = mts_input() if case == "mts" else write_input(case, forcefield)
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
    with pytest.raises(JitError, match="cannot be replayed") as info:
        replay()
    assert "host" in str(info.value.__cause__)


def test_staging_updates_or_checks_the_inputs_of_a_graph(fake_cuda):
    """Before a replay, scalars replaced by new objects are copied into those
    read by the graph, and the values fixed in the graph must be unchanged"""

    from ipi.utils.cudagraphs import _stage
    from ipi.utils.jit import JitError
    from ipi.utils.depend import depend_array, depend_value

    recorded = torch.tensor(1.0, dtype=torch.float64)
    scalar = depend_value(name="scalar", value=recorded)
    scalar.set(torch.tensor(2.0, dtype=torch.float64))
    _stage(scalar, recorded, None)
    scalar.set(3.0)
    _stage(scalar, recorded, None)
    assert float(recorded) == 3.0
    scalar.set(torch.zeros(2, dtype=torch.float64))
    with pytest.raises(JitError, match="changed shape"):
        _stage(scalar, recorded, None)

    host = depend_value(name="host", value=[1.0, 2.0])
    _stage(host, [1.0, 2.0], None)
    for value in [[1.0, 3.0], [1.0, 2.0, 3.0]]:
        host.set(value)
        with pytest.raises(JitError, match="has changed"):
            _stage(host, [1.0, 2.0], None)
    held = np.zeros(2)
    host.set(held)
    copy = held.copy()
    held[0] = 1.0
    with pytest.raises(JitError, match="has changed"):
        _stage(host, held, copy)

    array = depend_array(name="array", value=torch.zeros(3, dtype=torch.float64))
    _stage(array, array._wvalue, None)
    with pytest.raises(JitError, match="has been replaced"):
        _stage(array, torch.zeros(3, dtype=torch.float64), None)
    names = depend_array(name="names", value=np.zeros(3), host=True)
    copy = names._wvalue.copy()
    names[0] = 1.0
    with pytest.raises(JitError, match="has changed"):
        _stage(names, names._wvalue, copy)


@pytest.mark.parametrize("capability,is_compiled", [((6, 1), False), ((9, 0), True)])
def test_functions_are_compiled_on_recent_gpus(capability, is_compiled, monkeypatch):
    from ipi.utils.jit import _compiled

    monkeypatch.setattr(torch.cuda, "get_device_capability", lambda d: capability)
    monkeypatch.setattr(torch, "compile", lambda func, **kwargs: "compiled")

    def func(x):
        return x

    assert (_compiled(func) == "compiled") == is_compiled


def test_only_cuda_tensors_are_device_values():
    from ipi.utils.cudagraphs import _is_cuda_tensor

    assert not _is_cuda_tensor(torch.zeros(1)) and not _is_cuda_tensor(1.0)


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

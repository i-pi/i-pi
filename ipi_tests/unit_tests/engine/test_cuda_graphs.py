"""Replay of the MD steps as CUDA graphs (`--cuda-graphs`)."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

import json

import numpy as np
import pytest

import ipi.utils.array_backend as ab
from ipi_tests.unit_tests.engine.test_backend_equivalence import (
    HARMONIC,
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


def simulation(xml, cuda_graphs):
    from ipi.scripting import InteractiveSimulation

    return InteractiveSimulation(
        xml, backend="torch", device="cuda", cuda_graphs=cuda_graphs
    )


def graph_step(sim):
    """The object replaying the steps of the first system, if any."""

    from ipi.engine.cudagraphs import GraphStep

    step = sim.syslist[0].motion.integrator.step
    return step if isinstance(step, GraphStep) else None


def test_cuda_graphs_require_torch_on_cuda():
    with pytest.raises(ValueError):
        ab.set_backend("numpy", cuda_graphs=True)
    with pytest.raises(ValueError):
        ab.set_backend("torch", device="cpu", cuda_graphs=True)
    ab.set_backend("numpy")


@cuda
@pytest.mark.parametrize(
    "case,forcefield,nsegments",
    [
        ("nve", NO_FORCES, 1),
        ("nvt_langevin", HARMONIC, 2),
        ("pimd_nvt", DUMMY, 5),
        ("pimd_nvt", DUMMY_BATCHED, 2),
        ("mts", None, 14),
    ],
)
def test_replay_matches_eager_steps(
    case, forcefield, nsegments, tmp_path, monkeypatch, device_random
):
    """The replayed steps give the same trajectory and outputs as eager steps,
    with one CUDA graph per segment between eager force evaluations (of a
    force component, or only of the driver for direct forcefields)"""

    monkeypatch.chdir(tmp_path)
    xml = mts_input() if case == "mts" else write_input(case, forcefield)
    try:
        eager, graph = simulation(xml, False), simulation(xml, True)
        for step in range(8):
            eager.run_step(step)
            graph.run_step(step)
            for name in PROPERTIES:
                np.testing.assert_array_equal(
                    eager.properties(name), graph.properties(name)
                )
            a, b = eager.syslist[0].beads, graph.syslist[0].beads
            assert torch.equal(a.q, b.q) and torch.equal(a.p, b.p)
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
        assert torch.equal(a.q, b.q) and torch.equal(a.p, b.p)
    finally:
        ab.set_backend("numpy")


@cuda
def test_unsupported_steps_are_eager(tmp_path, monkeypatch, device_random):
    monkeypatch.chdir(tmp_path)
    try:
        sim = simulation(write_input("nvt_svr", HARMONIC), True)
        assert graph_step(sim) is None
        sim.run_step(0)
    finally:
        ab.set_backend("numpy")


@cuda
def test_changed_parameters_stop_the_replay(tmp_path, monkeypatch, device_random):
    from ipi.utils.cudagraphs import CUDAGraphError

    monkeypatch.chdir(tmp_path)
    try:
        sim = simulation(write_input("nvt_langevin", HARMONIC), True)
        for step in range(4):
            sim.run_step(step)
        sim.syslist[0].ensemble.temp *= 1.1
        with pytest.raises(CUDAGraphError, match="cannot change"):
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

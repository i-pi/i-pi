"""Tests that a soft exit leaves outputs and RESTART in a consistent state.

The contract is that restarting from RESTART continues the outputs with no
gaps, duplicates or partial frames, whenever the exit request arrives. The
request is injected at a chosen point of the run, so the tests do not depend
on timing.
"""

import io
import os
import signal
import subprocess
import sys
import threading
import time

import numpy as np
import pytest

import ipi
import ipi.utils.softexit as softexit_module
from ipi.engine.simulation import Simulation
from ipi.utils.softexit import softexit

NSTEPS = 6
EXIT_STEP = 3  # index of the step during which the exit is requested


def simulation_xml(threading_mode):
    """Returns a small self-contained simulation with two outputs."""

    return f"""<simulation verbosity="quiet" threading="{threading_mode}">
  <output prefix="sim">
    <properties stride="1" filename="out">[ step, potential ]</properties>
    <trajectory stride="1" filename="pos">positions</trajectory>
  </output>
  <total_steps>{NSTEPS}</total_steps>
  <prng><seed>31415</seed></prng>
  <ffdirect name="harmonic">
    <pes>harmonic</pes>
    <parameters>{{k1: 0.1}}</parameters>
  </ffdirect>
  <system>
    <forces><force forcefield="harmonic"/></forces>
    <motion mode="dynamics">
      <fixcom>False</fixcom>
      <dynamics mode="nve">
        <timestep units="femtosecond">0.5</timestep>
      </dynamics>
    </motion>
    <beads natoms="1" nbeads="1">
      <q shape="(1,3)">[0.1,0,0]</q>
      <p shape="(1,3)">[1,0.2,-0.1]</p>
      <m shape="(1)">[1837.36223469]</m>
      <names shape="(1)">[H]</names>
    </beads>
    <cell shape="(3,3)">
      [18.897261,0,0,0,18.897261,0,0,0,18.897261]
    </cell>
  </system>
</simulation>"""


@pytest.fixture
def workdir(tmp_path, monkeypatch):
    """Runs in a scratch directory, with a responsive soft-exit monitor."""

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(softexit_module, "SOFTEXITLATENCY", 0.01)
    softexit.reset()
    yield tmp_path
    softexit.reset()


def exit_from_another_thread():
    """Requests a soft exit as the monitoring thread would."""

    thread = threading.Thread(target=softexit.cleanup)
    thread.start()
    thread.join(timeout=10)
    assert not thread.is_alive()


def exit_from_signal():
    """Delivers a kill signal to the main thread and waits for the request."""

    signal.raise_signal(signal.SIGTERM)
    deadline = time.time() + 10
    while not softexit.triggered and time.time() < deadline:
        time.sleep(0.001)
    assert softexit.triggered


def once_at_exit_step(simulation, action):
    """Wraps an action so it fires a single time, during the chosen step."""

    fired = []

    def guarded():
        if simulation.step == EXIT_STEP and not fired:
            fired.append(True)
            action()

    return guarded


def inject_before(obj, name, action):
    """Calls action right before obj.name()."""

    original = getattr(obj, name)

    def wrapped(*args, **kwargs):
        action()
        return original(*args, **kwargs)

    setattr(obj, name, wrapped)


def inject_after(obj, name, action):
    """Calls action right after obj.name()."""

    original = getattr(obj, name)

    def wrapped(*args, **kwargs):
        result = original(*args, **kwargs)
        action()
        return result

    setattr(obj, name, wrapped)


class StreamProxy:
    """A stream that calls action in the middle of writing a row."""

    def __init__(self, stream, action):
        self._stream = stream
        self._action = action
        self._nwrites = 0

    def write(self, data):
        self._nwrites += 1
        if self._nwrites % 2 == 0:
            self._action()
        return self._stream.write(data)

    def __getattr__(self, name):
        return getattr(self._stream, name)


def finish(simulation):
    """Runs the simulation and the end-of-run soft exit, without sys.exit."""

    try:
        simulation.run()
    except SystemExit:
        pytest.fail("The soft exit unwound the main thread in the middle of a run")
    softexit.cleanup()
    softexit.reset()


def restart_step():
    """Returns the step RESTART would continue from."""

    with open("RESTART") as restart_file:
        return Simulation.load_from_xml(restart_file, read_only=True).step


def assert_continuous():
    """Restarts from RESTART and checks the outputs have every step once."""

    with open("RESTART") as restart_file:
        restarted = Simulation.load_from_xml(restart_file)
    assert restarted.step < NSTEPS
    finish(restarted)

    steps = np.loadtxt("sim.out")[:, 0]
    np.testing.assert_array_equal(steps, np.arange(NSTEPS + 1))

    with open("sim.pos_0.xyz") as trajectory:
        lines = trajectory.readlines()
    # one atom: each frame has a count, a comment and a coordinate line
    assert len(lines) == 3 * (NSTEPS + 1)
    frames = [int(line.split("Step:")[1].split()[0]) for line in lines[1::3]]
    assert frames == list(range(NSTEPS + 1))


@pytest.mark.parametrize("threading_mode", ["false", "true"])
def test_exit_during_step_rolls_back(workdir, threading_mode):
    """A step interrupted halfway is discarded and redone on restart."""

    simulation = Simulation.load_from_xml(io.StringIO(simulation_xml(threading_mode)))
    inject_before(
        simulation.syslist[0].motion,
        "step",
        once_at_exit_step(simulation, exit_from_another_thread),
    )
    finish(simulation)

    assert restart_step() == EXIT_STEP
    assert_continuous()


@pytest.mark.parametrize("threading_mode", ["false", "true"])
def test_exit_between_outputs_completes_the_step(workdir, threading_mode):
    """Once a step started writing outputs, all of them are written."""

    simulation = Simulation.load_from_xml(io.StringIO(simulation_xml(threading_mode)))
    inject_after(
        simulation.outputs[0],
        "write",
        once_at_exit_step(simulation, exit_from_another_thread),
    )
    finish(simulation)

    assert restart_step() == EXIT_STEP + 1
    assert_continuous()


@pytest.mark.parametrize("threading_mode", ["false", "true"])
def test_exit_during_write_completes_the_row(workdir, threading_mode):
    """An output being written is not closed under the writer's feet."""

    simulation = Simulation.load_from_xml(io.StringIO(simulation_xml(threading_mode)))
    properties = simulation.outputs[0]
    properties.out = StreamProxy(
        properties.out, once_at_exit_step(simulation, exit_from_another_thread)
    )
    finish(simulation)

    assert_continuous()


def test_signal_during_outputs_completes_the_step(workdir):
    """A kill signal does not run the cleanup on top of an output write."""

    simulation = Simulation.load_from_xml(io.StringIO(simulation_xml("false")))
    inject_after(
        simulation.outputs[0],
        "write",
        once_at_exit_step(simulation, exit_from_signal),
    )
    finish(simulation)

    assert_continuous()


def test_pending_exit_runs_callbacks_once(workdir):
    """A request made while the lock is held is served later, a single time."""

    calls = []
    softexit.register_function(lambda: calls.append(True))

    with softexit.lock:
        exit_from_another_thread()
        assert softexit.triggered
        assert calls == []

    softexit.cleanup()
    softexit.cleanup()
    assert calls == [True]


def run_script(script):
    """Runs a script that uses the soft exit in a separate interpreter."""

    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        [os.path.dirname(os.path.dirname(ipi.__file__)), env.get("PYTHONPATH", "")]
    )
    header = """
import atexit, threading, time
from ipi.utils.softexit import softexit

def slow():
    time.sleep(0.5)
    open("done", "w").close()

softexit.register_function(slow)
atexit.register(lambda: open("finalized", "w").close())
"""
    return subprocess.run([sys.executable, "-c", header + script], env=env, timeout=60)


@pytest.mark.parametrize("hard_exit", [False, True])
def test_process_waits_for_cleanup(workdir, hard_exit):
    """The interpreter does not terminate in the middle of a soft exit."""

    process = run_script(f"""
softexit.hard_exit = {hard_exit}
threading.Thread(target=softexit.trigger, daemon=True).start()
deadline = time.time() + 10
while not softexit.exiting and time.time() < deadline:
    time.sleep(0.001)
# the end of bin/i-pi
softexit.trigger()
""")

    assert process.returncode == 0
    assert os.path.exists("done")


def test_hard_exit_skips_finalization(workdir):
    """A hard exit ends the process right after the cleanup."""

    process = run_script("""
softexit.hard_exit = True
softexit.trigger()
""")

    assert process.returncode == 0
    assert os.path.exists("done")
    assert not os.path.exists("finalized")


def test_hard_exit_from_another_thread_stops_the_run(workdir):
    """The main thread does not carry on after another thread cleaned up."""

    start = time.time()
    process = run_script("""
softexit.hard_exit = True
threading.Thread(target=softexit.trigger, daemon=True).start()
# a step that would take a long time to notice the exit
time.sleep(30)
open("carried_on", "w").close()
""")

    assert process.returncode == 0
    assert os.path.exists("done")
    assert not os.path.exists("carried_on")
    assert time.time() - start < 20

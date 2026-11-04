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
import time
import warnings

import numpy as np
import pytest

import ipi
import ipi.utils.softexit as softexit_module
from ipi.engine.outputs import BaseOutput
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
    """Runs in a scratch directory, checking often for the EXIT file."""

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(softexit_module, "SOFTEXITLATENCY", 0.01)
    softexit.reset()
    yield tmp_path
    softexit.reset()


def request_exit():
    """Requests a soft exit, as the EXIT file or the watchdog would."""

    softexit.request(" @SOFTEXIT: test")


def make_exit():
    """Makes a soft exit now, as when it is made while waiting for forces."""

    softexit.trigger(message=" @SOFTEXIT: test")


def exit_from_signal():
    """Delivers a kill signal to the main thread, which records the request."""

    signal.raise_signal(signal.SIGTERM)
    assert softexit.requested is not None


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
    assert_outputs_continuous()


def assert_outputs_continuous():
    """Checks the outputs of a completed run have every step once."""

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
        once_at_exit_step(simulation, request_exit),
    )
    finish(simulation)

    assert restart_step() == EXIT_STEP
    assert_continuous()


@pytest.mark.parametrize("threading_mode", ["false", "true"])
def test_request_during_outputs_completes_the_step(workdir, threading_mode):
    """A request made while the outputs are written, which need no more
    forces, is made at the start of the next step."""

    simulation = Simulation.load_from_xml(io.StringIO(simulation_xml(threading_mode)))
    inject_after(
        simulation.outputs[0],
        "write",
        once_at_exit_step(simulation, request_exit),
    )
    finish(simulation)

    assert restart_step() == EXIT_STEP + 1
    assert_continuous()


@pytest.mark.parametrize("threading_mode", ["false", "true"])
def test_exit_between_outputs_rolls_back(workdir, threading_mode):
    """An exit made after some of the outputs of a step are written removes
    them, and the step is redone on restart."""

    simulation = Simulation.load_from_xml(io.StringIO(simulation_xml(threading_mode)))
    inject_after(
        simulation.outputs[0],
        "write",
        once_at_exit_step(simulation, make_exit),
    )
    finish(simulation)

    assert restart_step() == EXIT_STEP
    assert_continuous()


@pytest.mark.parametrize("threading_mode", ["false", "true"])
def test_exit_during_a_row_rolls_back(workdir, threading_mode):
    """An exit made in the middle of a row removes it."""

    simulation = Simulation.load_from_xml(io.StringIO(simulation_xml(threading_mode)))
    properties = simulation.outputs[0]
    properties.out = StreamProxy(
        properties.out, once_at_exit_step(simulation, make_exit)
    )
    finish(simulation)

    assert restart_step() == EXIT_STEP
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


def test_request_runs_callbacks_once(workdir):
    """A request is made at the next check, a single time."""

    calls = []
    softexit.register_function(lambda: calls.append(True))

    request_exit()
    assert calls == []
    softexit.check()
    softexit.check()
    softexit.cleanup()
    assert calls == [True]


def test_outputs_go_back_to_the_stored_state(workdir):
    """An exit that goes back to the stored state truncates the outputs
    written since then, and removes those that were empty."""

    kept, new = BaseOutput("kept"), BaseOutput("new")
    for output in (kept, new):
        output.bind()
    kept.write("before\n")
    softexit.stored()
    kept.write("after\n")
    new.write("after\n")
    softexit.rollback = True
    softexit.cleanup()

    with open("kept") as stream:
        assert stream.read() == "before\n"
    assert not os.path.exists("new")


def test_forked_child_is_detached(workdir):
    """The signals of a forked child (e.g. a driver) end it, and do not reach
    the watchdog of the parent through the pipe they share."""

    softexit.start()
    with warnings.catch_warnings():
        # forking a process with threads
        warnings.simplefilter("ignore", DeprecationWarning)
        pid = os.fork()
    if pid == 0:
        os.kill(os.getpid(), signal.SIGTERM)
        os._exit(3)
    _, status = os.waitpid(pid, 0)
    assert os.WIFSIGNALED(status) and os.WTERMSIG(status) == signal.SIGTERM
    # the watchdog reads its pipe every SOFTEXITLATENCY
    time.sleep(0.2)
    assert softexit.requested is None


def test_second_signal_forces_the_exit(workdir):
    """The first Ctrl-C asks for a soft exit, the second one interrupts."""

    # as when several simulations are created in the same process
    softexit.start()
    softexit.start()
    signal.raise_signal(signal.SIGINT)
    with pytest.raises(KeyboardInterrupt):
        signal.raise_signal(signal.SIGINT)


def test_output_written_after_the_exit_is_discarded(workdir):
    """Outputs written in the middle of a step do not hit a closed file."""

    output = BaseOutput("floating")
    output.bind()
    output.write("before\n")
    softexit.cleanup()
    output.write("after\n")

    with open("floating") as stream:
        assert stream.read() == "before\n"


def run_script(script):
    """Runs a script that uses the soft exit in a separate interpreter."""

    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        [os.path.dirname(os.path.dirname(ipi.__file__)), env.get("PYTHONPATH", "")]
    )
    header = """
import atexit, os, signal, time
import ipi.utils.softexit as softexit_module
from ipi.utils.softexit import softexit

def slow():
    time.sleep(0.5)
    open("done", "w").close()

softexit.register_function(slow)
atexit.register(lambda: open("finalized", "w").close())
"""
    return subprocess.run([sys.executable, "-c", header + script], env=env, timeout=60)


@pytest.mark.parametrize("hard_exit", [False, True])
def test_pending_request_is_made_at_exit(workdir, hard_exit):
    """A request that the main thread did not get to is made when the
    interpreter terminates."""

    process = run_script(f"""
softexit.hard_exit = {hard_exit}
softexit.request("test")
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


WATCHDOG = """
softexit_module.GRACE = 0.5
softexit_module.SOFTEXITLATENCY = 0.05
softexit.hard_exit = True
softexit.register_emergency(lambda: open("emergency", "w").close())
softexit.start()
"""


@pytest.mark.parametrize(
    "exit_request", ["os.kill(os.getpid(), signal.SIGTERM)", 'open("EXIT", "w")']
)
def test_watchdog_ends_a_hung_process(workdir, exit_request):
    """A main thread that hangs in native code, where it runs neither the
    signal handlers nor the checks, is ended by the watchdog."""

    start = time.time()
    process = run_script(f"""{WATCHDOG}
from ipi.engine.outputs import BaseOutput
output = BaseOutput("floating")
output.bind()
output.write("before\\n")
softexit.stored()
output.write("after\\n")
output.out.flush()
# the signal reaches the watchdog thread, as when the main thread is stuck
signal.pthread_sigmask(signal.SIG_BLOCK, [signal.SIGTERM])
{exit_request}
time.sleep(30)
open("carried_on", "w").close()
""")

    assert process.returncode == 1
    assert os.path.exists("emergency")
    assert not os.path.exists("done")
    assert not os.path.exists("carried_on")
    assert time.time() - start < 20
    # the output goes back to the stored state with the checkpoint
    with open("floating") as stream:
        assert stream.read() == "before\n"


def test_watchdog_leaves_a_process_that_cleans_up(workdir):
    """A cleanup that takes longer than GRACE is not cut short."""

    process = run_script(f"""{WATCHDOG}
softexit_module.GRACE = 0.2
softexit.request("test")
softexit.check()
softexit.kill()
""")

    assert process.returncode == 0
    assert os.path.exists("done")
    assert not os.path.exists("emergency")


def test_status_and_message_are_keywords():
    """A message given in place of the status is refused before anything is done."""

    softexit.reset()
    try:
        with pytest.raises(TypeError):
            softexit.trigger("@FIRE is stuck")
        assert not softexit.triggered
    finally:
        softexit.reset()

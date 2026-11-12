"""Review tests of the output rollback made by a soft exit (PR #617).

Every test compares an interrupted-and-restarted run with an uninterrupted
reference run of the same input, made in another directory. The criterion
(C) is that, once the run restarted from RESTART has finished, every output
has each expected step exactly once, in order, with no partial row, and is
byte-identical to that of the reference (the forces come from an in-process
`ffdirect`, and the random numbers from a seeded generator stored in RESTART).
"""

import io
import os
import re
import time

import numpy as np
import pytest

import ipi.utils.softexit as softexit_module
from ipi.engine.simulation import Simulation
from ipi.utils.softexit import softexit

from .test_softexit import (
    finish,
    inject_after,
    inject_before,
    make_exit,
    request_exit,
)

NATOMS = 2
PROP_STRIDE = 1
POS_STRIDE = 2
VEL_STRIDE = 3
# round-off of a restart, many orders below the effect of a step that is
# missing, repeated, or made with other random numbers
RTOL, ATOL = 1e-9, 1e-12
STRIDES = {".out": PROP_STRIDE, ".pos_": POS_STRIDE, ".vel_": VEL_STRIDE}

OUTPUTS = f"""<properties stride="{PROP_STRIDE}" filename="out">
      [ step, conserved, potential, kinetic_md ]</properties>
    <trajectory stride="{POS_STRIDE}" filename="pos">positions</trajectory>
    <trajectory stride="{VEL_STRIDE}" filename="vel">velocities</trajectory>"""

CHECKPOINT = '<checkpoint stride="5" filename="chk" overwrite="true"/>'

REMD = "<smotion mode='remd'><remd><stride>1</stride></remd></smotion>"


def simulation_xml(
    nsteps,
    safe_stride=1,
    interleave="false",
    nbeads=4,
    temperatures=(300,),
    prefixes=("",),
    outputs=OUTPUTS,
    smotion="",
    total_time=None,
):
    """A path integral Langevin simulation of harmonic two-atom systems,
    extending the input of test_softexit.py to several beads and systems."""

    systems = []
    for i, (prefix, temperature) in enumerate(zip(prefixes, temperatures)):
        q = np.tile([0.2, -0.3, 0.1, -0.4, 0.2, 0.3], (nbeads, 1)) * (1 + 0.1 * i)
        q[:, 0] += np.arange(nbeads) * 0.01
        p = np.tile([0.1, 0.2, -0.1, 0.3, -0.2, 0.1], (nbeads, 1))
        systems.append(f"""<system prefix="{prefix}">
    <forces><force forcefield="harmonic"/></forces>
    <ensemble><temperature units="kelvin">{temperature}</temperature></ensemble>
    <motion mode="dynamics">
      <fixcom>False</fixcom>
      <dynamics mode="nvt">
        <timestep units="femtosecond">0.5</timestep>
        <thermostat mode="langevin"><tau units="femtosecond">10</tau></thermostat>
      </dynamics>
    </motion>
    <beads natoms="{NATOMS}" nbeads="{nbeads}">
      <q shape="({nbeads},6)">{q.ravel().tolist()}</q>
      <p shape="({nbeads},6)">{p.ravel().tolist()}</p>
      <m shape="(2)">[1837.3, 3674.6]</m>
      <names shape="(2)">[H, D]</names>
    </beads>
    <cell shape="(3,3)">[20,0,0,0,20,0,0,0,20]</cell>
  </system>""")
    ttime = "" if total_time is None else f"<total_time>{total_time}</total_time>"
    return f"""<simulation verbosity="quiet" interleave="{interleave}"
    safe_stride="{safe_stride}" floatformat="%24.16e">
  <output prefix="sim">
    {outputs}
  </output>
  <total_steps>{nsteps}</total_steps>
  {ttime}
  <prng><seed>31415</seed></prng>
  <ffdirect name="harmonic">
    <pes>harmonic</pes>
    <parameters>{{k1: 0.1}}</parameters>
  </ffdirect>
  {"".join(systems)}
  {smotion}
</simulation>"""


# ---------------------------------------------------------------- running


@pytest.fixture
def dirs(tmp_path, monkeypatch):
    """The directories of the reference and of the interrupted run. The
    tests run in the latter."""

    ref, run = tmp_path / "ref", tmp_path / "run"
    ref.mkdir()
    run.mkdir()
    monkeypatch.chdir(run)
    monkeypatch.setattr(softexit_module, "SOFTEXITLATENCY", 0.01)
    softexit.reset()
    yield ref, run
    softexit.reset()


def load(xml=None):
    """Loads a simulation from an input, or from RESTART."""

    if xml is not None:
        return Simulation.load_from_xml(io.StringIO(xml))
    with open("RESTART") as restart_file:
        return Simulation.load_from_xml(restart_file)


def restart_step(filename="RESTART"):
    """The step a checkpoint continues from (also when it is the last one,
    which Simulation refuses to load)."""

    with open(filename) as checkpoint:
        steps = re.findall(r"<step>\s*(\d+)\s*</step>", checkpoint.read())
    assert len(steps) <= 1, steps
    return int(steps[0]) if steps else 0  # step 0 is the default, not written


def run_reference(ref, xml):
    """Makes the uninterrupted run."""

    cwd = os.getcwd()
    os.chdir(ref)
    try:
        finish(load(xml))
    finally:
        os.chdir(cwd)


def run_to_end(nsteps):
    """Restarts from RESTART until the run is complete."""

    for _ in range(20):
        if restart_step() >= nsteps:
            return
        finish(load())
    pytest.fail("The restarted run does not get to the end")


class MidWrite:
    """A stream that calls action after each of its writes, i.e. also
    between the writes of a row or of a frame."""

    def __init__(self, stream, action):
        self._stream = stream
        self._action = action

    def write(self, data):
        result = self._stream.write(data)
        self._action()
        return result

    def __getattr__(self, name):
        return getattr(self._stream, name)


def once(simulation, step, action):
    """Wraps action to fire once: during `step`, or at the first chance if
    `step` is None. Returns the wrapper, and a list telling whether it fired."""

    fired = []

    def guarded():
        if not fired and (step is None or simulation.step == step):
            fired.append(simulation.step)
            action()

    return guarded, fired


def arm(simulation, kind, step):
    """Sets up an exit at one of the injection points of test_softexit.py."""

    if kind == "step":  # requested in a step, made while it waits for forces
        action, fired = once(simulation, step, request_exit)
        inject_before(simulation.syslist[0].motion, "step", action)
    elif kind == "request_in_outputs":  # made at the top of the next step
        action, fired = once(simulation, step, request_exit)
        inject_after(simulation.outputs[0], "write", action)
    elif kind == "between_outputs":  # made after some outputs were written
        action, fired = once(simulation, step, make_exit)
        inject_after(simulation.outputs[0], "write", action)
    elif kind == "mid_row":  # made in the middle of a row of properties
        action, fired = once(simulation, step, make_exit)
        simulation.outputs[0].out = MidWrite(simulation.outputs[0].out, action)
    else:
        raise ValueError(kind)
    return fired


def interrupt(simulation, kind, step):
    """Runs a simulation up to an exit, and returns the step of RESTART."""

    fired = arm(simulation, kind, step)
    finish(simulation)
    assert fired, f"the exit ({kind} at step {step}) was never injected"
    return restart_step()


def stored_step(step, kind, safe_stride):
    """The step of the last stored state, for an exit injected at `step`."""

    if kind == "request_in_outputs":
        # review fix 3: the step is completed before the exit is made, which
        # then starts from the current state, without a rollback
        return step + 1
    return step // safe_stride * safe_stride


# --------------------------------------------------------------- checking


def is_output(name):
    return not (
        name.startswith("#")
        or name in ("RESTART", "EXIT")
        or ".chk" in name
        or os.path.isdir(name)
    )


def steps_in(path):
    """The steps of the rows or frames of an output, which must be whole."""

    name = os.path.basename(path)
    with open(path) as stream:
        lines = stream.readlines()
    if name.endswith(".xyz"):
        size = NATOMS + 2
        assert len(lines) % size == 0, f"{name}: partial frame ({len(lines)} lines)"
        steps = []
        for i in range(0, len(lines), size):
            assert lines[i].split() == [str(NATOMS)], f"{name}: bad frame at line {i}"
            assert "Step:" in lines[i + 1], f"{name}: bad frame at line {i}"
            steps.append(int(lines[i + 1].split("Step:")[1].split()[0]))
            for line in lines[i + 2 : i + size]:
                assert len(line.split()) == 4, f"{name}: partial atom line {line!r}"
                assert line.endswith("\n"), f"{name}: partial atom line {line!r}"
        return steps
    rows = [line for line in lines if not line.startswith("#")]
    ncols = {len(row.split()) for row in rows}
    assert len(ncols) <= 1, f"{name}: partial row (columns {ncols})"
    assert all(row.endswith("\n") for row in rows), f"{name}: partial row"
    return [int(float(row.split()[0])) for row in rows]


def stride_of(name):
    for key, stride in STRIDES.items():
        if key in name:
            return stride
    return None


def assert_rolled_back(step, names):
    """Checks that the outputs hold the steps up to `step` and no later one,
    right after a soft exit whose RESTART continues from `step`."""

    for name in names:
        stride = stride_of(name)
        if not os.path.exists(name):
            # removed with the initial frame: the restart writes it again
            assert step == 0, f"{name} is missing, RESTART is at step {step}"
            continue
        steps = steps_in(name)
        if stride is None:  # e.g. the swap file: only an upper bound is known
            assert all(s < step for s in steps), (name, steps, step)
            assert steps == sorted(set(steps)), (name, steps)
            continue
        expected = list(range(0, step + 1, stride))
        if step == 0 and not steps:
            continue
        assert steps == expected, (
            f"{name} has steps {steps} after an exit whose RESTART is at step "
            f"{step}: expected {expected}"
        )


def compare(path, reference):
    """Compares two outputs: the same lines and words, and the same numbers
    to RTOL. Returns a description of the first difference, or None.

    The outputs are written with 17 digits. They are not compared byte by
    byte because the values after a restart differ in the last bits (1e-16
    relative) from those of a continuous run, also when the exit is made at
    a step boundary with nothing to roll back: RESTART holds the state with
    all its digits, but derived quantities are computed again from it.
    """

    with open(path) as stream, open(reference) as reference_stream:
        lines, expected = stream.readlines(), reference_stream.readlines()
    if len(lines) != len(expected):
        return f"{len(lines)} lines instead of {len(expected)}"
    for number, (line, wanted) in enumerate(zip(lines, expected), 1):
        if line == wanted:
            continue
        words, wanted_words = line.split(), wanted.split()
        if len(words) != len(wanted_words):
            return f"line {number} is {line!r} instead of {wanted!r}"
        for word, wanted_word in zip(words, wanted_words):
            if word == wanted_word:
                continue
            try:
                close = np.isclose(
                    float(word), float(wanted_word), rtol=RTOL, atol=ATOL
                )
            except ValueError:
                close = False
            if not close:
                return f"line {number} has {word} instead of {wanted_word}"
    return None


def assert_criterion(ref, run, nsteps):
    """The criterion C: each expected step once and in order, whole rows, and
    the same values as the uninterrupted reference."""

    names = sorted(name for name in os.listdir(ref) if is_output(name))
    assert names, "the reference run wrote no outputs"
    problems = []
    for name in names:
        if not os.path.exists(run / name):
            problems.append(f"{name}: missing")
            continue
        steps, reference = steps_in(run / name), steps_in(ref / name)
        stride = stride_of(name)
        if stride is not None:
            # the reference itself has what is expected
            assert reference == list(range(0, nsteps + 1, stride)), name
        if steps != reference:
            problems.append(f"{name}: steps {steps}, expected {reference}")
        else:
            difference = compare(run / name, ref / name)
            if difference:
                problems.append(f"{name}: right steps, but {difference}")
    assert not problems, "\n".join(problems)
    return names


# ------------------------------------------------------------------ tests


@pytest.mark.parametrize("safe_stride", [1, 3, 5])
@pytest.mark.parametrize(
    "kind", ["step", "request_in_outputs", "between_outputs", "mid_row"]
)
def test_rollback_with_safe_stride(dirs, safe_stride, kind):
    """1. An exit at step 7 goes back to the last stored state, removing the
    rows written since then, and the restart satisfies C."""

    ref, run = dirs
    nsteps, exit_step = 12, 7
    xml = simulation_xml(nsteps, safe_stride, outputs=OUTPUTS + CHECKPOINT)
    run_reference(ref, xml)
    names = [name for name in os.listdir(ref) if is_output(name)]

    restart = interrupt(load(xml), kind, exit_step)
    assert restart == stored_step(exit_step, kind, safe_stride)
    assert_rolled_back(restart, names)

    run_to_end(nsteps)
    assert_criterion(ref, run, nsteps)


@pytest.mark.parametrize("interleave", ["false", "true"])
@pytest.mark.parametrize("safe_stride", [1, 3])
def test_rollback_multibead_multisystem(dirs, interleave, safe_stride):
    """2. An exit made while the second system writes its outputs, after the
    first one wrote its own, takes all the files of both systems back to the
    same step."""

    ref, run = dirs
    nsteps, exit_step = 10, 5  # the three outputs are all active at step 5
    xml = simulation_xml(
        nsteps,
        safe_stride,
        interleave,
        temperatures=(300, 350),
        prefixes=("a", "b"),
    )
    run_reference(ref, xml)
    names = sorted(name for name in os.listdir(ref) if is_output(name))
    # properties, and 4 beads of positions and of velocities, per system
    assert len(names) == 2 * (1 + 4 + 4)

    simulation = load(xml)
    second = simulation.syslist[1]
    written = []
    for output in simulation.outputs:
        if output.system is simulation.syslist[0]:
            inject_after(output, "write", lambda o=output: written.append(o.filename))
    positions = [
        o for o in simulation.outputs if o.system is second and ".pos" in o.filename
    ][0]
    action, fired = once(simulation, exit_step, make_exit)
    inject_after(positions, "write", action)
    first_outputs = []

    def record():
        if simulation.step == exit_step and not first_outputs:
            first_outputs.extend(written)

    inject_before(positions, "write", record)
    finish(simulation)
    assert fired
    if interleave == "false":
        # the first system had written properties and positions
        assert any(".pos" in name for name in first_outputs[-4:])

    restart = restart_step()
    assert restart == exit_step // safe_stride * safe_stride
    assert_rolled_back(restart, names)

    run_to_end(nsteps)
    assert_criterion(ref, run, nsteps)


@pytest.mark.parametrize("safe_stride", [1, 3])
def test_rollback_on_appended_files(dirs, safe_stride):
    """3. Exits of runs that append to the outputs of the previous ones,
    including an exit in the very first step after a restart."""

    ref, run = dirs
    nsteps = 12
    xml = simulation_xml(nsteps, safe_stride)
    run_reference(ref, xml)
    names = [name for name in os.listdir(ref) if is_output(name)]

    restart = interrupt(load(xml), "step", 4)
    assert restart == stored_step(4, "step", safe_stride)
    assert_rolled_back(restart, names)

    # in the first step after the restart, in the middle of a row
    first = restart
    restart = interrupt(load(), "mid_row", None)
    assert restart == first
    assert_rolled_back(restart, names)

    # between the outputs of the first step after the restart
    restart = interrupt(load(), "between_outputs", None)
    assert restart == first
    assert_rolled_back(restart, names)

    restart = interrupt(load(), "mid_row", 8)
    assert restart == stored_step(8, "mid_row", safe_stride)
    assert_rolled_back(restart, names)

    restart = interrupt(load(), "request_in_outputs", 9)
    assert restart == stored_step(9, "request_in_outputs", safe_stride)
    assert_rolled_back(restart, names)

    run_to_end(nsteps)
    assert_criterion(ref, run, nsteps)


@pytest.mark.parametrize(
    "kind,step",
    [
        ("mid_row", -1),
        ("between_outputs", -1),
        ("request_in_outputs", -1),
        ("step", 0),
        ("mid_row", 0),
    ],
)
def test_exit_in_first_step_and_initial_outputs(dirs, kind, step):
    """4. An exit while the initial outputs are written, or during step 0,
    restarts from step 0, and the initial frame ends up written once."""

    ref, run = dirs
    nsteps = 6
    xml = simulation_xml(nsteps, safe_stride=2)
    run_reference(ref, xml)
    names = [name for name in os.listdir(ref) if is_output(name)]

    restart = interrupt(load(xml), kind, step)
    assert restart == 0
    for name in names:
        if not os.path.exists(name):
            continue
        steps = steps_in(name)  # also checks there is no partial row or frame
        if step == -1 and kind != "request_in_outputs":
            # the initial frame was being written: nothing of it is left
            # (a trajectory not written yet is still the empty file of bind)
            assert steps == [], (name, steps)
        else:
            assert steps == [0], (name, steps)

    run_to_end(nsteps)
    assert_criterion(ref, run, nsteps)


def test_restart_matches_uninterrupted_run(dirs):
    """5. A 50-step Langevin run interrupted at three random points, with
    safe_stride=4, gives the same bytes as the uninterrupted one."""

    ref, run = dirs
    nsteps, safe_stride = 50, 4
    xml = simulation_xml(nsteps, safe_stride)
    run_reference(ref, xml)
    names = [name for name in os.listdir(ref) if is_output(name)]

    rng = np.random.default_rng(617)
    steps = sorted(int(s) for s in rng.choice(np.arange(2, 48, 3), 3, replace=False))
    kinds = rng.choice(
        ["step", "request_in_outputs", "between_outputs", "mid_row"], 3
    ).tolist()

    simulation = load(xml)
    for step, kind in zip(steps, kinds):
        restart = interrupt(simulation, kind, step)
        assert restart == stored_step(step, kind, safe_stride), (step, kind)
        assert_rolled_back(restart, names)
        simulation = load()
    finish(simulation)
    assert restart_step() == nsteps
    assert_criterion(ref, run, nsteps)


def checkpoint_run(dirs):
    """safe_stride=10, a checkpoint every 5 steps, and an exit in step 17."""

    ref, run = dirs
    nsteps = 24
    xml = simulation_xml(nsteps, safe_stride=10, outputs=OUTPUTS + CHECKPOINT)
    run_reference(ref, xml)
    names = [name for name in os.listdir(ref) if is_output(name)]
    restart = interrupt(load(xml), "step", 17)
    # review fix 1: the checkpoint written at step 15 is also a stored state
    assert restart == 15
    assert_rolled_back(restart, names)
    return ref, run, nsteps


def test_checkpoint_between_stores(dirs):
    """6a. With a checkpoint written between two stored states, the restart
    from RESTART satisfies C."""

    ref, run, nsteps = checkpoint_run(dirs)
    run_to_end(nsteps)
    assert_criterion(ref, run, nsteps)


def test_checkpoint_between_stores_checkpoint_file(dirs):
    """6b. The checkpoint file left by the exit is consistent with the
    outputs: restarting from it also satisfies C (risk 1 of the review)."""

    ref, run, nsteps = checkpoint_run(dirs)
    chk_step = restart_step("sim.chk")
    out_steps = steps_in("sim.out")
    print(f"checkpoint at step {chk_step}, sim.out up to step {out_steps[-1]}")

    with open("sim.chk") as checkpoint:
        simulation = Simulation.load_from_xml(checkpoint)
    finish(simulation)
    assert_criterion(ref, run, nsteps)


@pytest.mark.parametrize("interleave", ["false", "true"])
@pytest.mark.parametrize("kind", ["step", "after_swap", "between_outputs"])
def test_remd_swapfile_rollback(dirs, interleave, kind):
    """7. The swap file of a replica exchange run is rolled back with the
    other outputs: after the restart it has each swap once."""

    ref, run = dirs
    nsteps, safe_stride, exit_step = 12, 3, 7
    xml = simulation_xml(
        nsteps,
        safe_stride,
        interleave,
        nbeads=1,
        temperatures=(300, 305),
        prefixes=("a", "b"),
        smotion=REMD,
    )
    run_reference(ref, xml)
    names = [name for name in os.listdir(ref) if is_output(name)]
    assert "sim.remd_idx" in names
    swaps = steps_in(ref / "sim.remd_idx")
    # swaps are written both before and after the state that is restored
    assert any(6 <= s <= exit_step for s in swaps), swaps
    assert any(s < 6 for s in swaps) and any(s > exit_step for s in swaps), swaps

    simulation = load(xml)
    if kind == "after_swap":
        # after the swap file is written, before the outputs of the step
        action, fired = once(simulation, exit_step, make_exit)
        inject_after(simulation.smotion, "step", action)
        finish(simulation)
        assert fired
        restart = restart_step()
    else:
        restart = interrupt(simulation, kind, exit_step)
    assert restart == 6
    assert_rolled_back(restart, names)

    run_to_end(nsteps)
    assert_criterion(ref, run, nsteps)
    assert steps_in("sim.remd_idx") == swaps


@pytest.mark.parametrize("safe_stride", [1, 3])
@pytest.mark.parametrize(
    "point", ["between_beads", "within_frame", "between_trajectories"]
)
def test_exit_during_output_write(dirs, safe_stride, point):
    """10. An exit made between the beads of a trajectory frame, between the
    writes of one xyz frame, or between two trajectories, leaves no partial
    frame."""

    ref, run = dirs
    nsteps, exit_step = 12, 5  # positions and velocities are written at step 5
    xml = simulation_xml(nsteps, safe_stride)
    run_reference(ref, xml)
    names = [name for name in os.listdir(ref) if is_output(name)]

    simulation = load(xml)
    positions = simulation.outputs[1]
    assert ".pos" in positions.filename and len(positions.out) == 4
    action, fired = once(simulation, exit_step, make_exit)
    if point == "between_beads":
        # after the frames of beads 0 and 1 are written, before the others
        original = positions.write_traj

        def write_traj(data, what, stream, b=0, **kwargs):
            result = original(data, what, stream, b, **kwargs)
            if b == 1:
                action()
            return result

        positions.write_traj = write_traj
    elif point == "within_frame":
        # after the first write (the header) of the frame of bead 2
        positions.out[2] = MidWrite(positions.out[2], action)
    else:
        inject_after(positions, "write", action)
    finish(simulation)
    assert fired

    restart = restart_step()
    assert restart == exit_step // safe_stride * safe_stride
    assert_rolled_back(restart, names)

    run_to_end(nsteps)
    assert_criterion(ref, run, nsteps)


def test_walltime_at_step_boundary(dirs, capsys):
    """11a. The wall-clock time expires between two steps, with
    safe_stride=5: C holds after the restart. The steps that were complete
    and get discarded are reported (risk 3 of the review)."""

    ref, run = dirs
    nsteps, safe_stride, last_step = 14, 5, 7
    xml = simulation_xml(nsteps, safe_stride)
    run_reference(ref, xml)
    names = [name for name in os.listdir(ref) if is_output(name)]

    simulation = load(xml)

    def expire():
        softexit.timeout = time.time() - 1.0
        softexit._next_check = 0.0

    # after the last output of the step, i.e. at the step boundary
    action, fired = once(simulation, last_step, expire)
    inject_after(simulation.outputs[-1], "write", action)
    finish(simulation)
    assert fired

    restart = restart_step()
    completed = last_step + 1
    with capsys.disabled():
        print(
            f"\n[walltime at boundary] completed steps: {completed}, RESTART at "
            f"step {restart}, discarded: {completed - restart}"
        )
    assert restart <= completed
    assert_rolled_back(restart, names)

    run_to_end(nsteps)
    assert_criterion(ref, run, nsteps)


def test_walltime_total_time_end_to_end(dirs, capsys):
    """11b. The same with a real <total_time>, and steps made slow."""

    ref, run = dirs
    nsteps, safe_stride = 30, 5
    run_reference(ref, simulation_xml(nsteps, safe_stride))
    names = [name for name in os.listdir(ref) if is_output(name)]

    # the restarts have the same limit, and make fast steps
    simulation = load(simulation_xml(nsteps, safe_stride, total_time=0.6))
    started = []
    inject_before(
        simulation.syslist[0].motion,
        "step",
        lambda: (started.append(simulation.step), time.sleep(0.05)),
    )
    finish(simulation)

    restart = restart_step()
    with capsys.disabled():
        print(
            f"\n[total_time] last step started: {started[-1]}, RESTART at step "
            f"{restart}, sim.out up to {steps_in('sim.out')[-1]}"
        )
    assert 0 < restart < nsteps, "the wall-clock limit did not interrupt the run"
    assert_rolled_back(restart, names)

    run_to_end(nsteps)
    assert_criterion(ref, run, nsteps)

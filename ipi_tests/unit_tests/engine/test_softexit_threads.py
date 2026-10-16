"""Soft exits from the different thread roles in a two-system socket run.

Each simulation has seven threads: main, four shared system/output workers,
one forcefield polling thread, and the soft-exit monitor. The harmonic driver
runs in a separate process. Barriers ensure the workers actually exist;
interruptions are tied to simulation steps rather than elapsed time.
"""

import copy
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import threading
import time
import uuid
import xml.etree.ElementTree as ET

import numpy as np
import pytest

import ipi.utils.softexit as softexit_module
from ipi.engine.simulation import Simulation
from ipi.utils.softexit import softexit

from ipi_tests.unit_tests.engine.test_softexit import EXIT_STEP, NSTEPS, simulation_xml


def socket_input(address):
    root = ET.fromstring(simulation_xml("true"))
    root.remove(root.find("ffdirect"))
    forcefield = ET.SubElement(root, "ffsocket", name="harmonic", mode="unix")
    ET.SubElement(forcefield, "address").text = address
    ET.SubElement(forcefield, "latency").text = "0.001"
    ET.SubElement(forcefield, "timeout").text = "10"
    first = root.find("system")
    first.set("prefix", "a")
    second = copy.deepcopy(first)
    second.set("prefix", "b")
    root.append(second)
    return ET.tostring(root, encoding="unicode")


def wait_for(condition):
    deadline = time.monotonic() + 20
    while not condition():
        if time.monotonic() > deadline:
            raise TimeoutError("Timed out waiting for the coordinated soft exit")
        time.sleep(0.001)


def run_simulation(input_file, origin):
    """Subprocess entry point, using the same hard-exit mode as bin/i-pi."""
    softexit_module.SOFTEXITLATENCY = 0.01
    with open(input_file) as stream:
        simulation = Simulation.load_from_xml(stream)
    softexit.hard_exit = True
    Path("ready").touch()

    def count_cleanup():
        with open("cleanup_calls", "a") as stream:
            stream.write("cleanup\n")

    softexit.register_function(count_cleanup)
    original_trigger = softexit.trigger

    def recorded_trigger(*args, **kwargs):
        with open("trigger_threads", "a") as stream:
            stream.write(threading.current_thread().name + "\n")
        original_trigger(*args, **kwargs)

    softexit.trigger = recorded_trigger
    barrier = threading.Barrier(4, timeout=20)

    def output_start():
        if simulation.step == -1:
            ticket = barrier.wait()
            if ticket == 0:
                threads = threading.enumerate()
                Path("threads.json").write_text(
                    json.dumps([thread.name for thread in threads])
                )

    fired = False

    def request_exit():
        nonlocal fired
        if simulation.step != EXIT_STEP or fired:
            return
        fired = True
        Path("origin").write_text(threading.current_thread().name)
        if origin == "signal":
            os.kill(os.getpid(), signal.SIGTERM)
            wait_for(lambda: softexit.triggered)
        elif origin == "exit-file":
            Path("EXIT").touch()
            wait_for(lambda: softexit.triggered)
        else:
            softexit.trigger(message="Test exit from " + origin)

    def wrap(obj, name, before=None, after=None):
        original = getattr(obj, name)

        def wrapped(*args, **kwargs):
            if before is not None:
                before()
            result = original(*args, **kwargs)
            if after is not None:
                after()
            return result

        setattr(obj, name, wrapped)

    for output in simulation.outputs:
        wrap(output, "write", before=output_start)

    if origin == "main":
        wrap(simulation, "run_step", before=request_exit)
    elif origin == "system":
        wrap(simulation.syslist[0].motion, "step", before=request_exit)
    elif origin == "forcefield":
        forcefield = simulation.fflist["harmonic"]

        def exit_during_force_request():
            if forcefield.requests:
                request_exit()

        wrap(forcefield, "poll", before=exit_during_force_request)
    elif origin in ("properties", "signal", "exit-file"):
        wrap(simulation.outputs[0], "write", after=request_exit)
    elif origin == "trajectory":
        wrap(simulation.outputs[2], "write", after=request_exit)

    try:
        simulation.run()
    finally:
        softexit.trigger()


def run_with_driver(tmp_path, input_file, origin, address):
    """Keep hard exits and signals out of pytest, and reap both processes."""
    root = Path(__file__).resolve().parents[3]
    env = dict(os.environ)
    env["PYTHONPATH"] = str(root)
    env["PYTHONUNBUFFERED"] = "1"
    server = driver = None
    (tmp_path / "ready").unlink(missing_ok=True)
    with open(tmp_path / (origin + ".log"), "w") as log:
        try:
            server = subprocess.Popen(
                [sys.executable, __file__, input_file, origin],
                cwd=tmp_path,
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
            )
            wait_for(lambda: (tmp_path / "ready").exists() or server.poll() is not None)
            assert server.poll() is None, (tmp_path / (origin + ".log")).read_text()
            driver = subprocess.Popen(
                [
                    sys.executable,
                    str(root / "drivers/py/driver.py"),
                    "-u",
                    "-a",
                    address,
                    "-m",
                    "harmonic",
                    "-o",
                    "0.1",
                ],
                cwd=tmp_path,
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
            )
            server.wait(timeout=30)
            driver.wait(timeout=10)
            assert server.returncode == 0, (tmp_path / (origin + ".log")).read_text()
            assert driver.returncode == 0, (tmp_path / (origin + ".log")).read_text()
        finally:
            for process in (server, driver):
                if process is not None and process.poll() is None:
                    process.kill()
                    process.wait(timeout=5)
            Path("/tmp/ipi_" + address).unlink(missing_ok=True)


@pytest.mark.parametrize(
    "origin",
    ["main", "system", "properties", "trajectory", "forcefield", "signal", "exit-file"],
)
def test_seven_thread_soft_exit_and_restart(tmp_path, origin):
    address = "softexit_" + uuid.uuid4().hex
    (tmp_path / "input.xml").write_text(socket_input(address))
    run_with_driver(tmp_path, "input.xml", origin, address)

    names = json.loads((tmp_path / "threads.json").read_text())
    assert len(names) == 7, names
    assert "MainThread" in names
    assert "softexit" in names
    assert "poll_harmonic" in names
    assert sum(name.startswith("ThreadPoolExecutor-") for name in names) == 4
    caller = (tmp_path / "origin").read_text()
    if origin == "main":
        assert caller == "MainThread"
    elif origin == "forcefield":
        assert caller == "poll_harmonic"
    else:
        assert caller.startswith("ThreadPoolExecutor-")
    triggers = (tmp_path / "trigger_threads").read_text().splitlines()
    if origin in ("signal", "exit-file"):
        assert "softexit" in triggers
    else:
        assert caller in triggers
    assert (tmp_path / "cleanup_calls").read_text() == "cleanup\n"
    with open(tmp_path / "RESTART") as restart:
        saved = Simulation.load_from_xml(restart, read_only=True)
    assert 0 < saved.step < NSTEPS

    (tmp_path / "EXIT").unlink(missing_ok=True)
    run_with_driver(tmp_path, "RESTART", "resume", address)
    assert (tmp_path / "cleanup_calls").read_text() == "cleanup\ncleanup\n"
    for prefix in ("a", "b"):
        data = np.loadtxt(tmp_path / (prefix + "_sim.out"))
        assert data.shape == (NSTEPS + 1, 2)
        np.testing.assert_array_equal(data[:, 0], np.arange(NSTEPS + 1))
        lines = (tmp_path / (prefix + "_sim.pos_0.xyz")).read_text().splitlines()
        assert len(lines) == 3 * (NSTEPS + 1)
        steps = [int(line.split("Step:")[1].split()[0]) for line in lines[1::3]]
        assert steps == list(range(NSTEPS + 1))


if __name__ == "__main__":
    run_simulation(*sys.argv[1:])

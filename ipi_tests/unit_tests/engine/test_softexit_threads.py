"""Soft exit requested by the polling thread of a socket forcefield.

The exit is requested while a force request is pending, as when a client
fails, so the cleanup runs in the thread that the forcefield would otherwise
join when it stops. i-PI and the driver run as separate processes, with the
hard exit used by bin/i-pi.
"""

import os
from pathlib import Path
import subprocess
import sys
import threading
import time
import uuid
import xml.etree.ElementTree as ET

import ipi.utils.softexit as softexit_module
from ipi.engine.simulation import Simulation
from ipi.utils.softexit import softexit

from ipi_tests.unit_tests.engine.test_softexit import (
    EXIT_STEP,
    NSTEPS,
    assert_outputs_continuous,
    simulation_xml,
)


def socket_input(address):
    """Returns the test simulation, with forces from a socket client."""

    root = ET.fromstring(simulation_xml("false"))
    root.remove(root.find("ffdirect"))
    forcefield = ET.SubElement(root, "ffsocket", name="harmonic", mode="unix")
    ET.SubElement(forcefield, "address").text = address
    ET.SubElement(forcefield, "latency").text = "0.001"
    ET.SubElement(forcefield, "timeout").text = "10"
    return ET.tostring(root, encoding="unicode")


def wait_for(condition):
    """Waits until condition() is true."""

    deadline = time.monotonic() + 20
    while not condition():
        if time.monotonic() > deadline:
            raise TimeoutError("Timed out waiting for i-PI to start")
        time.sleep(0.001)


def run_simulation(input_file, request_exit):
    """Subprocess entry point, using the same hard-exit mode as bin/i-pi."""

    softexit_module.SOFTEXITLATENCY = 0.01
    with open(input_file) as stream:
        simulation = Simulation.load_from_xml(stream)
    softexit.hard_exit = True

    if request_exit == "true":
        forcefield = simulation.fflist["harmonic"]
        poll = forcefield.poll

        def poll_and_exit():
            if simulation.step == EXIT_STEP and forcefield.requests:
                Path("origin").write_text(threading.current_thread().name)
                softexit.trigger(message="Test exit from the polling thread")
            poll()

        forcefield.poll = poll_and_exit

    Path("ready").touch()
    try:
        simulation.run()
    finally:
        softexit.trigger()


def run_with_driver(tmp_path, input_file, request_exit, address):
    """Runs i-PI and the driver, and checks that both exit cleanly."""

    root = Path(__file__).resolve().parents[3]
    env = dict(os.environ)
    env["PYTHONPATH"] = str(root)
    env["PYTHONUNBUFFERED"] = "1"
    logfile = tmp_path / "log"
    server = driver = None
    (tmp_path / "ready").unlink(missing_ok=True)
    with open(logfile, "w") as log:
        try:
            server = subprocess.Popen(
                [sys.executable, __file__, input_file, request_exit],
                cwd=tmp_path,
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
            )
            wait_for(lambda: (tmp_path / "ready").exists() or server.poll() is not None)
            assert server.poll() is None, logfile.read_text()
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
            # the driver only exits if i-PI tells it to
            driver.wait(timeout=10)
            assert server.returncode == 0, logfile.read_text()
            assert driver.returncode == 0, logfile.read_text()
        finally:
            for process in (server, driver):
                if process is not None and process.poll() is None:
                    process.kill()
                    process.wait(timeout=5)
            Path("/tmp/ipi_" + address).unlink(missing_ok=True)


def test_exit_from_the_polling_thread(tmp_path, monkeypatch):
    """The driver is told to exit, and RESTART continues the outputs."""

    monkeypatch.chdir(tmp_path)
    address = "softexit_" + uuid.uuid4().hex
    (tmp_path / "input.xml").write_text(socket_input(address))
    run_with_driver(tmp_path, "input.xml", "true", address)

    assert (tmp_path / "origin").read_text() == "poll_harmonic"
    with open(tmp_path / "RESTART") as restart:
        assert 0 < Simulation.load_from_xml(restart, read_only=True).step < NSTEPS

    run_with_driver(tmp_path, "RESTART", "false", address)
    assert_outputs_continuous()


if __name__ == "__main__":
    run_simulation(*sys.argv[1:])

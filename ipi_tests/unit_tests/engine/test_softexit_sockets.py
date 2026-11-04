"""Soft exits of a run whose forces come from a socket client.

i-PI and the clients run as separate processes, with the hard exit used by
bin/i-pi. The exit is requested while a force request is pending, as when a
client fails or hangs, so the cleanup closes the socket that is being
polled.
"""

import os
from pathlib import Path
import signal
import socket
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

from ipi_tests.unit_tests.engine.test_softexit import (
    EXIT_STEP,
    NSTEPS,
    assert_outputs_continuous,
    simulation_xml,
)

ROOT = Path(__file__).resolve().parents[3]


def socket_input(address):
    """Returns the test simulation, with forces from a socket client."""

    root = ET.fromstring(simulation_xml("false"))
    root.remove(root.find("ffdirect"))
    forcefield = ET.SubElement(root, "ffsocket", name="harmonic", mode="unix")
    ET.SubElement(forcefield, "address").text = address
    ET.SubElement(forcefield, "latency").text = "0.001"
    ET.SubElement(forcefield, "timeout").text = "60"
    return ET.tostring(root, encoding="unicode")


def wait_for(condition):
    """Waits until condition() is true."""

    deadline = time.monotonic() + 20
    while not condition():
        if time.monotonic() > deadline:
            raise TimeoutError("Timed out waiting for i-PI")
        time.sleep(0.001)


def run_simulation(input_file, request_exit):
    """Subprocess entry point, using the same hard-exit mode as bin/i-pi."""

    softexit_module.SOFTEXITLATENCY = 0.01
    with open(input_file) as stream:
        simulation = Simulation.load_from_xml(stream)
    softexit.hard_exit = True
    # records the thread that makes the exit
    softexit.register_function(
        lambda: Path("origin").write_text(threading.current_thread().name)
    )

    if request_exit == "true":
        forcefield = simulation.fflist["harmonic"]
        poll = forcefield.poll

        def poll_and_exit(idle=False):
            if simulation.step == EXIT_STEP and forcefield.requests:
                softexit.request("Test exit while polling")
            return poll(idle)

        forcefield.poll = poll_and_exit

    Path("ready").touch()
    try:
        simulation.run()
    finally:
        softexit.trigger()


def start_ipi(tmp_path, input_file, request_exit, log):
    """Starts i-PI, and returns its process once it listens."""

    env = dict(os.environ, PYTHONPATH=str(ROOT), PYTHONUNBUFFERED="1")
    (tmp_path / "ready").unlink(missing_ok=True)
    server = subprocess.Popen(
        [sys.executable, __file__, input_file, request_exit],
        cwd=tmp_path,
        env=env,
        stdout=log,
        stderr=subprocess.STDOUT,
    )
    wait_for(lambda: (tmp_path / "ready").exists() or server.poll() is not None)
    return server


def stop(*processes):
    """Kills the processes that are still running."""

    for process in processes:
        if process is not None and process.poll() is None:
            process.kill()
            process.wait(timeout=5)


def run_with_driver(tmp_path, input_file, request_exit, address):
    """Runs i-PI and the driver, and checks that both exit cleanly."""

    logfile = tmp_path / "log"
    server = driver = None
    with open(logfile, "w") as log:
        try:
            server = start_ipi(tmp_path, input_file, request_exit, log)
            assert server.poll() is None, logfile.read_text()
            driver = subprocess.Popen(
                [
                    sys.executable,
                    str(ROOT / "drivers/py/driver.py"),
                    "-u",
                    "-a",
                    address,
                    "-m",
                    "harmonic",
                    "-o",
                    "0.1",
                ],
                cwd=tmp_path,
                env=dict(os.environ, PYTHONPATH=str(ROOT)),
                stdout=log,
                stderr=subprocess.STDOUT,
            )
            server.wait(timeout=30)
            # the driver only exits if i-PI tells it to
            driver.wait(timeout=10)
            assert server.returncode == 0, logfile.read_text()
            assert driver.returncode == 0, logfile.read_text()
        finally:
            stop(server, driver)
            Path("/tmp/ipi_" + address).unlink(missing_ok=True)


def test_exit_while_polling(tmp_path, monkeypatch):
    """The driver is told to exit, and RESTART continues the outputs."""

    monkeypatch.chdir(tmp_path)
    address = "softexit_" + uuid.uuid4().hex
    (tmp_path / "input.xml").write_text(socket_input(address))
    run_with_driver(tmp_path, "input.xml", "true", address)

    # the exit is made by the main thread
    assert (tmp_path / "origin").read_text() == "MainThread"
    with open(tmp_path / "RESTART") as restart:
        assert 0 < Simulation.load_from_xml(restart, read_only=True).step < NSTEPS

    run_with_driver(tmp_path, "RESTART", "false", address)
    assert_outputs_continuous()


def message(name):
    """A header of the i-PI protocol."""

    return name.upper().encode().ljust(12)


def receive(peer, length):
    """Receives exactly `length` bytes."""

    data = b""
    while len(data) < length:
        chunk = peer.recv(length - len(data))
        assert chunk, "i-PI closed the connection"
        data += chunk
    return data


def hang(address, nreplies):
    """Connects to i-PI as a client of a one-atom system, returns zero
    forces for `nreplies` requests, and then leaves a request unanswered."""

    peer = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    peer.settimeout(20)
    peer.connect("/tmp/ipi_" + address)
    for i in range(nreplies + 1):
        assert receive(peer, 12) == message("status")
        peer.sendall(message("ready"))
        assert receive(peer, 12) == message("posdata")
        # cell, inverse cell, number of atoms, positions, and GETFORCE
        assert receive(peer, 72 + 72 + 4 + 24 + 12).endswith(message("getforce"))
        if i == nreplies:
            return peer
        peer.sendall(
            message("forceready")
            + np.float64(0.0).tobytes()
            + np.int32(1).tobytes()
            + np.zeros(3 + 9).tobytes()
            + np.int32(0).tobytes()
        )


@pytest.mark.parametrize("held", [False, True], ids=["step", "outputs"])
def test_exit_while_a_client_hangs(tmp_path, monkeypatch, held):
    """SIGTERM ends a run whose client hangs, writing RESTART and telling the
    client to exit. A client that hangs while the outputs hold the exit back
    (here, the initial ones) needs a second SIGTERM."""

    monkeypatch.chdir(tmp_path)
    address = "softexit_" + uuid.uuid4().hex
    (tmp_path / "input.xml").write_text(socket_input(address))
    logfile = tmp_path / "log"
    server = peer = None
    with open(logfile, "w") as log:
        try:
            server = start_ipi(tmp_path, "input.xml", "false", log)
            # the initial outputs need the forces of the first request
            peer = hang(address, nreplies=0 if held else 1)
            time.sleep(0.2)
            server.send_signal(signal.SIGTERM)
            if held:
                time.sleep(1.0)
                assert server.poll() is None, logfile.read_text()
                server.send_signal(signal.SIGTERM)
            server.wait(timeout=10)
            assert server.returncode == 0, logfile.read_text()
            assert message("exit") in peer.recv(4096)
        finally:
            stop(server)
            if peer is not None:
                peer.close()
            Path("/tmp/ipi_" + address).unlink(missing_ok=True)

    assert (tmp_path / "origin").read_text() == "MainThread"
    with open(tmp_path / "RESTART") as restart:
        assert Simulation.load_from_xml(restart, read_only=True).step == 0


if __name__ == "__main__":
    run_simulation(*sys.argv[1:])

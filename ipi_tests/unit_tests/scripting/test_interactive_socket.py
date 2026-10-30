"""Runs an InteractiveSimulation with an external socket driver.

The driver (drivers/py/driver.py, the same script installed as
i-pi-py_driver) is launched as a separate process after the
InteractiveSimulation is created, which opens the socket, and before
the first call to run(), which waits for the client to connect.
"""

import os
import subprocess
import sys
import threading
import time
from pathlib import Path

import numpy as np
import pytest

from ipi.scripting import (
    InteractiveSimulation,
    forcefield_xml,
    motion_nvt_xml,
    simulation_xml,
)
from ipi.utils.softexit import softexit

_ROOT = Path(__file__).resolve().parents[3]
_DRIVER_PATH = _ROOT / "drivers" / "py" / "driver.py"

# generous upper bound: a few steps of a harmonic H2 take well under a second
_RUN_TIMEOUT = 60.0


def _abort(sim, thread, message):
    """Releases the i-PI threads stuck waiting for forces, then fails.

    Force requests only give up once softexit is triggered; cleanup() does
    that without calling sys.exit(). Otherwise the thread pool workers,
    which are joined at interpreter exit, would keep pytest from exiting.
    The softexit state is only reset once the stuck threads have noticed
    the trigger, so that it cannot leak into later tests.
    """

    softexit.cleanup(status="bad", message=message)
    thread.join(timeout=10.0)
    sim.stop()
    softexit.reset()
    pytest.fail(message)


def _run_with_deadline(sim, steps, driver, timeout=_RUN_TIMEOUT):
    """Calls sim.run(steps), failing instead of hanging if the driver never
    connects or dies. run() blocks until the forces are returned, so it is
    executed in a separate thread that can be abandoned on failure."""

    errors = []

    def target():
        try:
            sim.run(steps)
        except BaseException as e:
            errors.append(e)

    thread = threading.Thread(target=target, daemon=True)
    thread.start()
    deadline = time.monotonic() + timeout
    while thread.is_alive():
        if driver.poll() is not None:
            _abort(
                sim, thread, f"driver exited with code {driver.returncode} during run"
            )
        if time.monotonic() > deadline:
            _abort(sim, thread, f"sim.run({steps}) did not complete within {timeout} s")
        thread.join(0.1)
    if errors:
        raise errors[0]


def test_interactive_simulation_with_unix_socket_driver(tmp_path, monkeypatch):
    ase = pytest.importorskip("ase")
    monkeypatch.chdir(tmp_path)  # keeps the simulation outputs out of the way

    address = f"ipiut_interactive_{os.getpid()}"
    atoms = ase.Atoms(
        "H2",
        positions=[[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]],
        cell=5.0 * np.eye(3),
        pbc=True,
    )
    xml = simulation_xml(
        atoms,
        forcefield_xml(name="harmonic", mode="unix", address=address),
        motion_nvt_xml(timestep=0.5 * ase.units.fs),
        temperature=300,
    )

    sim = InteractiveSimulation(xml)

    # makes sure the driver imports this source tree
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        [str(_ROOT)] + [p for p in [env.get("PYTHONPATH")] if p]
    )
    driver = subprocess.Popen(
        [sys.executable, str(_DRIVER_PATH)]
        + ["-u", "-a", address, "-m", "harmonic", "-o", "1e-3"],
        env=env,
    )
    try:
        _run_with_deadline(sim, 5, driver)
        assert sim.step == 5
        assert np.isfinite(sim.properties("potential"))
        _run_with_deadline(sim, 5, driver)
        assert sim.step == 10

        # closing the socket sends the exit message to the driver
        sim.stop()
        assert driver.wait(timeout=30) == 0
    finally:
        if driver.poll() is None:
            driver.kill()

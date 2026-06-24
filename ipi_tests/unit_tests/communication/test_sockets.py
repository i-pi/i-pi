"""In-process transport tests for the i-PI socket server.

Stands up a real FFSocket server and connects the in-tree Python driver over
each transport -- plain unix, inet (TCP) and shared memory -- checking that a
harmonic force round-trips correctly. The driver runs as a forked subprocess,
as it would in practice (its own shared-memory resource tracker). Also covers
request batching and the non-consolidated (per-request worker thread) path.
"""

import importlib.util
import multiprocessing as mp
import os
import socket as pysocket
import time
from pathlib import Path

import numpy as np
import pytest

from ipi.engine.atoms import Atoms
from ipi.engine.cell import Cell
from ipi.engine.forcefields import FFSocket, ForceField
from ipi.interfaces.sockets import InterfaceSocket
from ipi.pes.harmonic import Harmonic_driver

# load the real Python driver client (drivers/py/driver.py) directly by path
_DRIVER_PATH = Path(__file__).resolve().parents[3] / "drivers" / "py" / "driver.py"
_spec = importlib.util.spec_from_file_location("ipi_py_driver", _DRIVER_PATH)
_py_driver = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_py_driver)
run_driver = _py_driver.run_driver

# the driver runs in its own process via fork (own shm resource tracker, like a
# real driver); spawn would need a picklable target
fork_only = pytest.mark.skipif(
    "fork" not in mp.get_all_start_methods(),
    reason="needs the fork start method to launch the driver subprocess",
)

K = 2.0  # harmonic force constant; force = -K*pos


def _free_port():
    s = pysocket.socket()
    s.bind(("localhost", 0))
    port = s.getsockname()[1]
    s.close()
    return port


def _client(mode, address, port):
    """Driver subprocess entry point."""
    run_driver(
        unix=(mode in ("unix", "shm")),
        address=address,
        port=port,
        driver=Harmonic_driver(K),
        shm=(mode == "shm"),
    )


def _exchange_one_force(mode, consolidate=True, batch_size=1):
    """Round-trips a single force evaluation through the chosen transport."""
    address = (
        "localhost"
        if mode == "inet"
        else "ipiut_{}_{}_{}_{}".format(mode, os.getpid(), batch_size, int(consolidate))
    )
    port = _free_port() if mode == "inet" else 0

    iface = InterfaceSocket(
        address=address,
        port=port,
        mode=mode,
        consolidate_messages=consolidate,
        batch_size=batch_size,
    )
    ff = FFSocket(latency=1e-4, name="ut", interface=iface, dopbc=False)

    # open (bind+listen) and fork the driver before the poll thread exists, so
    # the fork happens single-threaded; the client connection waits in the
    # listen backlog until the poll thread starts and accepts it
    iface.open()
    ctx = mp.get_context("fork")
    client = ctx.Process(target=_client, args=(mode, address, port), daemon=True)
    client.start()
    ForceField.start(ff)
    try:
        atoms = Atoms(2)
        atoms.q[:] = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
        cell = Cell(np.eye(3) * 10.0)
        req = ff.queue(atoms, cell, reqid=0)

        deadline = time.time() + 30
        while req["status"] != "Done":
            if time.time() > deadline:
                raise AssertionError("{}: no force returned within 30 s".format(mode))
            time.sleep(0.01)

        pot, force, vir, extra = req["result"]
        return float(pot), np.array(force), atoms.q.copy()
    finally:
        ff.stop()
        client.join(timeout=10)
        if client.is_alive():
            client.terminate()


@fork_only
@pytest.mark.parametrize("mode", ["unix", "inet", "shm"])
def test_transport_roundtrip(mode):
    pot, force, q = _exchange_one_force(mode)
    assert np.allclose(force, -K * q)
    assert np.isclose(pot, 0.5 * K * (q**2).sum())


@fork_only
def test_unix_batched():
    # batch_size>1 exercises the batched dispatch/receive path
    pot, force, q = _exchange_one_force("unix", batch_size=4)
    assert np.allclose(force, -K * q)


@fork_only
def test_shm_batched():
    # batch_size>1 over shared memory exercises the SHM batched payload path
    pot, force, q = _exchange_one_force("shm", batch_size=4)
    assert np.allclose(force, -K * q)


@fork_only
def test_unix_threaded_dispatch():
    # consolidate_messages=False exercises the per-request worker-thread path
    pot, force, q = _exchange_one_force("unix", consolidate=False)
    assert np.allclose(force, -K * q)

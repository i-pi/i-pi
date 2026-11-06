"""In-process transport tests for the i-PI socket server.

Stands up a real FFSocket server and connects the in-tree Python driver over
each transport -- plain unix, inet (TCP) and shared memory -- checking that a
harmonic force round-trips correctly. The driver runs as a forked subprocess,
as it would in practice (its own shared-memory resource tracker). Also covers
request batching, and the parallel evaluation of requests by several clients.
"""

import functools
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
from ipi.utils import scheduler

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


class SlowHarmonic(Harmonic_driver):
    """A harmonic potential that takes `delay` seconds per call."""

    def __init__(self, k, delay):
        super().__init__(k)
        self.delay = delay

    def __call__(self, cell, pos):
        time.sleep(self.delay)
        return super().__call__(cell, pos)


class DyingHarmonic(Harmonic_driver):
    """A harmonic potential whose process exits, without a reply, at the call
    number `ncalls`, after setting the event `died`."""

    def __init__(self, k, ncalls, died):
        super().__init__(k)
        self.ncalls = ncalls
        self.died = died

    def __call__(self, cell, pos):
        self.ncalls -= 1
        if self.ncalls == 0:
            self.died.set()
            os._exit(0)
        return super().__call__(cell, pos)


def _client(mode, address, port, driver):
    """Driver subprocess entry point."""
    run_driver(
        unix=(mode in ("unix", "shm")),
        address=address,
        port=port,
        driver=driver,
        shm=(mode == "shm"),
    )


def _wait_done(requests, timeout=30):
    """Waits until the requests are done, polling the forcefields."""

    deadline = time.time() + timeout

    def done():
        if time.time() > deadline:
            raise AssertionError("no force returned within %d s" % timeout)
        return all(req["status"] == "Done" for req in requests)

    scheduler.wait(done)


def _exchange_forces(
    mode,
    batch_size=1,
    match_mode="auto",
    nclients=1,
    nrequests=None,
    delay=0.0,
    clients=None,
):
    """Evaluates `nrequests` structures (one more than a batch by default)
    through the chosen transport, with `nclients` clients, and returns the
    potentials, the forces, the positions and the time the evaluation took.
    `clients` replaces the default clients with functions of the mode, the
    address and the port, each run in a process."""
    address = (
        "localhost"
        if mode == "inet"
        else "ipiut_{}_{}_{}_{}".format(mode, os.getpid(), batch_size, nclients)
    )
    port = _free_port() if mode == "inet" else 0
    if nrequests is None:
        nrequests = batch_size + 1

    iface = InterfaceSocket(
        address=address,
        port=port,
        mode=mode,
        batch_size=batch_size,
        match_mode=match_mode,
    )
    ff = FFSocket(latency=1e-4, name="ut", interface=iface, dopbc=False)

    # the clients wait in the listen queue until the forcefield is polled
    iface.open()
    ctx = mp.get_context("fork")
    if clients is None:
        driver = SlowHarmonic(K, delay) if delay > 0 else Harmonic_driver(K)
        clients = [functools.partial(_client, driver=driver)] * nclients
    clients = [
        ctx.Process(target=f, args=(mode, address, port), daemon=True) for f in clients
    ]
    for client in clients:
        client.start()
    ForceField.start(ff)
    try:
        atoms = Atoms(2)
        atoms.q_[:] = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
        cell = Cell(np.eye(3) * 10.0)
        positions = [atoms.q.copy() + i for i in range(nrequests)]
        requests = []
        start = time.time()
        for index, position in enumerate(positions):
            atoms.q_[:] = position
            requests.append(ff.queue(atoms, cell, reqid=index))
        _wait_done(requests)
        elapsed = time.time() - start

        potentials = np.array([req["result"][0] for req in requests])
        forces = np.array([req["result"][1] for req in requests])
        for req in requests:
            ff.release(req)
        assert not ff.requests
        return potentials, forces, np.array(positions), elapsed
    finally:
        ff.stop()
        for client in clients:
            client.join(timeout=10)
            if client.is_alive():
                client.terminate()
                client.join(timeout=5)


@fork_only
@pytest.mark.parametrize("mode", ["unix", "inet", "shm"])
@pytest.mark.parametrize("match_mode", ["auto", "any"])
def test_transport_roundtrip(mode, match_mode):
    pot, force, q, _ = _exchange_forces(mode, match_mode=match_mode)
    assert np.allclose(force, -K * q)
    np.testing.assert_allclose(pot, 0.5 * K * (q**2).sum(axis=1))


@fork_only
def test_unix_batched():
    # batch_size>1 exercises the batched dispatch/receive path
    pot, force, q, _ = _exchange_forces("unix", batch_size=4)
    assert np.allclose(force, -K * q)
    np.testing.assert_allclose(pot, 0.5 * K * (q**2).sum(axis=1))


@fork_only
def test_shm_batched():
    # batch_size>1 over shared memory exercises the SHM batched payload path
    pot, force, q, _ = _exchange_forces("shm", batch_size=4)
    assert np.allclose(force, -K * q)
    np.testing.assert_allclose(pot, 0.5 * K * (q**2).sum(axis=1))


@fork_only
@pytest.mark.parametrize("mode", ["unix", "shm"])
@pytest.mark.parametrize("match_mode", ["auto", "any"])
def test_clients_compute_in_parallel(mode, match_mode):
    # 8 requests over 4 clients that take 0.3 s each: about 0.6 s, against
    # 2.4 s if the clients worked one after the other
    pot, force, q, elapsed = _exchange_forces(
        mode, match_mode=match_mode, nclients=4, nrequests=8, delay=0.3
    )
    assert np.allclose(force, -K * q)
    np.testing.assert_allclose(pot, 0.5 * K * (q**2).sum(axis=1))
    assert elapsed < 1.5


@fork_only
@pytest.mark.parametrize("mode", ["unix", "inet", "shm"])
def test_client_replaced_mid_run(mode):
    """A client dies while it computes a request, and another one connects
    some time later: i-PI waits for it, and the request is completed."""
    died = mp.get_context("fork").Event()

    def dying(mode, address, port):
        _client(mode, address, port, DyingHarmonic(K, 2, died))

    def late(mode, address, port):
        # i-PI is left without clients until this one connects
        died.wait(30)
        time.sleep(0.2)
        _client(mode, address, port, Harmonic_driver(K))

    pot, force, q, _ = _exchange_forces(mode, nrequests=3, clients=[dying, late])
    assert died.is_set()
    assert np.allclose(force, -K * q)
    np.testing.assert_allclose(pot, 0.5 * K * (q**2).sum(axis=1))

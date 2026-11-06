"""i-PI socket protocol: wire format, malformed replies and recovery from failed clients.

The peer side of each connection is a plain socket that writes protocol bytes,
so these tests only rely on what a driver could observe: the bytes i-PI sends,
the errors it raises, and the state of the force requests.
"""

import contextlib
import gc
import os
import select
import socket
import struct
import tempfile
import threading
import time
import warnings
from multiprocessing import resource_tracker, shared_memory

import numpy as np
import pytest

from ipi.engine.forcefields import ForceRequest
from ipi.interfaces import sockets
from ipi.interfaces.sockets import Disconnected, Driver, InvalidSize, Status

# Protocol headers are fixed 12-byte ASCII fields, spelled out here rather
# than taken from sockets.MESSAGE so that a typo there is caught.
MESSAGE = {
    name: name.upper().encode().ljust(12)
    for name in (
        "status",
        "ready",
        "needinit",
        "havedata",
        "init",
        "posdata",
        "getforce",
        "forceready",
        "exit",
    )
}
STATUS, READY, HAVEDATA = MESSAGE["status"], MESSAGE["ready"], MESSAGE["havedata"]
DRIVERS = pytest.mark.parametrize(
    "connection", [Driver, sockets.SHMDriver], indirect=True, ids=["socket", "shm"]
)


def payload(natoms=1, header=None, extra=b"", energy=1.25, forces=None):
    """FORCEREADY reply with a unit virial and forces 0, 1, 2, ..."""
    forces = np.arange(3 * natoms, dtype=np.float64) if forces is None else forces
    return b"".join(
        [
            MESSAGE["forceready"] if header is None else header,
            np.float64(energy).tobytes(),
            np.int32(natoms).tobytes(),
            np.asarray(forces, dtype=np.float64).tobytes(),
            np.eye(3, dtype=np.float64).tobytes(),
            np.int32(len(extra)).tobytes(),
            extra,
        ]
    )


def request(rid=17, shift=0.0):
    """Two-atom request in which only the second atom is active."""
    return ForceRequest(
        id=rid,
        pars="",
        pos=np.arange(6, dtype=float) + shift,
        active=np.array([3, 4, 5]),
        cell=(np.eye(3), np.eye(3)),
        result=None,
        status="Queued",
        start=-1,
    )


def pending(peer):
    """Every byte the server has sent to `peer` so far."""
    peer.setblocking(False)
    data = b""
    try:
        while chunk := peer.recv(65536):
            data += chunk
    except BlockingIOError:
        pass
    finally:
        peer.setblocking(True)
    return data


@contextlib.contextmanager
def no_leaked_descriptors():
    """Fails if the block leaves descriptors open, or leaves them to the garbage collector."""
    # the tracker of the shared memory segments keeps a pipe open once started
    resource_tracker.ensure_running()
    gc.collect()
    before = len(os.listdir("/dev/fd"))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ResourceWarning)
        yield
        gc.collect()
    unclosed = [str(w.message) for w in caught if w.category is ResourceWarning]
    assert unclosed == []
    assert len(os.listdir("/dev/fd")) == before, "descriptors left open"


def announced_segments(wire):
    """Names of the shared-memory segments announced by a POSDATA trigger, and
    the length of the trigger."""
    names, offset = [], 16
    for _ in range(6):  # positions, cell, inverse cell, energy, forces, virial
        length = int(np.frombuffer(wire[offset : offset + 4], dtype=np.int32)[0])
        names.append(wire[offset + 4 : offset + 4 + length].decode())
        offset += 4 + length
    return names, offset


def assert_pending(r):
    assert r["status"] != "Done"
    assert r["result"] is None


def assert_completed(r, energy=1.0, forces=(0, 0, 0, 0, 1, 2)):
    assert r["status"] == "Done"
    assert r["result"][0] == pytest.approx(energy)
    np.testing.assert_array_equal(r["result"][1], forces)


@pytest.fixture
def connection(request):
    """A driver on one end of a socket pair, and the raw peer on the other."""
    server, peer = socket.socketpair()
    server.settimeout(0.05)
    driver = getattr(request, "param", Driver)(server)
    server.close()
    try:
        yield driver, peer
    finally:
        driver.disconnect()  # also releases shared-memory segments
        peer.close()


@pytest.fixture
def server(request):
    """A listening interface (TCP unless another mode is requested) and a
    function that connects handshaken clients."""
    mode = getattr(request, "param", "inet")
    directory = tempfile.TemporaryDirectory(prefix="ipi-socket-", dir="/tmp")
    address = dict(address="127.0.0.1", port=0)
    if mode != "inet":
        address = dict(address="server", sockets_prefix=directory.name + "/")
    interface = sockets.InterfaceSocket(mode=mode, timeout=0.3, **address)
    interface.open()
    interface.requests = []
    interface.offset = 0.25
    peers = []

    def connect(reply=READY):
        peer = socket.socket(interface.server.family, socket.SOCK_STREAM)
        peer.settimeout(2)
        peer.connect(interface.server.getsockname())
        peers.append(peer)
        interface.poll()  # accepts the client, and asks for its status
        assert peer.recv(12) == STATUS
        peer.sendall(reply)
        wait_readable(interface.clients[-1])
        interface.poll()  # reads the status
        return peer

    try:
        yield interface, connect
    finally:
        interface.close()
        for peer in peers:
            peer.close()
        directory.cleanup()


def wait_readable(client, timeout=2.0):
    """Waits until the reply of `client` can be read."""
    select.select([client], [], [], timeout)


def run_pool(interface, requests, timeout=2.0):
    """Polls until every request is done or the deadline passes."""
    deadline = time.time() + timeout
    while time.time() < deadline and not all(r["status"] == "Done" for r in requests):
        interface.poll()
        time.sleep(0.001)
    interface.poll()


# -- one client ---------------------------------------------------------------


@pytest.mark.parametrize(
    "reply,status",
    [
        (READY, Status.Up | Status.Ready),
        (MESSAGE["needinit"], Status.Up | Status.NeedsInit),
        (HAVEDATA, Status.Up | Status.HasData),
    ],
)
def test_status_replies(connection, reply, status):
    driver, peer = connection
    driver.ask_status()
    assert pending(peer) == STATUS
    assert driver.expect == "status" and not driver.free()
    peer.sendall(reply)
    driver.read_status()
    assert driver.status == status
    assert driver.expect is None
    assert driver.free() == (reply != HAVEDATA)


@pytest.mark.parametrize("reply", [b"UNKNOWN     ", b"BAD", b""])
def test_bad_status_reply_disconnects(connection, monkeypatch, reply):
    monkeypatch.setattr(sockets, "NTIMEOUT", 1)
    driver, peer = connection
    driver.ask_status()
    peer.sendall(reply)
    if reply != b"UNKNOWN     ":
        peer.close()
    with pytest.raises(Disconnected):
        driver.read_status()


@pytest.mark.parametrize("batch_size,pars", [(1, ""), (3, ""), (3, "model=harmonic")])
def test_initialization_is_sent_with_the_positions(connection, batch_size, pars):
    driver, peer = connection
    driver.status = Status.Up | Status.NeedsInit
    driver.batch_size = batch_size
    r = request()
    r["pars"] = pars
    driver.send_job([r])
    wire = pending(peer)
    assert wire[:12] == MESSAGE["init"]
    rid, length = np.frombuffer(wire[12:20], dtype=np.int32).tolist()
    expected = pars
    if batch_size > 1:
        expected += ("," if pars else "") + " batch_size:3"
    assert rid == 17 and wire[20 : 20 + length].decode() == expected
    assert wire[20 + length : 32 + length] == MESSAGE["posdata"]
    assert wire.endswith(MESSAGE["getforce"])
    assert driver.expect == "forces"


@pytest.mark.parametrize("extra", [b"", b'{"dipole":[1,2,3]}', b"x" * 2048])
def test_force_reply_is_read_completely(connection, extra):
    driver, peer = connection
    driver.status = Status.Up | Status.Ready
    driver.send_job([request()])
    assert pending(peer) == b"".join(
        [
            MESSAGE["posdata"],
            np.eye(3).tobytes(),
            np.eye(3).tobytes(),
            np.int32(1).tobytes(),
            np.arange(3.0, 6.0).tobytes(),
            MESSAGE["getforce"],
        ]
    )
    peer.sendall(payload(extra=extra))
    ((energy, force, virial, extras),) = driver.read_job()
    assert energy == 1.25
    np.testing.assert_array_equal(force, [0, 1, 2])
    np.testing.assert_array_equal(virial, np.eye(3))
    assert extras.get("raw", "") == extra.decode()
    assert driver.lastreq == 17 and driver.job is None and driver.expect is None


def test_batch_is_padded_and_results_keep_request_order(connection):
    driver, peer = connection
    driver.batch_size = 3
    driver.status = Status.Up | Status.Ready
    requests = [request(), request(rid=18, shift=10)]
    driver.send_job(requests)
    wire = pending(peer)
    assert wire[:12] == MESSAGE["posdata"]
    assert np.frombuffer(wire[12:16], dtype=np.int32)[0] == 1
    cells = np.frombuffer(wire[16:448], dtype=np.float64).reshape(3, 2, 3, 3)
    np.testing.assert_array_equal(cells, np.tile(np.eye(3), (3, 2, 1, 1)))
    positions = np.frombuffer(wire[448:-12], dtype=np.float64).reshape(3, 3)
    np.testing.assert_array_equal(positions, [[3, 4, 5], [13, 14, 15], [13, 14, 15]])
    assert wire[-12:] == MESSAGE["getforce"]

    # The padded third result is read from the wire but assigned to no request.
    peer.sendall(
        MESSAGE["forceready"]
        + np.array([1.0, 2.0, 999.0]).tobytes()
        + np.int32(1).tobytes()
        + np.arange(9, dtype=float).tobytes()
        + np.tile(np.eye(3), (3, 1, 1)).tobytes()
        + b"".join(
            np.int32(len(extra)).tobytes() + extra
            for extra in [b"first", b"second", b"padding"]
        )
    )
    results = driver.read_job()
    assert len(results) == 2
    for index, (energy, force, _, extras) in enumerate(results):
        assert energy == index + 1
        np.testing.assert_array_equal(force, [3 * index, 3 * index + 1, 3 * index + 2])
        assert extras["raw"] == ["first", "second"][index]


@pytest.mark.parametrize("batched", [False, True])
@pytest.mark.parametrize(
    "reply,error",
    [
        (payload(natoms=2), InvalidSize),
        (payload(header=b"BADHEADER   "), Disconnected),
        (payload()[:20], Disconnected),
        (payload(extra=b"abc")[:-2], Disconnected),
        (b"", Disconnected),
    ],
    ids=["wrong-size", "bad-header", "truncated", "truncated-extra", "silent"],
)
def test_malformed_force_reply_raises(connection, monkeypatch, batched, reply, error):
    monkeypatch.setattr(sockets, "NTIMEOUT", 1)
    driver, peer = connection
    driver.status = Status.Up | Status.Ready
    driver.batch_size = 2 if batched else 1
    driver.send_job([request()])
    pending(peer)
    if batched:
        reply = reply.replace(np.float64(1.25).tobytes(), np.zeros(2).tobytes(), 1)
    peer.sendall(reply)
    peer.shutdown(socket.SHUT_WR)
    with pytest.raises(error):
        driver.read_job()


@DRIVERS
def test_send_to_vanished_client_fails(connection):
    driver, peer = connection
    driver.status = Status.Up | Status.Ready
    peer.close()
    with pytest.raises(OSError):
        driver.send_job([request()])
    assert driver.job is None


def test_exit_on_disconnect_requests_soft_exit(connection):
    driver, peer = connection
    sockets.softexit.reset()
    try:
        driver.exit_on_disconnect = True
        driver.disconnect()
        assert sockets.softexit.requested == "Client shutdown."
        assert not sockets.softexit.triggered
    finally:
        sockets.softexit.reset()
    assert peer.recv(12) == MESSAGE["exit"]
    assert driver.status == Status.Disconnected
    assert driver.fileno() == -1


@pytest.mark.parametrize("connection", [sockets.SHMDriver], indirect=True)
def test_shared_memory_segments_are_announced_and_released(connection):
    driver, peer = connection
    driver.status = Status.Up | Status.Ready
    driver.send_job([request()])
    wire = pending(peer)
    assert wire[:12] == MESSAGE["posdata"]
    assert np.frombuffer(wire[12:16], dtype=np.int32)[0] == 1
    names, offset = announced_segments(wire)
    assert wire[offset:] == MESSAGE["getforce"]
    driver.disconnect()
    for name in names:
        with pytest.raises(FileNotFoundError):
            shared_memory.SharedMemory(name=name)


def test_shm_clients_use_distinct_segments():
    """Clients connected one after the other share no segment, also when the
    accepted sockets reuse the same descriptor numbers"""
    drivers, peers, names = [], [], []
    try:
        for _ in range(3):
            server, peer = socket.socketpair()
            server.settimeout(0.05)
            drivers.append(sockets.SHMDriver(server))
            server.close()
            peers.append(peer)
        for driver, peer in zip(drivers, peers):
            driver.status = Status.Up | Status.Ready
            driver.send_job([request()])
            names += announced_segments(pending(peer))[0]
        assert len(set(names)) == len(names)
    finally:
        for driver in drivers:
            driver.disconnect()
        for peer in peers:
            peer.close()


def test_unreferenced_shm_driver_closes_its_descriptor():
    server, peer = socket.socketpair()
    driver = sockets.SHMDriver(server)
    server.close()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ResourceWarning)
        del driver
        gc.collect()
    assert select.select([peer], [], [], 1)[0], "the driver descriptor is still open"
    assert peer.recv(12) == b""
    peer.close()


# -- the pool of clients -------------------------------------------------------

SERVER_MODES = pytest.mark.parametrize("server", ["unix", "inet", "shm"], indirect=True)


@pytest.mark.parametrize("reply", [READY, MESSAGE["needinit"]])
def test_handshake_and_exit_on_close(server, reply):
    interface, connect = server
    peer = connect(reply)
    assert len(interface.clients) == 1
    assert interface.clients[0].status == Status.Up | (
        Status.Ready if reply == READY else Status.NeedsInit
    )
    interface.close()
    assert peer.recv(12) == MESSAGE["exit"]
    assert interface.clients == []


@SERVER_MODES
@pytest.mark.parametrize("closed", [False, True], ids=["half-closed", "closed"])
def test_disconnected_handshake_is_not_accepted(server, closed):
    interface, _ = server
    with no_leaked_descriptors():
        peer = socket.socket(interface.server.family, socket.SOCK_STREAM)
        peer.connect(interface.server.getsockname())
        if closed:
            # e.g. a check that the server is listening
            peer.close()
        else:
            peer.shutdown(socket.SHUT_WR)
        deadline = time.time() + 2
        interface.poll()
        while interface.clients and time.time() < deadline:
            interface.poll()
        assert not interface.clients
        peer.close()


@SERVER_MODES
def test_dead_client_is_closed_when_removed(server):
    interface, connect = server
    with no_leaked_descriptors():
        peer = connect()
        driver = interface.clients[0]
        interface.requests.append(request())
        interface.poll()  # sends the request
        # a reset makes the exit message fail also over TCP
        peer.setsockopt(socket.SOL_SOCKET, socket.SO_LINGER, struct.pack("ii", 1, 0))
        peer.close()
        wait_readable(driver)
        interface.poll()
        assert interface.clients == []
        assert driver.fileno() == -1
        assert interface.requests[0]["status"] == "Queued"


def test_request_roundtrip_sends_active_atoms_and_expands_forces(server):
    interface, connect = server
    peer = connect()
    first, second = request(), request(rid=18)
    interface.requests[:] = [first]
    interface.poll()
    assert pending(peer) == b"".join(
        [
            MESSAGE["posdata"],
            np.eye(3).tobytes(),
            np.eye(3).tobytes(),
            np.int32(1).tobytes(),
            np.arange(3.0, 6.0).tobytes(),
            MESSAGE["getforce"],
        ]
    )
    assert first["status"] == "Running"
    peer.sendall(payload(extra=b'{"energy_component":2}'))
    run_pool(interface, [first])
    assert_completed(first)
    np.testing.assert_array_equal(first["result"][2], np.eye(3))
    assert first["result"][3]["energy_component"] == 2
    # the status of the client is asked after the results
    assert pending(peer) == STATUS
    # a later reply must not overwrite the results already handed out
    peer.sendall(READY)
    interface.requests[:] = [first, second]
    run_pool(interface, [second], timeout=0.2)
    peer.sendall(payload(energy=2.25, forces=[6.0, 7.0, 8.0]))
    run_pool(interface, [second])
    assert_completed(second, energy=2.0, forces=[0, 0, 0, 6, 7, 8])
    assert_completed(first)


def test_poll_does_not_wait_for_clients(server):
    """Requests are sent to the clients that become free while others still
    compute, and a poll returns at once, whatever the clients do."""
    interface, connect = server
    interface.timeout = 0.0
    slow, fast = connect(), connect()
    requests = [request(rid=i) for i in range(3)]
    interface.requests[:] = requests
    start = time.time()
    interface.poll()
    assert [r["status"] for r in requests] == ["Running", "Running", "Queued"]
    for peer in (slow, fast):
        assert MESSAGE["posdata"] in pending(peer)
    fast.sendall(payload() + READY)
    deadline = time.time() + 2
    while requests[2]["status"] == "Queued" and time.time() < deadline:
        interface.poll()
    assert requests[2]["status"] == "Running"
    assert MESSAGE["posdata"] in pending(fast)
    assert sum(r["status"] == "Done" for r in requests) == 1
    assert time.time() - start < 1.0


def test_pool_sends_each_request_to_one_client(server):
    interface, connect = server
    r = request()
    peers = [connect(), connect()]
    interface.requests[:] = [r]
    interface.poll()
    received = [MESSAGE["posdata"] in pending(peer) for peer in peers]
    assert sum(received) == 1
    peers[received.index(True)].sendall(payload() + READY)
    run_pool(interface, [r])
    assert_completed(r)


# the clients replaced over shared memory are tested with real drivers in test_sockets.py
@pytest.mark.parametrize("server", ["unix", "inet"], indirect=True)
@pytest.mark.parametrize("failure", ["in-flight", "silent"])
def test_request_of_failed_client_is_completed_by_another(server, failure):
    interface, connect = server
    failing = connect()
    r = request()
    interface.requests[:] = [r]
    interface.poll()
    assert r["status"] == "Running"
    if failure == "in-flight":
        failing.shutdown(socket.SHUT_WR)
    # a silent client is dropped after the time-out
    deadline = time.time() + 2
    while interface.clients and time.time() < deadline:
        interface.poll()
    assert interface.clients == []
    assert (r["status"], r["start"]) == ("Queued", -1)
    assert_pending(r)
    replacement = connect()
    replacement.sendall(payload() + READY)
    run_pool(interface, [r])
    assert_completed(r)


def test_wrong_size_reply_requeues_request(server):
    """A client that returns the forces of a different number of atoms is
    removed, and its request queued again"""
    interface, connect = server
    interface.timeout = 0.0
    peer = connect()
    r = request()
    interface.requests[:] = [r]
    interface.poll()
    peer.sendall(payload(natoms=2) + READY)
    deadline = time.time() + 2
    while interface.clients and time.time() < deadline:
        interface.poll()
    assert interface.clients == []
    assert (r["status"], r["start"]) == ("Queued", -1)
    assert_pending(r)


@pytest.mark.parametrize("reply", [b"UNKNOWN     ", HAVEDATA])
def test_client_in_unexpected_status_is_removed(server, reply):
    interface, connect = server
    r = request()
    interface.requests[:] = [r]
    peer = connect(reply)
    assert interface.clients == []
    assert r["status"] == "Queued"
    assert MESSAGE["posdata"] not in pending(peer)


def test_completed_request_of_failed_client_stays_done(server):
    """A request completed by a client that then disconnects (e.g. it exits
    after sending the forces) is not queued again"""
    interface, connect = server
    peer = connect()
    r = request()
    interface.requests[:] = [r]
    interface.poll()
    peer.sendall(payload())
    peer.close()
    deadline = time.time() + 2
    while interface.clients and time.time() < deadline:
        interface.poll()
    assert interface.clients == []
    assert_completed(r)


def test_locked_client_keeps_replica_with_large_id(server):
    interface, connect = server
    interface.match_mode = "lock"
    peer = connect()
    for _ in range(2):
        # Each round uses a new, equal integer id, as after a checkpoint restart.
        r = request(rid=int("1000"))
        interface.requests[:] = [r]
        interface.poll()
        assert MESSAGE["posdata"] in pending(peer)
        peer.sendall(payload() + READY)
        run_pool(interface, [r], timeout=0.5)
        assert_completed(r)


def test_locked_request_of_failed_client_goes_to_another(server):
    """The id bound to a client that fails is sent to another client in the
    poll that drops it."""
    interface, connect = server
    interface.match_mode = "lock"
    failing = connect()
    first = request()
    interface.requests[:] = [first]
    interface.poll()
    failing.sendall(payload() + READY)
    run_pool(interface, [first])
    assert_completed(first)
    replacement = connect()
    second = request()
    interface.requests[:] = [second]
    interface.poll()
    assert MESSAGE["posdata"] in pending(failing)
    assert MESSAGE["posdata"] not in pending(replacement)
    driver = interface.clients[0]
    failing.shutdown(socket.SHUT_WR)
    wait_readable(driver)
    interface.poll()
    assert interface.clients[0] is not driver
    assert second["status"] == "Running"
    assert MESSAGE["posdata"] in pending(replacement)
    replacement.sendall(payload() + READY)
    run_pool(interface, [second])
    assert_completed(second)


def test_close_waits_for_the_status_reply(server):
    """A client that answers the status request sent after its results is
    not interrupted while it writes the reply, which would kill it."""
    interface, connect = server
    peer = connect()
    r = request()
    interface.requests[:] = [r]
    interface.poll()
    peer.sendall(payload())
    run_pool(interface, [r])
    assert_completed(r)
    assert pending(peer).endswith(STATUS)
    late_reply = threading.Timer(0.05, peer.sendall, [READY])
    start = time.time()
    late_reply.start()
    interface.close()
    # the exit is sent once the reply has been read
    assert time.time() - start >= 0.05
    late_reply.join()
    assert peer.recv(12) == MESSAGE["exit"]


def test_closed_interface_stops_the_requests(server):
    interface, connect = server
    connect()
    r = request()
    interface.requests[:] = [r]
    interface.poll()
    assert r["status"] == "Running"
    interface.close()
    assert interface.fds() == []
    assert interface.poll()
    assert r["status"] == "Exit"

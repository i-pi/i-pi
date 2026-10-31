"""i-PI socket protocol: wire format, malformed replies and recovery from failed clients.

The peer side of each connection is a plain socket that writes protocol bytes,
so these tests only rely on what a driver could observe: the bytes i-PI sends,
the errors it raises, and the state of the force requests.
"""

import socket
import tempfile
import threading
import time
from concurrent.futures import Future
from multiprocessing import shared_memory

import numpy as np
import pytest

from ipi.engine.forcefields import ForceRequest
from ipi.interfaces import sockets
from ipi.interfaces.sockets import (
    Disconnected,
    Driver,
    InvalidSize,
    InvalidStatus,
    Status,
)

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
        offset=0.25,
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


def assert_pending(r):
    assert r["status"] != "Done"
    assert r["result"] is None
    assert not r._event_done.is_set()


def assert_completed(r, energy=1.0, forces=(0, 0, 0, 0, 1, 2)):
    assert r["status"] == "Done"
    assert r._event_done.is_set()
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
        try:
            driver.shutdown()  # also releases shared-memory segments
        except OSError:
            pass
        driver.close()
        peer.close()


@pytest.fixture
def server(monkeypatch):
    """A listening TCP interface and a function that connects handshaken clients."""
    monkeypatch.setattr(sockets, "SERVERTIMEOUT", 0.01)
    interface = sockets.InterfaceSocket(
        mode="inet", address="127.0.0.1", port=0, timeout=0.3
    )
    interface.open()
    interface.requests = []
    interface.offset = 0.25
    peers = []

    def connect(reply=READY):
        peer = socket.create_connection(interface.server.getsockname(), timeout=2)
        peers.append(peer)
        peer.sendall(reply)
        interface.pool_update()
        assert peer.recv(12) == STATUS
        return peer

    try:
        yield interface, connect
    finally:
        interface.close()
        interface.server.close()
        interface.executor.shutdown(wait=True)
        for peer in peers:
            peer.close()


def run_pool(interface, requests, timeout=2.0):
    """Distributes work until every request is done or the deadline passes."""
    deadline = time.time() + timeout
    while time.time() < deadline and not all(r._event_done.is_set() for r in requests):
        interface.pool_distribute()
    interface.pool_distribute()


@pytest.mark.parametrize(
    "reply,status",
    [
        (READY, Status.Up | Status.Ready),
        (MESSAGE["needinit"], Status.Up | Status.NeedsInit),
        (HAVEDATA, Status.Up | Status.HasData),
        (b"UNKNOWN     ", Status.Up),
        (b"BAD", Status.Disconnected),
    ],
)
def test_status_replies(connection, reply, status):
    driver, peer = connection
    peer.sendall(reply)
    assert driver.get_status() == status
    assert pending(peer) == STATUS


def test_status_wait_sends_a_single_request(connection, monkeypatch):
    """A slow client is asked once, not repeatedly, while i-PI waits."""
    monkeypatch.setattr(sockets, "SELECTTIMEOUT", 0.01)
    driver, peer = connection
    late_reply = threading.Timer(0.2, peer.sendall, [READY])
    late_reply.start()
    try:
        assert driver.get_status() == Status.Up | Status.Ready
    finally:
        late_reply.join()
    assert pending(peer) == STATUS


def test_closed_peer_is_disconnected(connection):
    driver, peer = connection
    peer.close()
    assert driver.get_status() == Status.Disconnected


@pytest.mark.parametrize("batch_size,pars", [(1, ""), (3, ""), (3, "model=harmonic")])
def test_initialization_wire_format(connection, batch_size, pars):
    driver, peer = connection
    peer.sendall(MESSAGE["needinit"])
    driver.get_status()
    driver.batch_size = batch_size
    driver.initialize(17, pars)
    wire = pending(peer)
    assert wire[:12] == STATUS
    init = wire[12:]
    assert init[:12] == MESSAGE["init"]
    assert np.frombuffer(init[12:20], dtype=np.int32).tolist() == [17, len(init) - 20]
    expected = pars
    if batch_size > 1:
        expected += ("," if pars else "") + " batch_size:3"
    assert init[20:].decode() == expected


@pytest.mark.parametrize("extra", [b"", b'{"dipole":[1,2,3]}', b"x" * 2048])
def test_getforce_reads_complete_reply(connection, extra):
    driver, peer = connection
    peer.sendall(HAVEDATA + payload(extra=extra))
    driver.get_status()
    energy, force, virial, extras = driver.getforce()
    assert pending(peer) == STATUS + MESSAGE["getforce"]
    assert energy == 1.25
    np.testing.assert_array_equal(force, [0, 1, 2])
    np.testing.assert_array_equal(virial, np.eye(3))
    assert extras.get("raw", "") == extra.decode()


def test_getforce_requires_data(connection):
    driver, _ = connection
    with pytest.raises(InvalidStatus):
        driver.getforce()


@DRIVERS
@pytest.mark.parametrize("method", ["initialize", "sendpos"])
def test_invalid_status_reports_protocol_error(connection, method):
    driver, peer = connection
    peer.sendall(HAVEDATA)
    driver.get_status()
    args = (17, "") if method == "initialize" else (np.zeros(3), (np.eye(3), np.eye(3)))
    with pytest.raises(InvalidStatus):
        getattr(driver, method)(*args)


def test_request_roundtrip_sends_active_atoms_and_expands_forces(connection):
    driver, peer = connection
    first, second = request(), request(rid=18)
    peer.sendall(READY)
    driver.get_status()
    assert driver.dispatch_send(first)
    assert pending(peer) == b"".join(
        [
            STATUS,
            MESSAGE["posdata"],
            np.eye(3).tobytes(),
            np.eye(3).tobytes(),
            np.int32(1).tobytes(),
            np.arange(3.0, 6.0).tobytes(),
            MESSAGE["getforce"],
        ]
    )
    peer.sendall(payload(extra=b'{"energy_component":2}') + READY)
    assert driver.dispatch_recv(first)
    assert_completed(first)
    np.testing.assert_array_equal(first["result"][2], np.eye(3))
    assert first["result"][3]["energy_component"] == 2
    # A later reply must not overwrite the results already handed out.
    assert driver.dispatch_send(second)
    peer.sendall(payload(energy=2.25, forces=[6.0, 7.0, 8.0]) + READY)
    assert driver.dispatch_recv(second)
    assert_completed(second, energy=2.0, forces=[0, 0, 0, 6, 7, 8])
    assert_completed(first)


def test_batch_is_padded_and_results_keep_request_order(connection):
    driver, peer = connection
    driver.batch_size = 3
    requests = [request(), request(rid=18, shift=10)]
    peer.sendall(READY)
    driver.get_status()
    assert driver.dispatch_send_batch(requests)
    wire = pending(peer)[12:]
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
        + MESSAGE["needinit"]
    )
    assert driver.dispatch_recv_batch(requests)
    for index, r in enumerate(requests):
        assert_completed(
            r,
            energy=index + 0.75,
            forces=[0, 0, 0, 3 * index, 3 * index + 1, 3 * index + 2],
        )
        assert r["result"][3]["raw"] == ["first", "second"][index]
    assert driver.status == Status.Up | Status.NeedsInit


@pytest.mark.parametrize("batched", [False, True])
@pytest.mark.parametrize(
    "reply,closed,error",
    [
        (payload(natoms=2), False, InvalidSize),
        (payload(header=b"BADHEADER   "), False, None),
        (payload()[:20], True, None),
        (payload(extra=b"abc")[:-2], True, None),
        (b"", False, None),
    ],
    ids=["wrong-size", "bad-header", "truncated", "truncated-extra", "silent"],
)
def test_malformed_force_reply_never_completes_request(
    connection, monkeypatch, batched, reply, closed, error
):
    monkeypatch.setattr(sockets, "NTIMEOUT", 1)
    driver, peer = connection
    r = request()
    peer.sendall(READY)
    driver.get_status()
    assert driver.dispatch_send_batch([r]) if batched else driver.dispatch_send(r)
    pending(peer)
    peer.sendall(reply)
    if closed:
        peer.shutdown(socket.SHUT_WR)
    receive = (
        (lambda: driver.dispatch_recv_batch([r]))
        if batched
        else (lambda: driver.dispatch_recv(r))
    )
    if error is None:
        assert receive() is False
        assert driver.status == Status.Disconnected
    else:
        with pytest.raises(error):
            receive()
    assert_pending(r)


@DRIVERS
@pytest.mark.parametrize("batched", [False, True])
@pytest.mark.parametrize("state", ["closed", "has-data"])
def test_unavailable_client_receives_no_coordinates(connection, batched, state):
    driver, peer = connection
    r = request()
    if state == "closed":
        peer.close()
    else:
        peer.sendall(HAVEDATA)
    assert not (driver.dispatch_send_batch([r]) if batched else driver.dispatch_send(r))
    assert_pending(r)
    if state == "closed":
        assert driver.status == Status.Disconnected
    else:
        assert pending(peer) == STATUS


@DRIVERS
@pytest.mark.parametrize("batched", [False, True])
def test_send_to_vanished_client_fails(connection, batched):
    driver, peer = connection
    r = request()
    peer.sendall(READY)
    driver.get_status()
    peer.close()
    assert not (driver.dispatch_send_batch([r]) if batched else driver.dispatch_send(r))
    assert driver.status == Status.Disconnected
    assert_pending(r)


@pytest.mark.parametrize("needs_initialization", [False, True])
def test_worker_dispatch_roundtrip(connection, needs_initialization):
    driver, peer = connection
    r = request()
    replies = READY
    if needs_initialization:
        replies = MESSAGE["needinit"] + replies
    peer.sendall(replies + HAVEDATA + payload(extra=b'{"energy_component":2}') + READY)
    driver.dispatch(r)
    assert_completed(r)
    np.testing.assert_array_equal(r["result"][2], np.eye(3))
    assert r["result"][3]["energy_component"] == 2
    assert driver.status == Status.Up | Status.Ready
    wire = pending(peer)
    assert (MESSAGE["init"] in wire) == needs_initialization
    assert MESSAGE["posdata"] in wire
    assert MESSAGE["getforce"] in wire


@pytest.mark.parametrize(
    "stage", ["closed", "not-ready", "no-data", "truncated", "wrong-size"]
)
def test_worker_dispatch_failure_leaves_request_pending(connection, monkeypatch, stage):
    monkeypatch.setattr(sockets, "NTIMEOUT", 1)
    driver, peer = connection
    r = request()
    if stage == "closed":
        peer.close()
    elif stage == "not-ready":
        peer.sendall(HAVEDATA)
    elif stage == "no-data":
        peer.sendall(READY + READY)
    elif stage == "truncated":
        peer.sendall(READY + HAVEDATA + payload()[:30])
        peer.shutdown(socket.SHUT_WR)
    else:
        r["active"] = np.arange(6)
        peer.sendall(READY + HAVEDATA + payload(natoms=1) + READY)
    if stage == "wrong-size":
        with pytest.raises(InvalidSize):
            driver.dispatch(r)
    else:
        driver.dispatch(r)
    assert r["status"] != "Done"
    assert not r._event_done.is_set()
    if stage == "truncated":
        assert driver.status == Status.Disconnected


def test_exit_on_disconnect_requests_soft_exit(connection, monkeypatch):
    driver, peer = connection
    messages = []
    called = threading.Event()

    def record_exit(**kwargs):
        messages.append(kwargs["message"])
        called.set()

    monkeypatch.setattr(sockets.softexit, "trigger", record_exit)
    driver.exit_on_disconnect = True
    driver.shutdown()
    assert called.wait(2)
    driver.exit_on_disconnect = False
    assert messages == ["Client shutdown."]
    assert peer.recv(12) == MESSAGE["exit"]
    assert driver.status == Status.Disconnected


@pytest.mark.parametrize("connection", [sockets.SHMDriver], indirect=True)
def test_shared_memory_segments_are_announced_and_released(connection):
    driver, peer = connection
    peer.sendall(READY)
    driver.get_status()
    assert driver.dispatch_send(request())
    wire = pending(peer)[12:]
    assert wire[:12] == MESSAGE["posdata"]
    assert np.frombuffer(wire[12:16], dtype=np.int32)[0] == 1
    names, offset = [], 16
    for _ in range(6):  # positions, cell, inverse cell, energy, forces, virial
        length = int(np.frombuffer(wire[offset : offset + 4], dtype=np.int32)[0])
        names.append(wire[offset + 4 : offset + 4 + length].decode())
        offset += 4 + length
    assert wire[offset:] == MESSAGE["getforce"]
    driver.shutdown()
    for name in names:
        with pytest.raises(FileNotFoundError):
            shared_memory.SharedMemory(name=name)


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


@pytest.mark.parametrize("closed", [False, True], ids=["half-closed", "closed"])
def test_disconnected_handshake_is_not_accepted(server, closed):
    interface, _ = server
    with socket.create_connection(interface.server.getsockname(), timeout=1) as peer:
        if closed:
            # e.g. a check that the server is listening
            peer.close()
        else:
            peer.shutdown(socket.SHUT_WR)
        interface.pool_update()
        assert not interface.clients


@pytest.mark.parametrize("consolidated", [True, False], ids=["select", "workers"])
def test_pool_sends_each_request_to_one_client(server, consolidated):
    interface, connect = server
    interface.consolidate_messages = consolidated
    r = request()
    interface.requests = [r]
    replies = payload() + READY
    if not consolidated:
        replies = READY + HAVEDATA + replies
    peers = [connect(), connect()]
    for peer in peers:
        peer.sendall(replies)
    run_pool(interface, [r])
    assert_completed(r)
    assert sum(MESSAGE["posdata"] in pending(peer) for peer in peers) == 1


@pytest.mark.parametrize("failure", ["in-flight", "silent"])
def test_request_of_failed_client_is_completed_by_another(server, failure):
    interface, connect = server
    r = request()
    interface.requests = [r]
    failing = connect()
    if failure == "in-flight":
        failing.shutdown(socket.SHUT_WR)
    interface.pool_distribute()  # a silent client is dropped after the timeout
    assert r["status"] == "Queued"
    assert_pending(r)
    replacement = connect()  # pool_update also removes the failed client
    assert len(interface.clients) == 1
    replacement.sendall(payload() + READY)
    run_pool(interface, [r])
    assert_completed(r)


def test_client_with_pending_data_gets_no_work(server):
    interface, connect = server
    r = request()
    interface.requests = [r]
    peer = connect(HAVEDATA)
    interface.pool_distribute()
    assert r["status"] == "Queued"
    assert_pending(r)
    assert pending(peer) == b""


def test_locked_client_keeps_replica_with_large_id(server):
    interface, connect = server
    interface.match_mode = "lock"
    peer = connect()
    for _ in range(2):
        # Each round uses a new, equal integer id, as after a checkpoint restart.
        r = request(rid=int("1000"))
        interface.requests = [r]
        peer.sendall(payload() + READY)
        run_pool(interface, [r], timeout=0.5)
        assert_completed(r)


def test_worker_timeout_does_not_lose_request(server, connection, monkeypatch):

    interface, _ = server
    interface.consolidate_messages = False
    driver, _ = connection
    driver.status = Status.Up | Status.Busy
    r = request()
    r["status"], r["start"] = "Running", 1
    future = Future()
    future.set_result(None)
    interface.requests = [r]
    interface.clients = [driver]
    interface.jobs = [[r, driver, future]]
    monkeypatch.setattr(sockets.time, "time", lambda: 3)
    interface.pool_distribute()
    interface.pool_update()
    assert r["status"] == "Queued"
    assert_pending(r)


def test_force_header_eof_is_disconnected(connection, monkeypatch):
    driver, peer = connection
    peer.sendall(HAVEDATA)
    driver.get_status()
    peer.shutdown(socket.SHUT_WR)
    receive = driver.recv_msg
    calls = []

    def bounded_receive(length=12):
        # Stop a receive loop deterministically instead of hanging the test.
        calls.append(1)
        if len(calls) > 1:
            raise ConnectionResetError("test guard: receive retried after EOF")
        return receive(length)

    monkeypatch.setattr(driver, "recv_msg", bounded_receive)
    with pytest.raises(Disconnected):
        driver.getforce()
    assert len(calls) == 1


def test_server_descriptor_closes_even_if_shutdown_fails(server, monkeypatch):
    interface, _ = server

    def fail_shutdown(*args):
        raise OSError("injected shutdown failure")

    monkeypatch.setattr(socket.socket, "shutdown", fail_shutdown)
    interface.close()
    assert interface.server.fileno() == -1


def test_duplicate_unix_address_does_not_replace_listener():
    with tempfile.TemporaryDirectory(prefix="ipi-socket-", dir="/tmp") as directory:
        first = sockets.InterfaceSocket(
            address="server", sockets_prefix=directory + "/"
        )
        second = sockets.InterfaceSocket(
            address="server", sockets_prefix=directory + "/"
        )
        first.open()
        try:
            with pytest.raises(RuntimeError, match="Error opening unix socket"):
                second.open()
            with socket.socket(socket.AF_UNIX, socket.SOCK_STREAM) as peer:
                peer.settimeout(1)
                peer.connect(directory + "/server")
                accepted, _ = first.server.accept()
                accepted.close()
        finally:
            second.server.close()
            first.close()
            first.server.close()
            first.executor.shutdown(wait=True)

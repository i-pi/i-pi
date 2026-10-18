"""Real socket pairs exercise malformed replies and interrupted force transfers."""

import socket

import numpy as np
import pytest
from ipi.interfaces import sockets
from ipi.interfaces.sockets import (
    Driver,
    Status,
    Disconnected,
    InvalidSize,
    InvalidStatus,
)

# Protocol headers are fixed 12-byte ASCII fields. Keep the peer's vocabulary
# independent of the production MESSAGE dictionary so spelling errors are caught.
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


@pytest.fixture
def connection(request):
    server, client = socket.socketpair()
    server.settimeout(0.05)
    driver = getattr(request, "param", Driver)(server)
    server.close()
    try:
        yield driver, client
    finally:
        if isinstance(driver, sockets.SHMDriver):
            driver._free_shm()
        driver.close()
        client.close()


@pytest.mark.parametrize("method", ["_getstatus_direct", "_getstatus_select"])
@pytest.mark.parametrize(
    "reply,status",
    [
        (MESSAGE["ready"], Status.Up | Status.Ready),
        (MESSAGE["needinit"], Status.Up | Status.NeedsInit),
        (MESSAGE["havedata"], Status.Up | Status.HasData),
        (b"UNKNOWN     ", Status.Up),
        (b"BAD", Status.Disconnected),
    ],
)
def test_status_reply_decoding(connection, method, reply, status):
    driver, client = connection
    client.sendall(reply)
    assert getattr(driver, method)() == status
    assert client.recv(12) == MESSAGE["status"]


def test_direct_status_timeout_preserves_pending_request(connection):
    driver, client = connection
    status = driver._getstatus_direct()
    assert status & Status.Timeout
    assert driver.waitstatus
    assert client.recv(12) == MESSAGE["status"]
    client.sendall(MESSAGE["ready"])
    assert driver._getstatus_direct() == Status.Up | Status.Ready
    assert not driver.waitstatus


def test_interrupted_array_transfer_is_disconnected(connection):
    driver, client = connection
    client.sendall(np.array([1.0], dtype=np.float64).tobytes())
    client.shutdown(socket.SHUT_WR)
    with pytest.raises(Disconnected):
        driver.recvall(np.zeros(3))


def test_unanswered_array_transfer_has_bounded_retries(connection, monkeypatch):
    driver, _ = connection
    monkeypatch.setattr(sockets, "NTIMEOUT", 1)
    with pytest.raises(Disconnected):
        driver.recvall(np.zeros(3))


def payload(natoms=1, header=None, extra=b""):
    return b"".join(
        [
            MESSAGE["forceready"] if header is None else header,
            np.float64(1.25).tobytes(),
            np.int32(natoms).tobytes(),
            np.arange(3, dtype=np.float64).tobytes(),
            np.eye(3, dtype=np.float64).tobytes(),
            np.int32(len(extra)).tobytes(),
            extra,
        ]
    )


@pytest.mark.parametrize("extra", [b"", b'{"dipole":[1,2,3]}'])
def test_force_payload_has_independent_storage(connection, extra):
    driver, client = connection
    client.sendall(payload(extra=extra))
    result = driver._recv_forces_bulk(1)
    assert result[0] == 1.25
    np.testing.assert_array_equal(result[1], np.arange(3))
    np.testing.assert_array_equal(result[2], np.eye(3))
    driver._buf[:] = 0
    np.testing.assert_array_equal(result[1], np.arange(3))
    np.testing.assert_array_equal(result[2], np.eye(3))


@pytest.mark.parametrize(
    "data,error",
    [(payload(natoms=2), InvalidSize), (payload(header=b"BADHEADER   "), Disconnected)],
)
@pytest.mark.parametrize("batched", [False, True])
def test_malformed_force_payload_is_rejected(connection, data, error, batched):
    driver, client = connection
    client.sendall(data)
    with pytest.raises(error):
        if batched:
            driver._recv_forces_bulk_batch(1, 1)
        else:
            driver._recv_forces_bulk(1)


def test_getforce_rejects_client_without_data(connection):
    driver, _ = connection
    with pytest.raises(InvalidStatus):
        driver.getforce()
    with pytest.raises(ValueError, match="empty"):
        driver.dispatch_send_batch([])


def request_for_active_atom():
    from ipi.engine.forcefields import ForceRequest

    return ForceRequest(
        {
            "id": 17,
            "pars": "",
            "pos": np.arange(6, dtype=float),
            "active": np.array([3, 4, 5]),
            "cell": (np.eye(3), np.eye(3)),
            "offset": 0.25,
            "result": None,
            "status": "Running",
        }
    )


@pytest.mark.parametrize("batched", [False, True])
@pytest.mark.parametrize(
    "connection", [Driver, sockets.SHMDriver], indirect=True, ids=["socket", "shm"]
)
def test_dispatch_to_disconnected_peer_leaves_request_unfinished(connection, batched):
    driver, client = connection
    request = request_for_active_atom()
    driver.status = Status.Up | Status.Ready
    client.close()
    success = (
        driver.dispatch_send_batch([request])
        if batched
        else driver.dispatch_send(request)
    )
    assert success is False
    assert driver.status == Status.Disconnected
    assert request["result"] is None
    assert not request._event_done.is_set()


@pytest.mark.parametrize("batched", [False, True])
@pytest.mark.parametrize(
    "connection", [Driver, sockets.SHMDriver], indirect=True, ids=["socket", "shm"]
)
def test_partial_force_reply_cannot_complete_request(connection, batched):
    driver, client = connection
    request = request_for_active_atom()
    if isinstance(driver, sockets.SHMDriver):
        driver._alloc_shm(1)
        # Shared memory carries the numerical result; extras still use the socket.
        client.sendall(MESSAGE["forceready"] + np.int32(3).tobytes() + b"x")
    else:
        client.sendall(MESSAGE["forceready"] + np.float64(1.0).tobytes())
    client.shutdown(socket.SHUT_WR)
    success = (
        driver.dispatch_recv_batch([request])
        if batched
        else driver.dispatch_recv(request)
    )
    assert success is False
    assert driver.status == Status.Disconnected
    assert request["result"] is None
    assert not request._event_done.is_set()


@pytest.mark.parametrize("batched", [False, True])
def test_active_atom_force_is_expanded_and_offset_applied(connection, batched):
    driver, client = connection
    request = request_for_active_atom()
    client.sendall(payload() + MESSAGE["ready"])
    success = (
        driver.dispatch_recv_batch([request])
        if batched
        else driver.dispatch_recv(request)
    )
    assert success
    assert request["status"] == "Done"
    assert request._event_done.is_set()
    assert request["result"][0] == pytest.approx(1.0)
    np.testing.assert_array_equal(request["result"][1], [0, 0, 0, 0, 1, 2])
    assert driver.lastreq == 17


@pytest.mark.parametrize("batched", [False, True])
@pytest.mark.parametrize(
    "connection", [Driver, sockets.SHMDriver], indirect=True, ids=["socket", "shm"]
)
def test_send_timeout_does_not_report_success(connection, monkeypatch, batched):
    driver, _ = connection
    request = request_for_active_atom()
    driver.status = Status.Up | Status.Ready
    original_timeout = sockets.TIMEOUT
    monkeypatch.setattr(sockets, "TIMEOUT", original_timeout)

    def timeout(data):
        raise socket.timeout("injected send timeout")

    monkeypatch.setattr(driver, "sendall", timeout)
    success = (
        driver.dispatch_send_batch([request])
        if batched
        else driver.dispatch_send(request)
    )
    assert success is False
    assert driver.status == Status.Timeout
    assert sockets.TIMEOUT == 2 * original_timeout
    assert not request._event_done.is_set()


@pytest.mark.parametrize("batch_size,pars", [(1, ""), (3, ""), (3, "model=harmonic")])
def test_initialization_wire_format(connection, batch_size, pars):
    driver, client = connection
    driver.batch_size = batch_size
    driver.status = Status.Up | Status.NeedsInit
    driver.initialize(17, pars)
    data = client.recv(1024)
    assert data[:12] == MESSAGE["init"]
    assert np.frombuffer(data[12:20], dtype=np.int32).tolist() == [17, len(data) - 20]
    expected = pars
    if batch_size > 1:
        expected += ("," if pars else "") + " batch_size:3"
    assert data[20:].decode() == expected


@pytest.mark.parametrize("method", ["initialize", "sendpos"])
@pytest.mark.xfail(
    strict=True,
    raises=TypeError,
    reason="InvalidStatus message concatenates an integer status with a string",
)
@pytest.mark.parametrize(
    "connection", [Driver, sockets.SHMDriver], indirect=True, ids=["socket", "shm"]
)
def test_invalid_status_reports_protocol_error(connection, method):
    driver, _ = connection
    driver.status = Status.Up | Status.Busy
    args = (17, "") if method == "initialize" else (np.zeros(3), (np.eye(3), np.eye(3)))
    with pytest.raises(InvalidStatus):
        getattr(driver, method)(*args)


@pytest.mark.parametrize("extra", [b"", b'{"dipole":[1,2,3]}', b"x" * 2048])
def test_legacy_force_receive_preserves_extras(connection, extra):
    driver, client = connection
    driver.status = Status.Up | Status.HasData
    client.sendall(payload(extra=extra))
    energy, force, virial, extras = driver.getforce()
    assert client.recv(12) == MESSAGE["getforce"]
    assert energy == 1.25
    np.testing.assert_array_equal(force, [0, 1, 2])
    np.testing.assert_array_equal(virial, np.eye(3))
    assert extras.get("raw", "") == extra.decode()


@pytest.mark.parametrize("part", ["prefix", "extras"])
@pytest.mark.parametrize("closed", [False, True])
def test_bulk_receive_incomplete_message_is_bounded(
    connection, monkeypatch, part, closed
):
    driver, client = connection
    monkeypatch.setattr(sockets, "NTIMEOUT", 1)
    if part == "extras":
        client.sendall(payload(extra=b"abc")[:-2])
    if closed:
        client.shutdown(socket.SHUT_WR)
    with pytest.raises(Disconnected):
        driver._recv_forces_bulk(1)


def test_batch_padding_and_result_order(connection):
    driver, client = connection
    driver.batch_size = 3
    driver.status = Status.Up | Status.Ready
    requests = [request_for_active_atom() for _ in range(2)]
    requests[1]["id"] = 18
    requests[1]["pos"] += 10
    assert driver.dispatch_send_batch(requests)
    wire = client.recv(4096)
    assert wire[:12] == MESSAGE["posdata"]
    assert np.frombuffer(wire[12:16], dtype=np.int32)[0] == 1
    cells = np.frombuffer(wire[16:448], dtype=np.float64).reshape(3, 2, 3, 3)
    np.testing.assert_array_equal(cells, np.tile(np.eye(3), (3, 2, 1, 1)))
    positions = np.frombuffer(wire[448:-12], dtype=np.float64).reshape(3, 3)
    np.testing.assert_array_equal(positions, [[3, 4, 5], [13, 14, 15], [13, 14, 15]])
    assert wire[-12:] == MESSAGE["getforce"]

    # The padded third result must be consumed but never assigned to a request.
    client.sendall(
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
    for index, request in enumerate(requests):
        assert request["result"][0] == index + 0.75
        np.testing.assert_array_equal(
            request["result"][1], [0, 0, 0, 3 * index, 3 * index + 1, 3 * index + 2]
        )
        assert request["result"][3]["raw"] == ["first", "second"][index]
        assert request._event_done.is_set()
    assert driver.lastreq == 18
    assert driver.status == Status.Up | Status.NeedsInit


@pytest.mark.parametrize("batched", [False, True])
def test_failed_dispatch_requeues_requests(connection, batched):
    driver, client = connection
    interface = sockets.InterfaceSocket(batch_size=3 if batched else 1)
    requests = [request_for_active_atom() for _ in range(2 if batched else 1)]
    interface.prlist = requests.copy()
    interface.clients = [driver]
    interface.jobs = []
    driver.status = Status.Up | Status.Ready
    client.close()
    assert not interface.dispatch_free_client(driver)
    assert not interface.jobs
    assert driver.status == Status.Disconnected
    assert interface.poll_iter == sockets.UPDATEFREQ
    for request in requests:
        assert request["status"] == "Queued"
        assert request["start"] == -1
        assert request["result"] is None
        assert not request._event_done.is_set()


@pytest.mark.parametrize("batched", [False, True])
@pytest.mark.parametrize("failure", ["disconnect", "timeout"])
def test_inflight_failure_requeues_and_replacement_completes(
    connection, monkeypatch, batched, failure
):
    driver, client = connection
    interface = sockets.InterfaceSocket(batch_size=2 if batched else 1, timeout=1)
    request = request_for_active_atom()
    request["start"] = 1
    job = [request] if batched else request
    driver.batch_size = interface.batch_size
    driver.status = Status.Up | Status.HasData
    interface.requests = [request]
    interface.prlist = []
    interface.clients = [driver]
    interface.jobs = [[job, driver, None]]
    if failure == "disconnect":
        client.shutdown(socket.SHUT_WR)
    else:
        # Advance the clock past the request deadline, without a timing race.
        monkeypatch.setattr(sockets.time, "time", lambda: 3)
    interface._pool_distribute_select()
    assert not interface.jobs
    assert request["status"] == "Queued"
    assert request["result"] is None
    assert not request._event_done.is_set()
    assert driver.status == Status.Disconnected

    # A replacement connection can return the original request's result.
    server, replacement = socket.socketpair()
    replacement.settimeout(1)
    recovered = Driver(server)
    server.close()
    try:
        replacement.sendall(payload() + MESSAGE["ready"])
        assert recovered.dispatch_recv(request)
        assert request["status"] == "Done"
        assert request._event_done.is_set()
        assert request["result"][0] == 1.0
        np.testing.assert_array_equal(request["result"][1], [0, 0, 0, 0, 1, 2])
    finally:
        recovered.close()
        replacement.close()


@pytest.mark.parametrize(
    "status", [Status.Disconnected, Status.Up | Status.HasData, Status.Up]
)
def test_busy_or_disconnected_client_is_not_assigned_work(connection, status):
    driver, _ = connection
    driver.status = status
    interface = sockets.InterfaceSocket()
    request = request_for_active_atom()
    request["status"] = "Queued"
    interface.prlist = [request]
    assert not interface.dispatch_free_client(driver)
    assert interface.prlist == [request]
    assert request["status"] == "Queued"


@pytest.mark.parametrize(
    "lastreq,match,locked,assigned",
    [
        (17, "match", False, True),
        (18, "match", False, False),
        (17, "none", False, False),
        (None, "none", False, True),
        (17, "free", True, False),
        (17, "free", False, True),
    ],
)
def test_request_matching_respects_client_assignment(
    connection, lastreq, match, locked, assigned
):
    driver, client = connection
    driver.status = Status.Up | Status.Ready
    driver.lastreq = lastreq
    driver.locked = locked
    interface = sockets.InterfaceSocket()
    request = request_for_active_atom()
    request["status"] = "Queued"
    interface.prlist = [request]
    interface.clients = [driver]
    interface.jobs = []
    assert interface.dispatch_free_client(driver, match_ids=match) == assigned
    assert request["status"] == ("Running" if assigned else "Queued")
    assert len(interface.jobs) == int(assigned)
    assert len(interface.prlist) == int(not assigned)
    if assigned:
        assert client.recv(4096).startswith(MESSAGE["posdata"])


def test_bulk_large_extra_does_not_overwrite_force_result(connection):
    driver, client = connection
    extra = b"x" * 2048
    client.sendall(payload(extra=extra))
    energy, force, virial, extras = driver._recv_forces_bulk(1)
    assert energy == 1.25
    np.testing.assert_array_equal(force, [0, 1, 2])
    np.testing.assert_array_equal(virial, np.eye(3))
    assert extras["raw"] == extra.decode()


@pytest.mark.parametrize("connection", [sockets.SHMDriver], indirect=True)
@pytest.mark.parametrize("batched", [False, True])
def test_shared_memory_rejects_invalid_force_header(connection, batched):
    driver, client = connection
    client.sendall(b"BADHEADER   ")
    with pytest.raises(Disconnected):
        if batched:
            driver._recv_forces_bulk_batch(1, 1)
        else:
            driver._recv_forces_bulk(1)


@pytest.fixture
def listening_interface(monkeypatch):
    # Use an OS-assigned TCP port; no port-selection race or socket file to leave behind.
    monkeypatch.setattr(sockets, "SERVERTIMEOUT", 0.01)
    interface = sockets.InterfaceSocket(mode="inet", address="127.0.0.1", port=0)
    interface.open()
    interface.requests = []
    try:
        yield interface
    finally:
        interface.close()
        interface.server.close()  # also release the descriptor if shutdown failed
        interface.executor.shutdown(wait=True)


@pytest.mark.parametrize("reply", [MESSAGE["ready"], MESSAGE["needinit"]])
def test_new_client_handshake_and_shutdown(listening_interface, reply):
    interface = listening_interface
    with socket.create_connection(interface.server.getsockname(), timeout=1) as client:
        client.sendall(reply)
        interface.pool_update()
        assert client.recv(12) == MESSAGE["status"]
        assert len(interface.clients) == 1
        assert interface.clients[0].status == Status.Up | (
            Status.Ready if reply == MESSAGE["ready"] else Status.NeedsInit
        )
        interface.close()
        assert client.recv(12) == MESSAGE["exit"]
        assert interface.clients == interface.jobs == []


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="Handshake uses bitwise OR with Up, accepting a disconnected client",
)
def test_disconnected_handshake_is_not_accepted(listening_interface):
    interface = listening_interface
    with socket.create_connection(interface.server.getsockname(), timeout=1) as client:
        client.shutdown(socket.SHUT_WR)
        interface.pool_update()
        assert not interface.clients


def test_disconnected_worker_job_is_requeued(listening_interface, connection):
    from concurrent.futures import Future

    interface = listening_interface
    driver, _ = connection
    request = request_for_active_atom()
    request["start"] = 1
    future = Future()
    future.set_result(None)
    driver.status = Status.Disconnected
    interface.clients = [driver]
    interface.jobs = [[request, driver, future]]
    interface.pool_update()
    assert interface.clients == interface.jobs == []
    assert request["status"] == "Queued"
    assert request["start"] == -1
    assert not request._event_done.is_set()


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="Worker-mode timeout removes the job before pool_update can requeue its request",
)
def test_worker_timeout_does_not_lose_request(
    listening_interface, connection, monkeypatch
):
    from concurrent.futures import Future

    interface = listening_interface
    interface.consolidate_messages = False
    driver, _ = connection
    driver.status = Status.Up | Status.Busy
    request = request_for_active_atom()
    request["start"] = 1
    future = Future()
    future.set_result(None)
    interface.requests = [request]
    interface.clients = [driver]
    interface.jobs = [[request, driver, future]]
    monkeypatch.setattr(sockets.time, "time", lambda: 3)
    interface.pool_distribute()
    interface.pool_update()
    assert request["status"] == "Queued"
    assert request["result"] is None
    assert not request._event_done.is_set()


def test_select_status_retry_keeps_one_pending_request(connection, monkeypatch):
    driver, client = connection
    monkeypatch.setattr(sockets, "SELECTTIMEOUT", 0.01)
    assert driver._getstatus_select() == Status.Busy
    assert driver.waitstatus
    assert client.recv(12) == MESSAGE["status"]
    client.sendall(MESSAGE["ready"])
    assert driver._getstatus_select() == Status.Up | Status.Ready
    assert not driver.waitstatus
    client.setblocking(False)
    with pytest.raises(BlockingIOError):
        client.recv(12)


@pytest.mark.parametrize("method", ["_getstatus_select", "_getstatus_direct"])
def test_closed_connection_status_is_disconnected(connection, method):
    driver, client = connection
    client.close()
    assert getattr(driver, method)() == Status.Disconnected


@pytest.mark.parametrize(
    "connection", [Driver, sockets.SHMDriver], indirect=True, ids=["socket", "shm"]
)
@pytest.mark.parametrize("failure", [socket.timeout, BrokenPipeError])
def test_legacy_position_send_failure(connection, monkeypatch, failure):
    driver, _ = connection
    driver.status = Status.Up | Status.Ready
    timeout = sockets.TIMEOUT
    monkeypatch.setattr(sockets, "TIMEOUT", timeout)

    def fail_send(data):
        raise failure("injected send failure")

    monkeypatch.setattr(driver, "sendall", fail_send)
    if failure is socket.timeout:
        driver.sendpos(np.zeros(3), (np.eye(3), np.eye(3)))
        assert driver.status == Status.Timeout
        assert sockets.TIMEOUT == 2 * timeout
    else:
        with pytest.raises(BrokenPipeError):
            driver.sendpos(np.zeros(3), (np.eye(3), np.eye(3)))
        assert sockets.TIMEOUT == timeout


@pytest.mark.parametrize("needs_initialization", [False, True])
def test_worker_dispatch_initializes_and_returns_active_atom_forces(
    connection, needs_initialization
):
    driver, client = connection
    request = request_for_active_atom()
    driver.status = Status.Up
    replies = MESSAGE["ready"]
    if needs_initialization:
        replies = MESSAGE["needinit"] + replies
    client.sendall(
        replies
        + MESSAGE["havedata"]
        + payload(extra=b'{"energy_component":2}')
        + MESSAGE["ready"]
    )
    driver.dispatch(request)
    assert request["status"] == "Done"
    assert request._event_done.is_set()
    assert request["result"][0] == 1.0
    np.testing.assert_array_equal(request["result"][1], [0, 0, 0, 0, 1, 2])
    np.testing.assert_array_equal(request["result"][2], np.eye(3))
    assert request["result"][3]["energy_component"] == 2
    assert driver.lastreq == 17
    assert driver.status == Status.Up | Status.Ready
    wire = client.recv(4096)
    assert (MESSAGE["init"] in wire) == needs_initialization
    assert MESSAGE["posdata"] in wire
    assert MESSAGE["getforce"] in wire


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="Request matching compares integer identity instead of value",
)
def test_matching_large_request_id_uses_value(connection):
    driver, _ = connection
    request = request_for_active_atom()
    request["id"] = int("1000")
    driver.lastreq = int("1000")
    assert driver.lastreq is not request["id"]
    driver.status = Status.Up | Status.Ready
    interface = sockets.InterfaceSocket(match_mode="lock")
    interface.prlist = [request]
    interface.clients = [driver]
    interface.jobs = []
    assert interface.dispatch_free_client(driver, match_ids="match")


@pytest.mark.parametrize("state", [Status.Disconnected, Status.Up | Status.Busy])
@pytest.mark.parametrize("batched", [False, True])
@pytest.mark.parametrize(
    "connection", [Driver, sockets.SHMDriver], indirect=True, ids=["socket", "shm"]
)
def test_unready_driver_cannot_send_coordinates(connection, state, batched):
    driver, client = connection
    driver.status = state
    request = request_for_active_atom()
    if state & Status.Up:
        client.sendall(MESSAGE["havedata"])
    send = driver.dispatch_send_batch if batched else driver.dispatch_send
    assert not send([request] if batched else request)
    assert request["result"] is None
    assert not request._event_done.is_set()
    client.setblocking(False)
    if state & Status.Up:
        assert client.recv(12) == MESSAGE["status"]
    with pytest.raises(BlockingIOError):
        client.recv(4096)


@pytest.mark.parametrize("method", ["_getstatus_direct", "_getstatus_select"])
@pytest.mark.parametrize("error", [socket.timeout, ConnectionResetError])
def test_status_receive_os_failure(connection, monkeypatch, method, error):
    driver, client = connection
    client.sendall(MESSAGE["ready"])  # select reports a readable descriptor

    def fail_receive(length=12):
        raise error("injected receive failure")

    monkeypatch.setattr(driver, "recv_msg", fail_receive)
    status = getattr(driver, method)()
    assert status == (
        Status.Up | Status.Busy | Status.Timeout
        if error is socket.timeout
        else Status.Disconnected
    )
    assert client.recv(12) == MESSAGE["status"]


@pytest.mark.parametrize("first_reply", [socket.timeout, b"UNEXPECTED  "])
def test_force_header_retry_accepts_subsequent_valid_reply(
    connection, monkeypatch, first_reply
):
    driver, client = connection
    client.sendall(MESSAGE["forceready"])
    receive = driver.recv_msg
    calls = []

    def receive_after_retry():
        calls.append(1)
        if len(calls) == 1:
            if first_reply is socket.timeout:
                raise socket.timeout("injected receive timeout")
            return first_reply
        return receive()

    monkeypatch.setattr(driver, "recv_msg", receive_after_retry)
    driver._recv_forceready()
    assert len(calls) == 2


def test_force_header_connection_reset_is_disconnected(connection, monkeypatch):
    driver, _ = connection

    def reset():
        raise ConnectionResetError("injected reset")

    monkeypatch.setattr(driver, "recv_msg", reset)
    with pytest.raises(Disconnected):
        driver._recv_forceready()


def test_status_waits_until_busy_client_replies(connection, monkeypatch):
    driver, client = connection
    client.sendall(MESSAGE["ready"])
    receive = driver._getstatus
    calls = []

    def busy_once():
        calls.append(1)
        return Status.Up | Status.Busy if len(calls) == 1 else receive()

    monkeypatch.setattr(driver, "_getstatus", busy_once)
    assert driver.get_status() == Status.Up | Status.Ready
    assert len(calls) == 2


def test_failed_initialization_refreshes_client_status(connection, monkeypatch):
    driver, client = connection
    driver.status = Status.Up | Status.NeedsInit
    send = driver.sendall

    def fail_init_only(data):
        if data.startswith(MESSAGE["init"]):
            raise BrokenPipeError("injected INIT failure")
        return send(data)

    monkeypatch.setattr(driver, "sendall", fail_init_only)
    client.sendall(MESSAGE["needinit"])
    driver.initialize(17, "")
    assert driver.status == Status.Up | Status.NeedsInit
    assert client.recv(12) == MESSAGE["status"]


@pytest.mark.parametrize(
    "stage", ["disconnected", "not-ready", "no-data", "force-reset", "wrong-size"]
)
def test_worker_dispatch_failure_never_marks_request_done(
    connection, monkeypatch, stage
):
    driver, client = connection
    request = request_for_active_atom()
    driver.status = Status.Disconnected if stage == "disconnected" else Status.Up
    if stage == "not-ready":
        client.sendall(MESSAGE["havedata"])
    elif stage != "disconnected":
        client.sendall(
            MESSAGE["ready"] + MESSAGE["ready" if stage == "no-data" else "havedata"]
        )
    if stage == "force-reset":

        def reset_force():
            raise Disconnected()

        monkeypatch.setattr(driver, "getforce", reset_force)
    if stage == "wrong-size":
        request["active"] = np.arange(6)
        client.sendall(payload())  # one atom returned for a two-atom request
        with pytest.raises(InvalidSize):
            driver.dispatch(request)
    else:
        driver.dispatch(request)
    assert request["status"] != "Done"
    assert not request._event_done.is_set()
    if stage == "force-reset":
        assert driver.status == Status.Disconnected


@pytest.mark.parametrize("consolidated", [False, True])
@pytest.mark.parametrize("assignment", ["different-id", "already-owned"])
def test_locked_matching_leaves_unmatched_requests_queued(
    connection, consolidated, assignment
):
    driver, _ = connection
    driver.status = Status.Up | Status.Ready
    driver.lastreq = 18 if assignment == "different-id" else None
    request = request_for_active_atom()
    request["status"] = "Queued"
    interface = sockets.InterfaceSocket(
        match_mode="lock", consolidate_messages=consolidated
    )
    interface.requests = [request]
    interface.clients = [driver]
    interface.jobs = []
    if assignment == "already-owned":
        interface.prlist = [request]
        from types import SimpleNamespace

        # A second client already owns this replica and has work in flight.
        owner = SimpleNamespace(lastreq=17)
        interface.clients.append(owner)
        assert not interface.dispatch_free_client(driver, match_ids="none")
    else:
        interface.pool_distribute()
    assert request["status"] == "Queued"
    assert not interface.jobs


@pytest.mark.parametrize("consolidated", [False, True])
def test_idle_pool_flags_disconnected_client_for_removal(connection, consolidated):
    driver, _ = connection
    driver.status = Status.Disconnected
    interface = sockets.InterfaceSocket(consolidate_messages=consolidated)
    interface.clients = [driver]
    interface.requests = [request_for_active_atom()]
    interface.jobs = []
    interface.poll_iter = 0
    interface.pool_distribute()
    assert interface.poll_iter == sockets.UPDATEFREQ
    interface.requests = []
    interface.pool_distribute()
    assert not interface.jobs


def test_exit_on_disconnect_requests_soft_exit_from_background_thread(
    connection, monkeypatch
):
    import threading

    driver, client = connection
    completed = threading.Event()
    callers = []

    def record_exit(**kwargs):
        callers.append((threading.current_thread(), kwargs["message"]))
        completed.set()

    monkeypatch.setattr(sockets.softexit, "trigger", record_exit)
    driver.exit_on_disconnect = True
    driver.shutdown()
    assert completed.wait(2)
    callers[0][0].join(timeout=2)
    assert callers[0][0] is not threading.current_thread()
    assert callers[0][1] == "Client shutdown."
    assert client.recv(12) == MESSAGE["exit"]
    assert driver.status == Status.Disconnected


@pytest.mark.parametrize("connection", [sockets.SHMDriver], indirect=True)
def test_shared_memory_reclaims_stale_segment_and_unlinks_resources(connection):
    from multiprocessing import shared_memory

    driver, _ = connection
    name = driver._shm_prefix + "pos"
    stale = shared_memory.SharedMemory(create=True, size=24, name=name)
    try:
        driver._alloc_shm(1)
        names = [segment.name for segment in driver._shm.values()]
        assert driver._shm["pos"].name == name
        # A segment may already have been removed by another process.
        driver._shm["force"].unlink()
        driver._free_shm()
        assert not driver._shm and not driver._view
        for segment_name in names:
            with pytest.raises(FileNotFoundError):
                shared_memory.SharedMemory(name=segment_name)
    finally:
        stale.close()


@pytest.mark.parametrize("connection", [sockets.SHMDriver], indirect=True)
def test_shared_memory_empty_extra_and_empty_batch(connection):
    driver, client = connection
    client.sendall(np.int32(0).tobytes())
    assert driver._recv_extra() == ""
    with pytest.raises(ValueError, match="empty"):
        driver.dispatch_send_batch([])


def test_duplicate_unix_address_does_not_replace_listener():
    import tempfile

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
            with socket.socket(socket.AF_UNIX, socket.SOCK_STREAM) as client:
                client.settimeout(1)
                client.connect(directory + "/server")
                accepted, _ = first.server.accept()
                accepted.close()
        finally:
            second.server.close()
            first.close()
            first.server.close()
            first.close()  # repeated cleanup tolerates an already removed socket file
            first.executor.shutdown(wait=True)


def test_unknown_transport_is_rejected():
    with pytest.raises(NameError, match="not implemented"):
        sockets.InterfaceSocket(mode="unknown").open()


def test_tcp_option_failure_does_not_prevent_connection(monkeypatch):
    original = sockets.socket.socket

    class SocketWithoutOptions(original):
        def setsockopt(self, *args):
            raise OSError("injected unsupported socket option")

    monkeypatch.setattr(sockets.socket, "socket", SocketWithoutOptions)
    interface = sockets.InterfaceSocket(mode="inet", address="127.0.0.1", port=0)
    interface.open()
    try:
        with original() as client:
            client.settimeout(1)
            client.connect(interface.server.getsockname())
            accepted, _ = interface.server.accept()
            accepted.close()
    finally:
        interface.close()
        interface.server.close()  # also release the descriptor if shutdown failed
        interface.executor.shutdown(wait=True)


@pytest.mark.parametrize("consolidated", [False, True])
def test_timeout_survives_shutdown_error(connection, monkeypatch, consolidated):
    driver, _ = connection
    request = request_for_active_atom()
    request["start"] = 1
    driver.status = Status.Up | Status.Busy
    interface = sockets.InterfaceSocket(timeout=1)
    interface.clients = [driver]
    interface.requests = [request]
    interface.jobs = [[request, driver, None]]

    def fail_shutdown(*args):
        raise OSError("injected shutdown failure")

    monkeypatch.setattr(driver, "shutdown", fail_shutdown)
    monkeypatch.setattr(sockets.time, "time", lambda: 3)
    if consolidated:
        interface._pool_distribute_select()
        assert not interface.jobs
        assert request["status"] == "Queued"
    else:
        assert interface.check_job_finished(request, driver, None) == 0
    assert driver.status == Status.Disconnected
    assert driver.fileno() == -1
    assert not request._event_done.is_set()


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="Force-header receive loops on EOF because bytes are compared to an empty string",
)
def test_force_header_eof_is_disconnected(connection, monkeypatch):
    driver, client = connection
    client.shutdown(socket.SHUT_WR)
    receive = driver.recv_msg
    calls = []

    def bounded_receive():
        # The first EOF is real. Stop a broken loop deterministically instead of hanging pytest.
        calls.append(1)
        if len(calls) > 1:
            raise ConnectionResetError("test guard: receive retried after EOF")
        return receive()

    monkeypatch.setattr(driver, "recv_msg", bounded_receive)
    with pytest.raises(Disconnected):
        driver._recv_forceready()

    assert len(calls) == 1


@pytest.mark.parametrize("error", [OSError, ValueError])
def test_select_failure_does_not_complete_or_drop_jobs(connection, monkeypatch, error):
    driver, _ = connection
    driver.status = Status.Up | Status.HasData
    request = request_for_active_atom()
    interface = sockets.InterfaceSocket()
    interface.clients = [driver]
    interface.requests = [request]
    interface.jobs = [[request, driver, None]]

    def fail_select(*args):
        raise error("injected select failure")

    monkeypatch.setattr(sockets.select, "select", fail_select)
    interface._pool_distribute_select()
    assert interface.jobs == [[request, driver, None]]
    assert request["result"] is None
    assert not request._event_done.is_set()


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="Interface close skips closing the server descriptor when shutdown raises",
)
def test_server_descriptor_closes_even_if_shutdown_fails(
    listening_interface, monkeypatch
):
    interface = listening_interface

    def fail_shutdown(*args):
        raise OSError("injected shutdown failure")

    monkeypatch.setattr(socket.socket, "shutdown", fail_shutdown)
    try:
        interface.close()
        assert interface.server.fileno() == -1
    finally:
        interface.server.close()


def test_poll_collects_only_readable_clients(connection, monkeypatch):
    first, first_peer = connection
    server, second_peer = socket.socketpair()
    second = Driver(server)
    server.close()
    requests = [request_for_active_atom(), request_for_active_atom()]
    interface = sockets.InterfaceSocket()
    interface.clients = [first, second]
    interface.requests = requests
    interface.jobs = [[r, c, None] for r, c in zip(requests, interface.clients)]
    for request in requests:
        request["start"] = -1
    first_peer.sendall(payload() + MESSAGE["ready"])
    second_peer.sendall(payload() + MESSAGE["ready"])
    real_select = sockets.select.select
    rounds = []

    def stagger_results(readable, writable, exceptional, timeout):
        # Control readiness order, but let Driver's status exchanges use real select.
        if timeout == 0.01 and readable:
            rounds.append(1)
            if len(rounds) == 1:
                return [first], [], []
            assert requests[0]._event_done.is_set()
            assert not requests[1]._event_done.is_set()
        return real_select(readable, writable, exceptional, timeout)

    monkeypatch.setattr(sockets.select, "select", stagger_results)
    try:
        interface._pool_distribute_select()
        assert len(rounds) == 2
        assert not interface.jobs
        for request in requests:
            assert request._event_done.is_set()
            assert request["result"][0] == 1.0
            np.testing.assert_array_equal(request["result"][1], [0, 0, 0, 0, 1, 2])
    finally:
        second.close()
        second_peer.close()


def test_collection_mode_stops_dispatch_when_no_client_can_accept_work(connection):
    driver, _ = connection
    driver.status = Status.Up | Status.HasData
    request = request_for_active_atom()
    request["status"] = "Queued"
    interface = sockets.InterfaceSocket()
    interface.clients = [driver]
    interface.requests = [request]
    interface.jobs = []
    interface._pool_distribute_select()
    assert request["status"] == "Queued"
    assert interface.prlist == [request]
    assert not interface.jobs
    interface.prlist = []
    assert not interface._dispatch_free_client_batch(driver)


def test_interface_close_releases_server_after_successful_shutdown(
    listening_interface, monkeypatch
):
    # Some platforms allow shutdown on a listening socket; others raise ENOTCONN.
    monkeypatch.setattr(socket.socket, "shutdown", lambda *args: None)
    listening_interface.close()
    assert listening_interface.server.fileno() == -1


def test_worker_pool_does_not_send_one_request_to_two_clients(
    listening_interface, connection
):
    interface = listening_interface
    interface.consolidate_messages = False
    driver, client = connection
    server, spare_peer = socket.socketpair()
    spare = Driver(server)
    server.close()
    request = request_for_active_atom()
    request["status"] = "Queued"
    request["start"] = -1  # ForceField.queue supplies this before worker submission
    interface.offset = 0.25
    driver.status = spare.status = Status.Up | Status.Ready
    interface.clients = [driver, spare]
    interface.requests = [request]
    client.sendall(
        MESSAGE["ready"] + MESSAGE["havedata"] + payload() + MESSAGE["ready"]
    )
    try:
        interface.pool_distribute()
        # The worker may have completed before the first pool check; either way wait with a bound.
        assert request._event_done.wait(2)
        for _, _, future in interface.jobs:
            future.result(timeout=2)
        interface.pool_distribute()
        assert not interface.jobs
        assert request["result"][0] == 1.0
        np.testing.assert_array_equal(request["result"][1], [0, 0, 0, 0, 1, 2])
        spare_peer.setblocking(False)
        with pytest.raises(BlockingIOError):
            spare_peer.recv(12)
    finally:
        interface.clients = []  # connection fixture owns driver; close spare here
        spare.close()
        spare_peer.close()

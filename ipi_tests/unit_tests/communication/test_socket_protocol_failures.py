"""Real socket pairs exercise malformed replies and interrupted force transfers."""

import socket

import numpy as np
import pytest
from ipi.interfaces import sockets
from ipi.interfaces.sockets import (
    Driver,
    Status,
    MESSAGE,
    Disconnected,
    InvalidSize,
    InvalidStatus,
)


@pytest.fixture
def connection():
    server, client = socket.socketpair()
    server.settimeout(0.05)
    driver = Driver(server)
    server.close()
    try:
        yield driver, client
    finally:
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
def test_malformed_force_payload_is_rejected(connection, data, error):
    driver, client = connection
    client.sendall(data)
    with pytest.raises(error):
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
def test_partial_force_reply_cannot_complete_request(connection, batched):
    driver, client = connection
    request = request_for_active_atom()
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

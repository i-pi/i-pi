"""Force requests preserve input snapshots and complete before being released."""

import numpy as np
import pytest

from ipi.engine.atoms import Atoms
from ipi.engine.cell import Cell
from ipi.engine.forcefields import ForceField, FFEval
from ipi.utils.softexit import Softexit


@pytest.mark.parametrize("provider", [ForceField, FFEval])
@pytest.mark.parametrize("threaded", [False, True])
def test_request_lifecycle(provider, threaded, monkeypatch):
    exit_state = Softexit()
    monkeypatch.setattr("ipi.engine.forcefields.softexit", exit_state)
    ff = provider(
        threaded=threaded,
        offset=0.25,
        dopbc=True,
        active=np.array([1]),
        pars={"label": "reference"},
    )
    atoms = Atoms(2)
    atoms.q[:] = [11, 0, 0, 0, -12, 0]
    cell = Cell(np.eye(3) * 10)
    ff.start()
    try:
        if threaded:
            with pytest.raises(NameError, match="already started"):
                ff.start()
        request = ff.queue(atoms, cell, reqid=7, template={"extra": "preserved"})
        duplicate = ff.queue(atoms, cell, reqid=7)
        assert not (request == duplicate)
        atoms.q[:] = 0
        cell.h[:] *= 2
        assert request._event_done.wait(timeout=5)
        # release must join an evaluation worker before its result is consumed.
        ff.release(request)
        assert ff.requests == [duplicate]
        assert duplicate._event_done.wait(timeout=5)
        ff.release(duplicate)
        assert not ff.requests
        ff.release(request, lock=False)  # releasing twice is harmless
        np.testing.assert_allclose(
            request["pos"], [1, 0, 0, 0, -2, 0], rtol=0, atol=3e-15
        )
        np.testing.assert_array_equal(request["active"], [3, 4, 5])
        np.testing.assert_array_equal(request["cell"][0], np.eye(3) * 10)
        assert request["id"] == 7 and request["extra"] == "preserved"
        assert "label : reference" in request["pars"]
        energy, force, virial, extras = request["result"]
        assert energy == pytest.approx(-0.25)
        np.testing.assert_array_equal(force, np.zeros(6))
        np.testing.assert_array_equal(virial, np.zeros((3, 3)))
        assert extras == {"raw": ""}
    finally:
        ff.stop()
        if ff._thread is not None:
            ff._thread.join(timeout=5)
            assert not ff._thread.is_alive()


@pytest.mark.parametrize("active", [[2], [0, 1, 2]])
def test_request_rejects_active_atoms_outside_system(active):
    ff = ForceField(active=np.array(active))
    with pytest.raises(ValueError, match="more active atoms"):
        ff.queue(Atoms(2), Cell(np.eye(3)))
    assert not ff.requests


def test_stop_marks_pending_requests_for_exit():
    ff = ForceField(threaded=True)
    request = ff.queue(Atoms(1), Cell(np.eye(3)))
    ff.stop()
    ff.poll()
    assert request["status"] == "Exit"
    assert request["result"] is None
    assert not request._event_done.is_set()
    ff.release(request)
    assert not ff.requests


@pytest.mark.parametrize("kind", ["committee", "rotations"])
def test_composite_forcefields_complete_threaded_requests(kind, monkeypatch):
    from ipi.engine.forcefields import FFCommittee, FFDirect, FFRotations
    from ipi.utils.prng import Random

    exit_state = Softexit()
    monkeypatch.setattr("ipi.engine.forcefields.softexit", exit_state)
    if kind == "committee":
        members = [
            FFDirect(pes="harmonic", pars={"k1": k}, threaded=True) for k in [1.0, 3.0]
        ]
        ff = FFCommittee(fflist=members)
    else:
        member = FFDirect(
            pes="harmonic", pars={"k1": 1.0, "k2": 2.0, "k3": 3.0}, threaded=True
        )
        ff = FFRotations(ffdirect=member, grid_order=3, prng=Random(seed=31))
    ff.bind(None)
    ff.start()
    try:
        atoms = Atoms(1)
        atoms.q[:] = [0.2, -0.3, 0.4]
        request = ff.queue(atoms, Cell(np.eye(3) * 10))
        assert request._event_done.wait(timeout=5)
        assert request["result"][0] == pytest.approx(0.29, abs=1e-12)
        np.testing.assert_allclose(request["result"][1], [-0.4, 0.6, -0.8], atol=1e-12)
        ff.release(request)
    finally:
        for _, control in exit_state.tlist:
            control[0] = False
        for thread, _ in exit_state.tlist:
            thread.join(timeout=5)
            assert not thread.is_alive()


def test_batched_socket_rejects_unconsolidated_protocol():
    from ipi.engine.forcefields import FFSocket
    from ipi.interfaces.sockets import InterfaceSocket

    ff = FFSocket(interface=InterfaceSocket(batch_size=2, consolidate_messages=False))
    with pytest.raises(ValueError, match="consolidate_messages"):
        ff.start()
    assert ff._thread is None


def test_socket_forcefield_stops_polling_and_removes_endpoint(monkeypatch):
    import socket
    import tempfile
    from pathlib import Path
    from ipi.engine.forcefields import FFSocket
    from ipi.interfaces.sockets import InterfaceSocket

    monkeypatch.setattr("ipi.engine.forcefields.softexit", Softexit())
    # A short path also respects the Unix socket path limit on macOS.
    with tempfile.TemporaryDirectory(prefix="ipi-ff-", dir="/tmp") as directory:
        endpoint = Path(directory) / "socket"
        ff = FFSocket(
            interface=InterfaceSocket(address="socket", sockets_prefix=directory + "/")
        )
        ff.start()
        try:
            assert endpoint.exists()
            with socket.socket(socket.AF_UNIX) as client:
                client.settimeout(5)
                client.connect(str(endpoint))
        finally:
            ff.stop()
        assert not ff._thread.is_alive()
        assert not endpoint.exists()

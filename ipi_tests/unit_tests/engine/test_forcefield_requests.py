"""Force requests preserve input snapshots, and are completed by the polls of
the scheduler."""

import socket
import tempfile
from pathlib import Path

import numpy as np
import pytest

from ipi.engine.atoms import Atoms
from ipi.engine.cell import Cell
from ipi.engine.forcefields import (
    FFCommittee,
    FFDirect,
    FFEval,
    FFRotations,
    FFSocket,
    ForceField,
)
from ipi.interfaces.sockets import InterfaceSocket
from ipi.utils import scheduler
from ipi.utils.prng import Random
from ipi.utils.softexit import Softexit


def done(request):
    return lambda: request["status"] == "Done"


@pytest.mark.parametrize("provider", [ForceField, FFEval])
def test_request_lifecycle(provider, monkeypatch):
    exit_state = Softexit()
    monkeypatch.setattr("ipi.engine.forcefields.softexit", exit_state)
    ff = provider(
        offset=0.25,
        dopbc=True,
        active=np.array([1]),
        pars={"label": "reference"},
    )
    atoms = Atoms(2)
    atoms.q_[:] = [11, 0, 0, 0, -12, 0]
    cell = Cell(np.eye(3) * 10)
    ff.start()
    try:
        request = ff.queue(atoms, cell, reqid=7, template={"extra": "preserved"})
        duplicate = ff.queue(atoms, cell, reqid=7)
        assert not (request == duplicate)
        # the requests hold copies of the positions and of the cell
        atoms.q_[:] = 0
        cell.h_[:] *= 2
        assert request["status"] == "Queued"
        scheduler.wait(done(request))
        ff.release(request)
        assert ff.requests == [duplicate]
        scheduler.wait(done(duplicate))
        ff.release(duplicate)
        assert not ff.requests
        ff.release(request)  # releasing twice is harmless
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


@pytest.mark.parametrize("active", [[2], [0, 1, 2]])
def test_request_rejects_active_atoms_outside_system(active):
    ff = ForceField(active=np.array(active))
    with pytest.raises(ValueError, match="more active atoms"):
        ff.queue(Atoms(2), Cell(np.eye(3)))
    assert not ff.requests


def test_socket_rejects_systems_of_different_sizes():
    """The active atoms are set at the first request, for one system size."""
    ff = FFSocket(interface=InterfaceSocket())
    ff.queue(Atoms(2), Cell(np.eye(3)))
    with pytest.raises(ValueError, match="systems with 2 and 3 atoms"):
        ff.queue(Atoms(3), Cell(np.eye(3)))


def test_stop_marks_pending_requests_for_exit():
    ff = ForceField()
    request = ff.queue(Atoms(1), Cell(np.eye(3)))
    ff.stop()
    assert not ff.poll(idle=True)
    assert request["status"] == "Exit"
    assert request["result"] is None
    ff.release(request)
    assert not ff.requests


@pytest.mark.parametrize("kind", ["committee", "rotations"])
def test_composite_forcefields_complete_requests(kind, monkeypatch):

    monkeypatch.setattr("ipi.engine.forcefields.softexit", Softexit())
    if kind == "committee":
        members = [FFDirect(pes="harmonic", pars={"k1": k}) for k in [1.0, 3.0]]
        ff = FFCommittee(fflist=members)
    else:
        member = FFDirect(pes="harmonic", pars={"k1": 1.0, "k2": 2.0, "k3": 3.0})
        ff = FFRotations(ffdirect=member, grid_order=3, prng=Random(seed=31))
    ff.bind(None)
    ff.start()
    atoms = Atoms(1)
    atoms.q_[:] = [0.2, -0.3, 0.4]
    request = ff.queue(atoms, Cell(np.eye(3) * 10))
    scheduler.wait(done(request))
    assert request["result"][0] == pytest.approx(0.29, abs=1e-12)
    np.testing.assert_allclose(request["result"][1], [-0.4, 0.6, -0.8], atol=1e-12)
    ff.release(request)


def test_socket_forcefield_stops_and_removes_endpoint(monkeypatch):

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
        assert ff.socket.server is None
        assert not endpoint.exists()

"""Unit tests for the FFMPI forcefield.

The MPI transport is exercised against a fake communicator that emulates a
harmonic driver, so these run without mpi4py or a launched MPI runtime. An
end-to-end test against a real `mpirun` lives in test_ffmpi_mpirun.py.
"""

import sys
import types

import numpy as np

from ipi.engine.forcefields import FFMPI, ForceRequest
from ipi.interfaces.utils import parse_extra
from ipi.interfaces.mpi import (
    MPI_TAG_INIT,
    MPI_TAG_WORK,
    MPI_TAG_RESULT,
    MPI_TAG_EXTRA,
    MPI_TAG_EXIT,
)


class _FakeStatus:
    def __init__(self):
        self._count = 0

    def Get_count(self, dtype):
        return self._count

    def Get_source(self):
        return 0


class _FakeMPI:
    """Minimal stand-in for the mpi4py.MPI namespace used by FFMPI.poll()."""

    INT = "int"
    DOUBLE = "double"
    BYTE = "byte"
    ANY_SOURCE = -1
    Status = _FakeStatus


class _FakeComm:
    """Emulates a pool of harmonic driver roots over the FFMPI protocol."""

    def __init__(self, k):
        self.k = k
        self._work = {}  # dest -> partial work (header then body)
        self._ready = {}  # source -> staged result buffers

    @staticmethod
    def _arr(buf):
        return buf[0] if isinstance(buf, list) else buf

    def Send(self, buf, dest=None, tag=None):
        if tag in (MPI_TAG_INIT, MPI_TAG_EXIT):
            return
        arr = self._arr(buf)
        if tag == MPI_TAG_WORK:
            st = self._work.get(dest)
            if st is None:
                self._work[dest] = {"nat": int(arr[0]), "batch": int(arr[1])}
            else:
                self._compute(dest, st["nat"], st["batch"], np.array(arr, float))
                del self._work[dest]

    def _compute(self, dest, nat, batch, body):
        si, so = 18 + 3 * nat, 1 + 3 * nat + 9
        numeric = np.zeros(batch * so)
        for b in range(batch):
            pos = body[b * si + 18 : b * si + 18 + 3 * nat]
            numeric[b * so] = 0.5 * self.k * (pos**2).sum()
            numeric[b * so + 1 : b * so + 1 + 3 * nat] = -self.k * pos
        # one zero-length extra string per structure
        ebuf = np.zeros(batch, np.int32).view(np.uint8).copy()
        self._ready[dest] = {"numeric": numeric, "ebuf": ebuf}

    def Iprobe(self, source=None, tag=None, status=None):
        return source in self._ready

    def Probe(self, source=None, tag=None, status=None):
        status._count = len(self._ready[source]["ebuf"])

    def Recv(self, buf, source=None, tag=None):
        arr = self._arr(buf)
        r = self._ready[source]
        if tag == MPI_TAG_RESULT:
            arr[:] = r["numeric"]
        elif tag == MPI_TAG_EXTRA:
            arr[:] = r["ebuf"]
            del self._ready[source]


def _make_request(pos):
    return ForceRequest(
        {
            "id": 0,
            "pos": np.asarray(pos, float),
            "cell": (np.eye(3), np.eye(3)),
            "pars": " ",
            "result": None,
            "status": "Queued",
            "t_dispatched": 0,
            "t_finished": 0,
        }
    )


def _make_ffmpi(monkeypatch, k=2.0, ranks=(1, 2), batch_size=1):
    fake = types.ModuleType("mpi4py")
    fake.MPI = _FakeMPI
    fake.rc = types.SimpleNamespace(thread_level="multiple")
    monkeypatch.setitem(sys.modules, "mpi4py", fake)

    obj = FFMPI(name="mpi", batch_size=batch_size)
    # wire the fake communicator and pre-built rank pool onto the interface,
    # bypassing the MPI handshake that interface.open() would otherwise run
    obj.interface.comm = _FakeComm(k)
    obj.interface.driver_ranks = list(ranks)
    obj.interface._ranks = {
        r: {"rank": r, "state": "IDLE", "requests": [], "nat": None} for r in ranks
    }
    return obj


def test_ffmpi_single_requests(monkeypatch):
    k = 2.0
    obj = _make_ffmpi(monkeypatch, k=k, ranks=(1, 2))
    positions = [np.array([0.1, -0.2, 0.3]), np.array([0.4, 0.5, -0.6])]
    reqs = [_make_request(p) for p in positions]
    obj.requests.extend(reqs)

    obj.poll()  # dispatch to the two idle roots
    obj.poll()  # collect the staged results

    for req, pos in zip(reqs, positions):
        assert req["status"] == "Done"
        energy, force, virial, extra = req["result"]
        assert np.isclose(energy, 0.5 * k * (pos**2).sum())
        assert np.allclose(force, -k * pos)
        assert np.allclose(virial, 0.0)


def test_ffmpi_batched_request(monkeypatch):
    k = 1.5
    obj = _make_ffmpi(monkeypatch, k=k, ranks=(1,), batch_size=3)
    positions = [
        np.array([0.1, 0.0, 0.0]),
        np.array([0.0, 0.2, 0.0]),
        np.array([0.0, 0.0, 0.3]),
    ]
    reqs = [_make_request(p) for p in positions]
    obj.requests.extend(reqs)

    obj.poll()  # all three bundled into one request to the single root
    obj.poll()  # collect

    assert obj.interface._ranks[1]["state"] == "IDLE"
    for req, pos in zip(reqs, positions):
        assert req["status"] == "Done"
        energy, force, _, _ = req["result"]
        assert np.isclose(energy, 0.5 * k * (pos**2).sum())
        assert np.allclose(force, -k * pos)


def test_ffmpi_offset_applied(monkeypatch):
    obj = _make_ffmpi(monkeypatch, k=1.0, ranks=(1,))
    obj.interface.offset = 5.0
    req = _make_request(np.array([1.0, 0.0, 0.0]))
    obj.requests.append(req)

    obj.poll()
    obj.poll()

    energy = req["result"][0]
    assert np.isclose(energy, 0.5 * 1.0 * 1.0 - 5.0)


def test_parse_extra():
    extra = parse_extra('{"a": 1}')
    assert extra["a"] == 1
    assert extra["raw"] == '{"a": 1}'

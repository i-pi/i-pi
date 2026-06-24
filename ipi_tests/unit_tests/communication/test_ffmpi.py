"""Unit tests for the FFMPI forcefield.

The MPI transport is exercised against a fake communicator that emulates a
harmonic driver, so these run without mpi4py or a launched MPI runtime. An
end-to-end test against a real `mpirun` lives in test_ffmpi_mpirun.py.
"""

import sys
import types

import numpy as np
import pytest

import ipi.interfaces.mpi as mpi_mod
from ipi.engine.forcefields import FFMPI, ForceRequest
from ipi.inputs.forcefields import InputFFMPI
from ipi.inputs.simulation import InputSimulation
from ipi.interfaces.utils import parse_extra
from ipi.utils.io.inputs.io_xml import xml_parse_string
from ipi.interfaces.mpi import (
    MPIWorldManager,
    InterfaceMPI,
    get_mpi_world_manager,
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
    # valid JSON: fields returned plus the literal under 'raw'
    extra = parse_extra('{"a": 1}')
    assert extra["a"] == 1
    assert extra["raw"] == '{"a": 1}'

    # empty / whitespace -> empty dict
    assert parse_extra("") == {}
    assert parse_extra("   ") == {}

    # non-JSON string -> empty dict that still exposes the literal under 'raw'
    extra = parse_extra("not json")
    assert extra == {"raw": "not json"}

    # a JSON object that itself carries a 'raw' key is rejected
    with pytest.raises(ValueError):
        parse_extra('{"raw": 1}')


# --- MPIWorldManager: rank partitioning logic (no MPI runtime needed) ---


def test_assign_lone_ff_claims_all_drivers():
    """A single ffmpi claims every driver root, tagged or not."""
    mgr = MPIWorldManager()
    mgr.register("solo")
    assigned = mgr._assign({"": [3, 1], "ignored-tag": [2]})
    assert assigned == {"solo": [1, 2, 3]}


def test_assign_multiple_ff_tagged():
    mgr = MPIWorldManager()
    mgr.register("ffA")
    mgr.register("ffB")
    assigned = mgr._assign({"ffA": [1], "ffB": [3, 2]})
    assert assigned == {"ffA": [1], "ffB": [2, 3]}


def test_assign_multiple_ff_untagged_raises():
    mgr = MPIWorldManager()
    mgr.register("ffA")
    mgr.register("ffB")
    with pytest.raises(ValueError, match="every driver must pass --mpi-name"):
        mgr._assign({"ffA": [1], "": [2]})


def test_assign_unknown_name_raises():
    mgr = MPIWorldManager()
    mgr.register("ffA")
    mgr.register("ffB")
    with pytest.raises(ValueError, match="not an"):
        mgr._assign({"ffA": [1], "ffC": [2]})


def test_assign_missing_driver_raises():
    mgr = MPIWorldManager()
    mgr.register("ffA")
    mgr.register("ffB")
    with pytest.raises(ValueError, match="no driver ranks"):
        mgr._assign({"ffA": [1]})


def test_register_idempotent_and_roots_for():
    mgr = MPIWorldManager()
    mgr.register("ffA")
    mgr.register("ffA")
    assert mgr.registered_names == {"ffA"}
    mgr.roots_by_name = {"ffA": [1, 2]}
    assert mgr.roots_for("ffA") == [1, 2]
    # an unknown name yields an empty (fresh) list
    assert mgr.roots_for("ffZ") == []
    mgr.roots_for("ffA").append(99)
    assert mgr.roots_for("ffA") == [1, 2]  # returns a copy, not the internal list


def test_get_mpi_world_manager_singleton(monkeypatch):
    monkeypatch.setattr(mpi_mod, "_mpi_world_manager", None)
    mgr1 = get_mpi_world_manager()
    mgr2 = get_mpi_world_manager()
    assert mgr1 is mgr2


# --- MPIWorldManager.setup / _handshake driven by a fake communicator ---


class _HandshakeStatus:
    def __init__(self):
        self.src = 0
        self.n_int = 0
        self.n_byte = 0

    def Get_source(self):
        return self.src

    def Get_count(self, dtype):
        return self.n_int if dtype == "int" else self.n_byte


class _HandshakeComm:
    """Replays one HELLO (int payload + name bytes) per driver root."""

    def __init__(self, size, drivers):
        self._size = size
        self._drivers = drivers  # list of (src, [root, gsize, *members], name)
        self._i = 0
        self._stage = "int"
        self.split_called = False

    def Get_rank(self):
        return 0

    def Get_size(self):
        return self._size

    def Split(self, color=0, key=0):
        self.split_called = True

    def Probe(self, source=None, tag=None, status=None):
        src, ints, name = self._drivers[self._i]
        status.src = src
        status.n_int = len(ints)
        status.n_byte = len(name.encode("utf-8"))

    def Recv(self, buf, source=None, tag=None):
        arr = buf[0] if isinstance(buf, list) else buf
        _, ints, name = self._drivers[self._i]
        if self._stage == "int":
            arr[:] = np.array(ints, dtype=np.int32)
            self._stage = "byte"
        else:
            arr[:] = np.frombuffer(name.encode("utf-8"), dtype=np.uint8)
            self._stage = "int"
            self._i += 1


def _fake_mpi(comm):
    return types.SimpleNamespace(
        INT="int",
        DOUBLE="double",
        BYTE="byte",
        ANY_SOURCE=-1,
        COMM_WORLD=comm,
        Status=_HandshakeStatus,
    )


def test_setup_handshake_partitions_drivers():
    comm = _HandshakeComm(
        size=3, drivers=[(1, [1, 1, 1], "ffA"), (2, [2, 1, 2], "ffB")]
    )
    mgr = MPIWorldManager()
    mgr.register("ffA")
    mgr.register("ffB")
    mgr.setup(_fake_mpi(comm))
    assert comm.split_called
    assert mgr.roots_for("ffA") == [1]
    assert mgr.roots_for("ffB") == [2]
    # a second setup() is a no-op (the collective split must happen only once)
    comm.split_called = False
    mgr.setup(_fake_mpi(comm))
    assert not comm.split_called


def test_setup_requires_rank_zero():
    comm = _HandshakeComm(size=2, drivers=[])
    comm.Get_rank = lambda: 1
    mgr = MPIWorldManager()
    mgr.register("ffA")
    with pytest.raises(ValueError, match="rank 0"):
        mgr.setup(_fake_mpi(comm))


def test_setup_requires_driver_ranks():
    comm = _HandshakeComm(size=1, drivers=[])
    mgr = MPIWorldManager()
    mgr.register("ffA")
    with pytest.raises(ValueError, match="no driver ranks"):
        mgr.setup(_fake_mpi(comm))


# --- InterfaceMPI / FFMPI construction and shutdown ---


def test_interface_registers_name(monkeypatch):
    monkeypatch.setattr(mpi_mod, "_mpi_world_manager", None)
    iface = InterfaceMPI(name="ffX", batch_size=2)
    assert "ffX" in get_mpi_world_manager().registered_names
    assert iface._mpi_lock is get_mpi_world_manager().mpi_lock


class _ExitComm:
    def __init__(self):
        self.exits = []

    def Send(self, buf, dest=None, tag=None):
        if tag == MPI_TAG_EXIT:
            self.exits.append(dest)


def test_interface_close_sends_exit(monkeypatch):
    monkeypatch.setattr(mpi_mod, "_mpi_world_manager", None)
    monkeypatch.setitem(sys.modules, "mpi4py", types.SimpleNamespace(MPI=_FakeMPI))
    iface = InterfaceMPI(name="ffX")
    comm = _ExitComm()
    iface.comm = comm
    iface.driver_ranks = [1, 2]
    iface.close()
    assert comm.exits == [1, 2]


def test_interface_close_without_comm_is_noop(monkeypatch):
    monkeypatch.setattr(mpi_mod, "_mpi_world_manager", None)
    iface = InterfaceMPI(name="ffX")
    iface.comm = None
    iface.close()  # must not raise


def test_ffmpi_requires_threaded():
    with pytest.raises(ValueError, match="threaded=True"):
        FFMPI(name="mpi", threaded=False)


def test_ffmpi_rejects_unknown_mode():
    with pytest.raises(ValueError, match="Unknown ffmpi mode"):
        FFMPI(name="mpi", mode="bogus")


def test_ffmpi_uses_injected_interface():
    iface = types.SimpleNamespace(requests=None, offset=None)
    obj = FFMPI(name="mpi", offset=3.0, interface=iface)
    assert obj.interface is iface
    assert iface.requests is obj.requests
    assert iface.offset == obj.offset


def test_input_ffmpi_parses_xml_parameters():
    root = xml_parse_string("""
        <ffmpi name='mpi_xml' mode='mpi' threaded='true'>
          <latency>0.2</latency>
          <offset>1.5</offset>
          <batch_size>3</batch_size>
        </ffmpi>
        """)
    node = root.fields[0][1]
    parsed = InputFFMPI()

    parsed.parse(node)
    fetched = parsed.fetch()

    assert fetched.name == "mpi_xml"
    assert fetched.mode == "mpi"
    assert fetched.interface.batch_size == 3
    assert fetched.latency == 0.2
    assert fetched.offset == 1.5
    assert fetched.threaded


def test_simulation_input_routes_ffmpi_xml():
    root = xml_parse_string("""
        <simulation>
          <ffmpi name='mpi_xml' mode='mpi'>
            <batch_size>3</batch_size>
          </ffmpi>
        </simulation>
        """)
    node = root.fields[0][1]
    parsed = InputSimulation()

    parsed.parse(node)

    tag, forcefield = parsed.extra[0]
    assert tag == "ffmpi"
    assert isinstance(forcefield, InputFFMPI)
    assert forcefield.name.fetch() == "mpi_xml"
    assert forcefield.batch_size.fetch() == 3

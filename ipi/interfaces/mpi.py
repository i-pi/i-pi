"""Interface between i-PI and force drivers that communicate over MPI.

This is the MPI analogue of :mod:`ipi.interfaces.sockets`: the :class:`FFMPI`
forcefield owns an :class:`InterfaceMPI` that exchanges positions and forces
with a pool of co-launched driver group roots. A process-wide
:class:`MPIWorldManager` performs the one-time MPI_COMM_WORLD setup shared by all
the ``<ffmpi>`` forcefields and partitions the driver ranks among them by name.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

import sys
import time
import threading

import numpy as np

from ipi.utils.messages import info, verbosity
from ipi.utils.depend import dstrip
from ipi.interfaces.utils import parse_extra

# MPI message tags for the ffmpi protocol (must match the driver side, see
# drivers/py/driver.py and drivers/f90/fmpi.f90)
MPI_TAG_HELLO = 1
MPI_TAG_INIT = 2
MPI_TAG_WORK = 3
MPI_TAG_RESULT = 4
MPI_TAG_EXTRA = 5
MPI_TAG_EXIT = 6


class MPIWorldManager:
    """Process-wide owner of the one-time MPI_COMM_WORLD setup shared by every
    FFMPI forcefield.

    Several ``<ffmpi>`` forcefields live in the same MPMD job (i-PI is rank 0,
    the driver ranks are co-launched). This singleton performs, exactly once:
    the single collective ``COMM_WORLD.Split(color=0)`` (so the drivers can form
    their own sub-communicators), and the handshake that discovers every driver
    group root together with the forcefield name it was launched for. Each
    interface then claims the roots tagged with its own name. A shared lock
    serialises all MPI calls across the per-forcefield poll threads, so
    MPI_THREAD_SERIALIZED is enough and the (disjoint) point-to-point traffic
    never collides.
    """

    def __init__(self):
        self.comm = None
        self._setup_done = False
        self.roots_by_name = {}  # forcefield name -> sorted list of root ranks
        self.registered_names = set()
        self.mpi_lock = threading.Lock()
        self._setup_lock = threading.Lock()

    def register(self, name):
        """Records a forcefield name before the handshake (idempotent)."""
        self.registered_names.add(name)

    def setup(self, MPI):
        """One-time collective split + handshake; later calls are no-ops."""
        with self._setup_lock:
            if self._setup_done:
                return
            self.comm = MPI.COMM_WORLD
            if self.comm.Get_rank() != 0:
                raise ValueError("ffmpi expects i-PI to be rank 0 of MPI_COMM_WORLD.")
            if self.comm.Get_size() < 2:
                raise ValueError(
                    "ffmpi found no driver ranks: launch with "
                    "`mpirun -n 1 i-pi input.xml : -n N <driver> --mpi ...`."
                )
            # i-PI sits alone in color 0; the driver ranks build their own pools
            MPI.COMM_WORLD.Split(color=0, key=0)
            self._handshake(MPI)
            self._setup_done = True

    def _handshake(self, MPI):
        """Collects one HELLO per driver group root: an int payload
        [root_rank, group_size, members...] followed by a byte message with the
        forcefield name (empty = untagged). Buckets the roots by name."""

        expected = set(range(1, self.comm.Get_size()))
        covered = set()
        roots_by_name = {}
        status = MPI.Status()
        while covered != expected:
            self.comm.Probe(source=MPI.ANY_SOURCE, tag=MPI_TAG_HELLO, status=status)
            src = status.Get_source()
            buf = np.empty(status.Get_count(MPI.INT), dtype=np.int32)
            self.comm.Recv([buf, MPI.INT], source=src, tag=MPI_TAG_HELLO)
            # the forcefield name follows from the same source, same tag
            self.comm.Probe(source=src, tag=MPI_TAG_HELLO, status=status)
            nbuf = np.empty(status.Get_count(MPI.BYTE), dtype=np.uint8)
            self.comm.Recv([nbuf, MPI.BYTE], source=src, tag=MPI_TAG_HELLO)
            name = nbuf.tobytes().decode("utf-8") if nbuf.size else ""
            roots_by_name.setdefault(name, []).append(int(buf[0]))
            covered |= set(buf[2:].tolist())

        self.roots_by_name = self._assign(roots_by_name)

    def _assign(self, roots_by_name):
        """Resolves the driver-name buckets against the registered forcefields."""
        if len(self.registered_names) == 1:
            # a lone ffmpi claims every driver, tagged or not (back-compat)
            only = next(iter(self.registered_names))
            allroots = sorted(r for rs in roots_by_name.values() for r in rs)
            return {only: allroots}
        if "" in roots_by_name:
            raise ValueError(
                "with multiple ffmpi forcefields every driver must pass --mpi-name; "
                "found untagged driver ranks."
            )
        for name in roots_by_name:
            if name not in self.registered_names:
                raise ValueError(
                    "drivers were launched with --mpi-name '%s', which is not an "
                    "ffmpi forcefield in the input." % name
                )
        for name in self.registered_names:
            if not roots_by_name.get(name):
                raise ValueError(
                    "ffmpi '%s' was given no driver ranks; check --mpi-name and the "
                    "mpirun launch line." % name
                )
        return {k: sorted(v) for k, v in roots_by_name.items()}

    def roots_for(self, name):
        return list(self.roots_by_name.get(name, []))


_mpi_world_manager = None


def get_mpi_world_manager():
    """Returns the process-wide MPIWorldManager, creating it on first use."""
    global _mpi_world_manager
    if _mpi_world_manager is None:
        _mpi_world_manager = MPIWorldManager()
    return _mpi_world_manager


class InterfaceMPI:
    """Pool of MPI driver group roots serving one FFMPI forcefield.

    Mirrors :class:`ipi.interfaces.sockets.InterfaceSocket`: it owns the
    communicator and the per-root slot pool, dispatches queued requests and
    collects the results. The ``requests`` list and ``offset`` are linked from
    the owning FFMPI. Each root is a slot that handles up to ``batch_size``
    structures per request.
    """

    def __init__(self, name="", batch_size=1):
        self.name = name
        self.batch_size = batch_size
        self.requests = None  # linked to the FFMPI request list
        self.offset = 0.0  # linked from the FFMPI
        self.comm = None
        self.driver_ranks = []
        self._ranks = {}
        self._initialized = False
        # register the name now, before any open()/handshake, so the manager
        # knows the full set of forcefields when it partitions the drivers; the
        # shared MPI lock is also taken from the manager (usable before open())
        mgr = get_mpi_world_manager()
        mgr.register(self.name)
        self._mpi_lock = mgr.mpi_lock

    def open(self):
        """Joins the shared MPI world and claims this forcefield's driver ranks."""

        try:
            import mpi4py
        except ImportError:
            info("ffmpi requires a functioning mpi4py installation")
            raise

        # MPI calls happen on the main thread (setup) and then on the poll
        # thread, one at a time; request a thread level that allows that.
        mpi4py.rc.thread_level = "multiple"
        from mpi4py import MPI

        mgr = get_mpi_world_manager()
        try:
            mgr.setup(MPI)  # one collective split + handshake, shared by all FFMPI
            self.comm = mgr.comm
            self._mpi_lock = mgr.mpi_lock
            self.driver_ranks = mgr.roots_for(self.name)
            if not self.driver_ranks:
                raise ValueError(
                    "ffmpi '%s' was given no driver ranks; launch its driver(s) "
                    "with `--mpi --mpi-name %s`." % (self.name, self.name)
                )
        except Exception as err:
            # a setup/launch mismatch is fatal and leaves the driver ranks
            # blocked; abort the whole MPI job so it exits (with the message)
            # instead of hanging. Print to stderr so it shows at any verbosity.
            sys.stderr.write(
                " @ffmpi (%s): fatal MPI setup error: %s\n" % (self.name, err)
            )
            sys.stderr.flush()
            MPI.COMM_WORLD.Abort(1)
            raise
        self._ranks = {
            r: {"rank": r, "state": "IDLE", "requests": [], "nat": None}
            for r in self.driver_ranks
        }
        info(
            " @ForceField (%s): connected to MPI driver roots %s"
            % (self.name, self.driver_ranks),
            verbosity.low,
        )

    def poll(self):
        """Dispatches queued requests to idle driver roots and collects the
        results of the ones that have finished."""

        from mpi4py import MPI

        # the shared lock serialises MPI across the poll threads of every FFMPI
        with self._mpi_lock:
            if not self._initialized:
                self._send_init(MPI)
                self._initialized = True
            self._collect(MPI)
            self._dispatch(MPI)

    def _send_init(self, MPI):
        """Sends the INIT to each driver root, mirroring the socket protocol:
        the request id and the parameter string (the same one ForceField.queue
        built, with the batch_size token appended). The bundled drivers parse the
        batch size; a custom client can use the rest as initialization data."""

        pars = self.requests[0]["pars"] if self.requests else " "
        if self.batch_size > 1:
            sep = "," if pars.strip() else ""
            pars = pars + "%s batch_size:%d" % (sep, self.batch_size)
        rid = int(self.requests[0]["id"]) if self.requests else 0
        pbytes = np.frombuffer(pars.encode("utf-8"), dtype=np.uint8)
        header = np.array([rid, pbytes.size], dtype=np.int32)
        for r in self.driver_ranks:
            self.comm.Send([header, MPI.INT], dest=r, tag=MPI_TAG_INIT)
            self.comm.Send([pbytes, MPI.BYTE], dest=r, tag=MPI_TAG_INIT)

    def _collect(self, MPI):
        """Receives results from any busy slot whose driver has replied."""

        status = MPI.Status()
        for slot in self._ranks.values():
            if slot["state"] != "BUSY":
                continue
            if not self.comm.Iprobe(source=slot["rank"], tag=MPI_TAG_RESULT):
                continue

            reqs = slot["requests"]
            nat = slot["nat"]
            stride = 1 + 3 * nat + 9
            numeric = np.empty(len(reqs) * stride, dtype=np.float64)
            self.comm.Recv(
                [numeric, MPI.DOUBLE], source=slot["rank"], tag=MPI_TAG_RESULT
            )

            # the per-structure 'extra' strings follow in a single byte message
            self.comm.Probe(source=slot["rank"], tag=MPI_TAG_EXTRA, status=status)
            ebuf = np.empty(status.Get_count(MPI.BYTE), dtype=np.uint8)
            self.comm.Recv([ebuf, MPI.BYTE], source=slot["rank"], tag=MPI_TAG_EXTRA)

            self._unpack_results(reqs, nat, numeric, ebuf)
            slot["state"] = "IDLE"
            slot["requests"] = []

    def _unpack_results(self, reqs, nat, numeric, ebuf):
        """Splits the packed numeric buffer and extra strings back into the
        individual force requests."""

        stride = 1 + 3 * nat + 9
        epos = 0
        for i, req in enumerate(reqs):
            base = i * stride
            pot = float(numeric[base])
            force = numeric[base + 1 : base + 1 + 3 * nat].copy()
            vir = numeric[base + 1 + 3 * nat : base + stride].reshape(3, 3).copy()
            elen = int(np.frombuffer(ebuf[epos : epos + 4], dtype=np.int32)[0])
            epos += 4
            extra = ebuf[epos : epos + elen].tobytes().decode("utf-8") if elen else ""
            epos += elen

            req["result"] = [pot - self.offset, force, vir, parse_extra(extra)]
            req["status"] = "Done"
            req["t_finished"] = time.time()
            req._event_done.set()

    def _dispatch(self, MPI):
        """Sends queued requests to idle slots, up to batch_size per slot."""

        queued = [r for r in self.requests if r["status"] == "Queued"]
        if not queued:
            return

        qi = 0
        for slot in self._ranks.values():
            if slot["state"] != "IDLE" or qi >= len(queued):
                continue
            batch = queued[qi : qi + self.batch_size]
            qi += len(batch)

            nat = len(batch[0]["pos"]) // 3
            slot["nat"] = nat
            header = np.array([nat, len(batch)], dtype=np.int32)
            body = np.empty(len(batch) * (18 + 3 * nat), dtype=np.float64)
            for i, req in enumerate(batch):
                h, ih = req["cell"]
                off = i * (18 + 3 * nat)
                body[off : off + 9] = np.asarray(h).reshape(-1)
                body[off + 9 : off + 18] = np.asarray(ih).reshape(-1)
                body[off + 18 : off + 18 + 3 * nat] = dstrip(req["pos"]).reshape(-1)
                req["status"] = "Running"
                req["t_dispatched"] = time.time()

            self.comm.Send([header, MPI.INT], dest=slot["rank"], tag=MPI_TAG_WORK)
            self.comm.Send([body, MPI.DOUBLE], dest=slot["rank"], tag=MPI_TAG_WORK)
            slot["requests"] = batch
            slot["state"] = "BUSY"

    def close(self):
        """Tells every driver root to exit.

        At normal completion all slots are idle (the run ends only once every
        force request is collected); a driver still computing on an abnormal
        exit may hang, in which case the job has to be aborted."""

        if self.comm is None:
            return
        from mpi4py import MPI

        buf = np.zeros(1, dtype=np.int32)
        # other forcefields' poll threads may still be live, so serialise the
        # EXIT sends through the shared MPI lock
        with self._mpi_lock:
            for r in self.driver_ranks:
                self.comm.Send([buf, MPI.INT], dest=r, tag=MPI_TAG_EXIT)

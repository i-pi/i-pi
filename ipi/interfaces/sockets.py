"""Deals with the socket communication between i-PI and the drivers.

i-PI listens on a socket, to which the drivers connect as clients. The
exchange with a client follows the i-PI protocol, made of 12-byte headers
followed by data:

    STATUS    -> READY, NEEDINIT or HAVEDATA
    INIT      + request id + length + initialization string
    POSDATA   + cell + inverse cell + number of atoms + positions
    GETFORCE  -> FORCEREADY + energy + number of atoms + forces + virial + extras

i-PI never waits for a client: it reads a reply only once the socket of the
client can be read, and `Driver.expect` records the reply that is owed. A
client goes through these states:

    event          i-PI sends                       expect
    accepted       STATUS                           "status"
    status read    nothing; free if READY/NEEDINIT  None
    job sent       [INIT +] POSDATA + GETFORCE      "forces"
    forces read    STATUS                           "status"

A client that closes the connection, breaks the protocol, returns the forces
of a different number of atoms or exceeds the time-out is dropped, and the
requests it was computing are queued again for the other clients.

`InterfaceSocket.poll` makes all the progress that is possible without
waiting. `InterfaceSocket.fds` returns the sockets on which a connection or a
reply is expected: the scheduler (`ipi.utils.scheduler`) waits on them when
nothing else can progress.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.


import os
import socket
import select
import time

import numpy as np

from multiprocessing import shared_memory

from ipi.utils.messages import verbosity, warning, info
from ipi.utils.softexit import softexit
from ipi.interfaces.utils import parse_extra

__all__ = ["InterfaceSocket"]


HDRLEN = 12
# time-out of the operations on the socket of a client: the sends, and the
# reads of a reply that is being received
TIMEOUT = 0.1
# number of time-outs after which a reply that is being received is abandoned
NTIMEOUT = 20
# length of the queue of the clients waiting to be accepted
BACKLOG = 64


def Message(mystr):
    """Returns a header of standard length HDRLEN."""

    # convert to bytestream since we'll be sending this over a socket
    return str.ljust(str.upper(mystr), HDRLEN).encode()


MESSAGE = {
    msg: Message(msg)
    for msg in [
        "exit",
        "status",
        "ready",
        "havedata",
        "init",
        "needinit",
        "posdata",
        "getforce",
        "forceready",
    ]
}


class Disconnected(Exception):
    """Raised if a client has disconnected, or has broken the protocol."""

    pass


class InvalidSize(Exception):
    """Raised if a client returns forces for a different number of atoms."""

    pass


class Status(object):
    """The status of a client, as bit flags.

    Attributes:
       Disconnected: The client is not connected.
       Up: The client is connected.
       Ready: The client can receive positions.
       NeedsInit: The client must be initialized before receiving positions.
       HasData: The client has results to send.
    """

    Disconnected = 0
    Up = 1
    Ready = 2
    NeedsInit = 4
    HasData = 8


class DriverSocket(socket.socket):
    """A socket with the helpers of the i-PI protocol.

    Attributes:
       _buf: A buffer for the data received from the other end.
    """

    def __init__(self, sock):
        """Takes over the descriptor of `sock`, which is left detached.

        Args:
           sock: A connected socket.
        """

        timeout = sock.gettimeout()
        super(DriverSocket, self).__init__(
            sock.family, sock.type, sock.proto, fileno=sock.detach()
        )
        self.settimeout(timeout)
        self._buf = np.zeros(0, np.byte)
        try:
            self.peername = self.getpeername()
        except OSError:
            # a client that closed the connection before it was accepted (macOS
            # raises EINVAL for unix sockets): its handshake fails, and drops it
            self.peername = "disconnected client"

    def sendall(self, data):
        """Sends all of `data`, giving up only if the other end reads nothing
        for NTIMEOUT consecutive time-outs of the socket.

        `socket.sendall` bounds the total duration of the send by the
        time-out: it would give up on a client that is still reading a message
        too large to be sent within it (large systems, slow networks), and
        leave part of the message on the socket.

        Raises:
           socket.timeout: If the other end reads nothing for NTIMEOUT
              time-outs.
        """

        view = memoryview(data).cast("B")
        ntimeout = 0
        while len(view) > 0:
            try:
                sent = self.send(view)
            except socket.timeout:
                ntimeout += 1
                if ntimeout > NTIMEOUT:
                    raise
                continue
            ntimeout = 0
            view = view[sent:]

    def send_msg(self, msg):
        """Sends the header `msg`."""

        return self.sendall(MESSAGE[msg])

    def recv_msg(self, length=HDRLEN):
        """Receives a header."""

        return bytes(self.recvall(np.zeros(length, np.dtype("S1"))))

    def recvall(self, dest):
        """Fills the array `dest` with data received from the socket, and
        returns it (or its value, for a 0-dimensional array).

        Raises:
           Disconnected: If the other end has closed the connection, or sends
              nothing for NTIMEOUT consecutive time-outs.
        """

        blen = dest.itemsize * dest.size
        if blen > len(self._buf):
            self._buf = np.zeros(blen, np.byte)
        view = memoryview(self._buf)
        bpos = 0
        ntimeout = 0
        while bpos < blen:
            try:
                bpart = self.recv_into(view[bpos:blen], blen - bpos)
            except socket.timeout:
                ntimeout += 1
                if ntimeout > NTIMEOUT:
                    warning(
                        " @SOCKET:  Couldn't receive within %5d attempts. Time to give up!"
                        % (NTIMEOUT),
                        verbosity.low,
                    )
                    raise Disconnected()
                continue
            except OSError:
                raise Disconnected()
            if bpart == 0:
                raise Disconnected()
            # a reply that keeps arriving is not cut by the time-outs
            ntimeout = 0
            bpos += bpart

        if dest.ndim > 0:
            dest[:] = np.frombuffer(self._buf[0:blen], dest.dtype).reshape(dest.shape)
            return dest
        return np.frombuffer(self._buf[0:blen], dest.dtype)[0]


class Driver(DriverSocket):
    """A client connected to i-PI. Messages are sent in one go, and
    `InterfaceSocket` reads a reply only once it starts to arrive, so that
    i-PI never waits for the client to compute.

    Attributes:
       status: The status of the client (`Status` flags), as last reported.
       expect: The reply that the client owes, "status" or "forces", or None
          if it owes nothing.
       job: The requests that the client is computing, or None.
       start: The time at which the job was sent.
       lastreq: The id of the last request computed by the client.
       locked: Whether the client has been computing the same request id.
       batch_size: The number of structures sent together to the client.
       exit_on_disconnect: Whether i-PI exits when the client disconnects.
    """

    def __init__(self, sock):
        super(Driver, self).__init__(sock)
        self.status = Status.Up
        self.expect = None
        self.job = None
        self.start = 0.0
        self.lastreq = None
        self.locked = False
        self.batch_size = 1
        self.exit_on_disconnect = False

    def free(self):
        """Whether the client can receive a job: it owes no reply, and its
        last status was READY or NEEDINIT."""

        return self.expect is None and bool(
            self.status & (Status.Ready | Status.NeedsInit)
        )

    def ask_status(self):
        """Sends STATUS. The reply is read by `read_status`, once it arrives."""

        self.send_msg("status")
        self.expect = "status"

    def read_status(self):
        """Reads the reply to STATUS, and records the status.

        Raises:
           Disconnected: If the reply is not a known status.
        """

        reply = self.recv_msg()
        self.expect = None
        if reply == MESSAGE["ready"]:
            self.status = Status.Up | Status.Ready
        elif reply == MESSAGE["needinit"]:
            self.status = Status.Up | Status.NeedsInit
        elif reply == MESSAGE["havedata"]:
            self.status = Status.Up | Status.HasData
        else:
            raise Disconnected("unrecognized status reply %s" % reply)

    def send_job(self, requests):
        """Sends the structures of `requests` to the client, preceded by the
        initialization if the client needs it, and asks for the forces, all
        in one message. The forces are read by `read_job`, once they arrive."""

        message = b""
        if self.status & Status.NeedsInit:
            message += self._init_message(requests[0]["id"], requests[0]["pars"])
        message += self._posdata_message(requests) + MESSAGE["getforce"]
        self.sendall(message)
        self.status = Status.Up | Status.HasData
        self.expect = "forces"
        self.job = requests
        self.start = time.time()

    def read_job(self):
        """Reads the reply to GETFORCE, and returns the results of the job as
        a list of [pot, f, vir, extras], one for each request. The client
        owes nothing afterwards, and is bound to the id of the last request."""

        requests = self.job
        natoms = len(requests[0]["active"]) // 3
        if self.recv_msg() != MESSAGE["forceready"]:
            raise Disconnected("unexpected reply to GETFORCE")
        if self.batch_size > 1:
            results = self._read_batch(natoms, len(requests))
        else:
            results = [self._read_single(natoms)]
        self.expect = None
        self.job = None
        self.lastreq = requests[-1]["id"]
        return results

    def _init_message(self, rid, pars):
        """Returns the INIT message, with the batch size in the parameters of
        a batched client."""

        if self.batch_size > 1:
            sep = "," if pars.strip() else ""
            pars = pars + "%s batch_size:%d" % (sep, self.batch_size)
        return MESSAGE["init"] + np.int32(rid) + np.int32(len(pars)) + pars.encode()

    def _posdata_message(self, requests):
        """Returns the POSDATA message with the structures of `requests`.

        A batch is padded to `batch_size` structures by repeating the last
        one, so that the client always receives the same number. All the
        structures have the same number of atoms, sent once."""

        if self.batch_size == 1:
            r = requests[0]
            pos = r["pos"][r["active"]]
            return (
                MESSAGE["posdata"]
                + r["cell"][0].tobytes()
                + r["cell"][1].tobytes()
                + np.int32(len(pos) // 3).tobytes()
                + pos.tobytes()
            )
        padded = list(requests) + [requests[-1]] * (self.batch_size - len(requests))
        message = MESSAGE["posdata"] + np.int32(len(padded[0]["active"]) // 3).tobytes()
        for r in padded:
            message += r["cell"][0].tobytes() + r["cell"][1].tobytes()
        for r in padded:
            message += r["pos"][r["active"]].tobytes()
        return message

    def _read_extra(self):
        """Reads an extra string, preceded by its length."""

        elen = int(self.recvall(np.int32()))
        if elen > 0:
            return bytearray(self.recvall(np.zeros(elen, np.dtype("S1")))).decode(
                "utf-8"
            )
        return ""

    def _read_single(self, natoms):
        """Reads the results of one structure, after FORCEREADY."""

        pot = float(self.recvall(np.float64()))
        if int(self.recvall(np.int32())) != natoms:
            raise InvalidSize
        forces = self.recvall(np.zeros(3 * natoms, np.float64))
        vir = self.recvall(np.zeros((3, 3), np.float64))
        return [pot, forces, vir, parse_extra(self._read_extra())]

    def _read_batch(self, natoms, count):
        """Reads the results of a batch, after FORCEREADY, and returns those of
        the first `count` structures (the others are padding)."""

        nbatch = self.batch_size
        pots = self.recvall(np.zeros(nbatch, np.float64))
        if int(self.recvall(np.int32())) != natoms:
            raise InvalidSize
        forces = self.recvall(np.zeros((nbatch, 3 * natoms), np.float64))
        virs = self.recvall(np.zeros((nbatch, 3, 3), np.float64))
        # the extras of the padding structures are read, and discarded
        extras = [self._read_extra() for _ in range(nbatch)]
        return [
            [float(pots[i]), forces[i].copy(), virs[i].copy(), parse_extra(extras[i])]
            for i in range(count)
        ]

    def disconnect(self):
        """Asks the client to exit, and closes the connection. Does nothing
        if the client is already disconnected."""

        if self.status == Status.Disconnected:
            return
        self.status = Status.Disconnected
        if self.exit_on_disconnect:
            # made once the poll that is running here is over
            softexit.request("Client shutdown.")
        try:
            # a best effort, within one time-out: a client that is dropped
            # because it stopped reading must not hold up the others
            socket.socket.sendall(self, MESSAGE["exit"])
            self.shutdown(socket.SHUT_RDWR)
        except OSError:
            pass  # the client is already gone
        self.close()


class SHMDriver(Driver):
    """A client that exchanges positions, cells, energies, forces and virials
    through shared memory.

    The socket carries the headers, the initialization, the names of the
    segments (with the first POSDATA) and the extra strings, and the order of
    the messages on it orders the accesses to the segments. The segments are
    owned by i-PI, named after the process and the descriptor of the socket,
    and sized for `batch_size` structures when the first job is sent.
    """

    _SHM_KINDS = ("pos", "h", "ih", "pot", "force", "vir")

    def __init__(self, sock):
        super(SHMDriver, self).__init__(sock)
        self._first_posdata = True
        self._shm = {}  # kind -> SharedMemory segment
        self._view = {}  # kind -> numpy view, shape (batch_size, ...)
        self._shm_prefix = f"ipi_shm_{os.getpid()}_{self.fileno()}_"

    def _create_segment(self, kind, nbytes):
        name = self._shm_prefix + kind
        try:
            shm = shared_memory.SharedMemory(create=True, size=nbytes, name=name)
        except FileExistsError:
            # left over by a crashed run with the same pid and descriptor
            shared_memory.SharedMemory(name=name).unlink()
            shm = shared_memory.SharedMemory(create=True, size=nbytes, name=name)
        self._shm[kind] = shm
        return shm

    def _alloc_shm(self, natoms):
        n = self.batch_size
        shapes = {
            "pos": (n, 3 * natoms),
            "h": (n, 3, 3),
            "ih": (n, 3, 3),
            "pot": (n,),
            "force": (n, 3 * natoms),
            "vir": (n, 3, 3),
        }
        for kind in self._SHM_KINDS:
            nbytes = int(np.prod(shapes[kind])) * np.dtype(np.float64).itemsize
            shm = self._create_segment(kind, nbytes)
            self._view[kind] = np.ndarray(
                shapes[kind], dtype=np.float64, buffer=shm.buf
            )

    def _free_shm(self):
        self._view.clear()
        for shm in self._shm.values():
            try:
                shm.close()
                shm.unlink()
            except Exception:
                pass
        self._shm.clear()

    def _posdata_message(self, requests):
        """Copies the structures to the segments, and returns the POSDATA
        trigger: header, number of atoms and, the first time, the names of
        the segments. A batch is padded as for sockets."""

        padded = list(requests) + [requests[-1]] * (self.batch_size - len(requests))
        natoms = len(padded[0]["active"]) // 3
        if not self._shm:
            self._alloc_shm(natoms)
        for i, r in enumerate(padded):
            self._view["pos"][i] = r["pos"][r["active"]]
            self._view["h"][i] = r["cell"][0]
            self._view["ih"][i] = r["cell"][1]

        message = MESSAGE["posdata"] + np.int32(natoms).tobytes()
        if self._first_posdata:
            for kind in self._SHM_KINDS:
                name = self._shm[kind].name.encode("utf-8")
                message += np.int32(len(name)).tobytes() + name
            self._first_posdata = False
        return message

    def _shm_results(self, i, extra):
        return [
            float(self._view["pot"][i]),
            self._view["force"][i].copy(),
            self._view["vir"][i].copy(),
            parse_extra(extra),
        ]

    def _read_single(self, natoms):
        return self._shm_results(0, self._read_extra())

    def _read_batch(self, natoms, count):
        # the extras of the padding structures are read, and discarded
        extras = [self._read_extra() for _ in range(self.batch_size)]
        return [self._shm_results(i, extras[i]) for i in range(count)]

    def disconnect(self):
        try:
            super(SHMDriver, self).disconnect()
        finally:
            self._free_shm()

    def __del__(self):
        # also closes the socket, as the finalizer of the socket would
        if hasattr(self, "_shm"):
            self._free_shm()
        super(SHMDriver, self).__del__()


class InterfaceSocket(object):
    """The server socket of a forcefield, and its clients.

    Attributes:
       address: The address of the socket: a host name, or the name of a unix
          socket (in the folder `sockets_prefix`).
       port: The port of an inet socket.
       slots: The minimum length of the queue of clients waiting to connect.
       mode: The type of socket: "unix", "inet" or "shm".
       timeout: The time after which a client that has not returned its
          results is disconnected, and its requests sent to another one. Zero
          for no time-out.
       match_mode: How requests are assigned to clients: "auto" prefers the
          client that computed the same request id at the last step, then new
          clients, then those that are not bound to an id, then any; "any"
          takes any free client; "lock" keeps each id on one client.
       exit_on_disconnect: Whether i-PI exits when a client disconnects.
       batch_size: The number of structures sent together to a client.
       requests: The requests of the forcefield, shared with it.
       offset: An energy subtracted from the results of the clients.
       server: The listening socket, or None when closed.
       clients: The connected clients.
       client_timeout: The time-out of the sockets of the clients, which is
          doubled when a client takes longer to receive the positions.
    """

    def __init__(
        self,
        address="localhost",
        port=31415,
        slots=4,
        mode="unix",
        timeout=1.0,
        match_mode="auto",
        exit_on_disconnect=False,
        sockets_prefix="/tmp/ipi_",
        batch_size=1,
    ):
        self.address = address
        self.port = port
        self.slots = slots
        self.mode = mode
        self.timeout = timeout
        self.sockets_prefix = sockets_prefix
        self.match_mode = match_mode
        self.requests = None  # linked to the requests of the forcefield
        self.exit_on_disconnect = exit_on_disconnect
        self.batch_size = batch_size
        self.offset = 0.0
        self.server = None
        self.clients = []
        self.client_timeout = TIMEOUT

    def open(self):
        """Creates the listening socket."""

        # "shm" uses a unix socket for the messages; only the arrays are moved
        # into shared memory by SHMDriver.
        if self.mode in ("unix", "shm"):
            self.server = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
            try:
                self.server.bind(self.sockets_prefix + self.address)
                info(
                    " @interfacesocket.open: Created unix socket with address "
                    + self.address,
                    verbosity.medium,
                )
            except socket.error:
                self.server.close()
                self.server = None
                raise RuntimeError(
                    "Error opening unix socket. Check if a file "
                    + (self.sockets_prefix + self.address)
                    + " exists, and remove it if unused."
                )

        elif self.mode == "inet":
            self.server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

            try:
                self.server.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
                # TCP_NODELAY is set because Nagle's algorithm slows down a lot
                # the communication pattern of i-PI
                self.server.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
            except OSError as e:
                warning(f"Error setting socket options {e}")

            self.server.bind((self.address, self.port))
            info(
                " @interfacesocket.open: Created inet socket with address "
                + self.address
                + " and port number "
                + str(self.port),
                verbosity.medium,
            )
        else:
            raise NameError(
                "InterfaceSocket mode "
                + self.mode
                + " is not implemented (should be unix/inet/shm)"
            )

        # clients are accepted while i-PI waits for forces: the queue holds
        # those that connect in between
        self.server.listen(max(self.slots, BACKLOG))
        self.server.setblocking(False)
        self.clients = []

    def close(self):
        """Disconnects the clients and closes the listening socket."""

        info(
            " @interfacesocket.close: Shutting down the driver interface.",
            verbosity.low,
        )
        # the clients that answer a status request are not interrupted
        # while they write the reply (which would kill them with SIGPIPE)
        answering = {c.fileno(): c for c in self.clients if c.expect == "status"}
        poller = select.poll()
        for c in answering.values():
            poller.register(c, select.POLLIN)
        deadline = time.time() + self.client_timeout
        while answering and time.time() < deadline:
            # a negative time-out would wait forever
            for fd, _ in poller.poll(max(0.0, deadline - time.time()) * 1e3):
                poller.unregister(fd)
                try:
                    answering.pop(fd).read_status()
                except (Disconnected, OSError):
                    pass
        for c in self.clients:
            c.disconnect()
        self.clients = []
        if self.server is None:
            return
        self.server.close()
        self.server = None
        if self.mode in ("unix", "shm"):
            try:
                os.unlink(self.sockets_prefix + self.address)
            except OSError:
                pass

    def fds(self):
        """Returns the sockets that `poll` acts on once they can be read: the
        listening socket, for new connections, and the clients that owe a
        reply. Free clients send nothing, and are not watched."""

        if self.server is None:
            return []
        return [self.server] + [c for c in self.clients if c.expect is not None]

    def poll(self, idle=False):
        """Makes one round of progress without waiting for any client, and
        returns True if a request was completed or sent.

        The round accepts the new clients and reads the replies that have
        arrived, drops the clients that are late, and sends the queued
        requests to the free clients: an incomplete batch only if no task can
        queue more requests (`idle`)."""

        # closed by a soft exit: the requests will never be computed
        if self.server is None:
            return self._stop_requests()
        progress = self._receive()
        self._check_timeouts()
        progress = self._dispatch(idle) or progress
        # forgets the clients dropped in this round
        self.clients = [c for c in self.clients if c.status != Status.Disconnected]
        return progress

    def _stop_requests(self):
        """Marks the pending requests as stopped, as they cannot be computed
        once the socket is closed, and returns True if there were any."""

        pending = [r for r in self.requests if r["status"] in ("Queued", "Running")]
        for r in pending:
            r["status"] = "Exit"
        return bool(pending)

    def _accept(self):
        """Empties the queue of connections of the listening socket. Each one
        becomes a client, which is sent STATUS: reading the reply completes
        the handshake."""

        while True:
            try:
                client, address = self.server.accept()
            except (BlockingIOError, InterruptedError):
                return  # the queue is empty
            except OSError as err:
                # e.g. a client that reset the connection before it was
                # accepted, or no descriptors left: the others are accepted
                # at a later poll, rather than ending the run
                warning(
                    " @SOCKET:   Could not accept a connection (%s)." % err,
                    verbosity.low,
                )
                return
            client.settimeout(self.client_timeout)
            driver = SHMDriver(client) if self.mode == "shm" else Driver(client)
            driver.batch_size = self.batch_size
            driver.exit_on_disconnect = self.exit_on_disconnect
            info(
                " @interfacesocket.poll: Client asked for connection from "
                + str(address)
                + ". Now hand-shaking.",
                verbosity.low,
            )
            try:
                driver.ask_status()
            except OSError:
                warning(" @SOCKET:   Handshaking failed. Dropping connection.")
                driver.disconnect()
                continue
            self.clients.append(driver)

    def _receive(self):
        """Handles every socket of `fds` that can be read: accepts the new
        connections, and reads the replies of the clients, dropping those
        whose reply fails. Returns True if a job was completed."""

        waiting = {}
        poller = select.poll()
        for c in self.fds():
            waiting[c.fileno()] = c
            poller.register(c, select.POLLIN)
        progress = False
        # a zero time-out only reports the sockets that are ready now
        for fd, _ in poller.poll(0):
            c = waiting[fd]
            if c is self.server:
                self._accept()
                continue
            # readable also when the client has closed: reading then fails
            try:
                if c.expect == "status":
                    self._read_status(c)
                else:
                    self._read_job(c)
                    progress = True
            except (Disconnected, InvalidSize, OSError, UnicodeDecodeError) as err:
                self._drop(c, "failed to reply (%s)" % (str(err) or type(err).__name__))
        return progress

    def _read_status(self, c):
        """Reads the reply to STATUS, which makes the client free if it is
        READY or NEEDINIT. A client with results that i-PI did not ask for
        (HAVEDATA) is dropped."""

        # a client that never had a job is completing its handshake
        known = c.lastreq is not None or c.job is not None
        c.read_status()
        if c.status & Status.HasData:
            self._drop(c, "has results that were not asked for")
        elif not known:
            info(
                " @interfacesocket.poll:   Handshaking was successful. "
                "Added to the client list.",
                verbosity.low,
            )

    def _read_job(self, c):
        """Reads the forces computed by a client, finishes its requests, and
        sends STATUS: the client is free again once it answers."""

        requests = c.job
        results = c.read_job()
        for r, result in zip(requests, results):
            result[0] -= self.offset
            # forces on the inactive atoms are zero
            if len(r["active"]) != len(r["pos"]):
                forces = np.zeros(len(r["pos"]), dtype=np.float64)
                forces[r["active"]] = result[1]
                result[1] = forces
            r.finish(result)
        c.ask_status()

    def _check_timeouts(self):
        """Drops the clients that have computed a job for longer than
        `timeout`, so that other clients can take their requests."""

        if self.timeout <= 0:
            return
        now = time.time()
        for c in self.clients:
            if c.job is not None and now - c.start > self.timeout:
                self._drop(c, "has been computing for %f s" % (now - c.start))

    def _drop(self, c, reason):
        """Disconnects a client, and queues again the requests that it was
        computing. `poll` removes it from the clients at the end of the
        round."""

        warning(
            " @SOCKET:   Client %s %s: disconnecting." % (c.peername, reason),
            verbosity.low,
        )
        # requests that are already done (e.g. the client exited after
        # sending the forces) are left alone
        if c.job is not None:
            for r in c.job:
                if r["status"] == "Running":
                    r["status"] = "Queued"
                    r["start"] = -1
            c.job = None
        c.disconnect()

    def _dispatch(self, idle):
        """Sends the queued requests to the free clients, and returns True if
        any was sent.

        Batches go to the clients in order, and an incomplete batch only if
        `idle`, as other requests may complete it. Single requests are matched to
        clients by passes of decreasing preference, set by `match_mode`:
        "match" pairs a request with the client that computed its id last
        (which may keep state for it, e.g. wavefunctions), "none" with a
        client that has computed nothing yet, "free" with a client that is
        not locked to another id, and "any" with any client."""

        queued = [r for r in self.requests if r["status"] == "Queued"]
        free = [c for c in self.clients if c.free()]
        if not (queued and free):
            return False
        sent = False
        if self.batch_size > 1:
            # batches are taken in order, regardless of the request ids
            if not idle:
                queued = queued[: len(queued) - len(queued) % self.batch_size]
            for c in free:
                if not queued:
                    break
                sent = self._send(c, queued[: self.batch_size]) or sent
                queued = queued[self.batch_size :]
            return sent

        # the passes of each mode, by decreasing preference
        if self.match_mode == "auto":
            match_seq = ["match", "none", "free", "any"]
        elif self.match_mode == "any":
            match_seq = ["any"]
        elif self.match_mode == "lock":
            match_seq = ["match", "none"]
        else:
            raise ValueError("Unknown matching mode " + self.match_mode)
        # in lock mode, requests whose id is bound to a connected client wait for it
        bound = {c.lastreq for c in self.clients if c.status != Status.Disconnected}
        for match in match_seq:
            for c in free[:]:
                for r in queued:
                    if match == "match" and c.lastreq != r["id"]:
                        continue
                    if match == "none" and c.lastreq is not None:
                        continue
                    if (
                        match == "none"
                        and self.match_mode == "lock"
                        and r["id"] in bound
                    ):
                        continue
                    if match == "free" and c.locked:
                        continue
                    # a client that gets its last id again is locked to it
                    c.locked = c.lastreq == r["id"]
                    queued.remove(r)
                    free.remove(c)
                    sent = self._send(c, [r]) or sent
                    break
        return sent

    def _send(self, c, requests):
        """Sends a job to a client, and marks its requests as running. If the
        send fails, the client is dropped and the requests stay queued; if it
        times out, the clients that connect later get a longer time-out."""

        if verbosity.high:
            info(
                " @interfacesocket.poll: %s Assigning %d request(s) with ids %s to "
                "client %s with last id %s"
                % (
                    time.strftime("%y/%m/%d-%H:%M:%S"),
                    len(requests),
                    [r["id"] for r in requests],
                    str(c.peername),
                    str(c.lastreq),
                ),
            )
        try:
            c.send_job(requests)
        except socket.timeout:
            self.client_timeout *= 2
            self._drop(
                c,
                "stopped reading the positions for %d time-outs: the time-out is "
                "%f s for the clients that connect later"
                % (NTIMEOUT + 1, self.client_timeout),
            )
            return False
        except OSError as err:
            self._drop(c, "failed to receive the positions (%s)" % err)
            return False
        now = time.time()
        for r in requests:
            r["status"] = "Running"
            r["t_dispatched"] = now
            r["start"] = now
        return True

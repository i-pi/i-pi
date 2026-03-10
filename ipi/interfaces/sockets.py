"""Deals with the socket communication between the i-PI and drivers.

Deals with creating the socket, transmitting and receiving data, accepting and
removing different driver routines and the parallelization of the force
calculation.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import os
import socket
import select
import time
import threading

import numpy as np
import json

from ipi.utils.messages import verbosity, warning, info
from ipi.utils.softexit import softexit

from concurrent.futures import ThreadPoolExecutor
from ipi.utils.timing_manager import timers
from multiprocessing import shared_memory

__all__ = ["InterfaceSocket"]


HDRLEN = 12
UPDATEFREQ = 10
TIMEOUT = 0.1
SERVERTIMEOUT = 5.0 * TIMEOUT
NTIMEOUT = 20
SELECTTIMEOUT = 60


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
    """Disconnected: Raised if client has been disconnected."""

    pass


class InvalidSize(Exception):
    """Disconnected: Raised if client returns forces with inconsistent number of atoms."""

    pass


class InvalidStatus(Exception):
    """InvalidStatus: Raised if client has the wrong status.

    Shouldn't have to be used if the structure of the program is correct.
    """

    pass


class Status(object):
    """Simple class used to keep track of the status of the client.

    Uses bitwise or to give combinations of different status options.
    i.e. Status.Up | Status.Ready would be understood to mean that the client
    was connected and ready to receive the position and cell data.

    Attributes:
       Disconnected: Flag for if the client has disconnected.
       Up: Flag for if the client is running.
       Ready: Flag for if the client has ready to receive position and cell data.
       NeedsInit: Flag for if the client is ready to receive forcefield
          parameters.
       HasData: Flag for if the client is ready to send force data.
       Busy: Flag for if the client is busy.
       Timeout: Flag for if the connection has timed out.
    """

    Disconnected = 0
    Up = 1
    Ready = 2
    NeedsInit = 4
    HasData = 8
    Busy = 16
    Timeout = 32


class DriverSocket(socket.socket):
    """Deals with communication between the client and driver code.

    Deals with sending and receiving the data between the client and the driver
    code. This class holds common functions which are used in the driver code,
    but can also be used to directly implement a python client.
    Basically it's just a wrapper around socket to simplify some of the
    specific needs of i-PI communication pattern.

    Attributes:
       _buf: A string buffer to hold the reply from the other connection.
    """

    def __init__(self, sock):
        """Initialises DriverSocket.

        Args:
           socket: A socket through which the communication should be done.
        """

        super(DriverSocket, self).__init__(
            sock.family, sock.type, sock.proto, fileno=socket.dup(sock.fileno())
        )
        self.settimeout(sock.gettimeout())

        self._buf = np.zeros(0, np.byte)
        if socket:
            self.peername = self.getpeername()
        else:
            self.peername = "no_socket"

    def send_msg(self, msg):
        """Send the next message through the socket.

        Args:
           msg: The message to send through the socket.
        """
        return self.sendall(MESSAGE[msg])

    def recv_msg(self, length=HDRLEN):
        """Get the next message send through the socket.

        Args:
           l: Length of the accepted message. Defaults to HDRLEN.
        """
        return self.recv(length)

    def recvall(self, dest):
        """Gets the potential energy, force and virial from the driver.

        Args:
           dest: Object to be read into.

        Raises:
           Disconnected: Raised if client is disconnected.

        Returns:
           The data read from the socket to be read into dest.
        """

        blen = dest.itemsize * dest.size
        if blen > len(self._buf):
            self._buf = np.zeros(blen, np.byte)
        bpos = 0
        ntimeout = 0

        while bpos < blen:
            timeout = False

            try:
                bpart = 0
                bpart = self.recv_into(self._buf[bpos:], blen - bpos)
            except socket.timeout:
                # warning(" @SOCKET:   Timeout in recvall, trying again!", verbosity.low)
                timeout = True
                ntimeout += 1
                if ntimeout > NTIMEOUT:
                    warning(
                        " @SOCKET:  Couldn't receive within %5d attempts. Time to give up!"
                        % (NTIMEOUT),
                        verbosity.low,
                    )
                    raise Disconnected()
                pass

            if not timeout and bpart == 0:
                raise Disconnected()
            bpos += bpart

        if dest.ndim > 0:
            dest[:] = np.frombuffer(self._buf[0:blen], dest.dtype).reshape(dest.shape)
            return dest  # tmp.copy() #np.frombuffer(self._buf[0:blen], dest.dtype).reshape(dest.shape).copy()
        else:
            return np.frombuffer(self._buf[0:blen], dest.dtype)[0]

class SHMDriver(DriverSocket):
    """Shared-memory driver client with socket-based control.

    Uses the socket only for control/status messages and to exchange the
    shared-memory buffer names during initialization. Bulk data (positions,
    cell, forces, potential, virial) are transferred through named
    ``multiprocessing.shared_memory`` blocks created on the first dispatch.

    Attributes:
       waitstatus: Boolean giving whether the driver is waiting to get a status answer.
       status: Keeps track of the status of the driver.
       lastreq: The ID of the last request processed by the client.
       locked: Flag to mark if the client has been working consistently on one image.
       first_dispatch: True until shared-memory buffers are allocated.
       nat: Number of atoms inferred from the first request.
       nbeads: Number of beads inferred from the first request.
       *_bufname: Names of shared-memory blocks (pos/h/ih/pot/f/vir).
       *_shm: SharedMemory handles for each buffer.
       *_snp: Numpy views into shared-memory buffers.
    """


    def __init__(self, sock):
        """Initialises an SHM-backed driver client.

        Args:
           socket: Control socket used for status/commands and buffer name exchange.
              Shared-memory blocks are allocated on the first dispatch.
        """

        super(SHMDriver, self).__init__(sock)
        self.waitstatus = False
        self.status = Status.Up
        self.lastreq = None
        self.locked = False
        self.exit_on_disconnect = False
        self.first_dispatch = True

        self.nat = None
        self.nbeads = None

        self.id = sock.fileno()
        # must be uppercase!
        self.pos_bufname = f"IPI-POS-{self.id}"
        self.h_bufname = f"IPI-H-{self.id}"
        self.ih_bufname = f"IPI-IH-{self.id}"
        
        self.pot_bufname = f"IPI-POT-{self.id}"
        self.f_bufname = f"IPI-F-{self.id}"
        self.vir_bufname = f"IPI-VIR-{self.id}"

        info(
            f" @SOCKET: SHMDriver initied with buffer: {self.pos_bufname}",
            verbosity.low,
        )

        self.pos_shm = None
        self.h_shm = None
        self.ih_shm = None
        self.pot_shm = None
        self.f_shm = None
        self.vir_shm = None

        self.pos_snp = None
        self.h_snp = None
        self.ih_snp = None
        self.pot_snp = None
        self.f_snp = None
        self.vir_snp = None


    def shutdown(self, how=socket.SHUT_RDWR):
        """Tries to send an exit message to clients to let them exit gracefully."""

        if self.exit_on_disconnect:
            trd = threading.Thread(
                target=softexit.trigger, kwargs={"message": "Client shutdown."}
            )
            trd.daemon = True
            trd.start()

        self.send_msg("exit")
        self.status = Status.Disconnected

        self.pos_shm.close()
        self.h_shm.close()
        self.ih_shm.close()
        self.pot_shm.close()
        self.f_shm.close()
        self.vir_shm.close()

        self.pos_shm.unlink()
        self.h_shm.unlink()
        self.ih_shm.unlink()
        self.pot_shm.unlink()
        self.f_shm.unlink()
        self.vir_shm.unlink()



        super(DriverSocket, self).shutdown(how)

    def _getstatus_select(self):
        """Gets driver status. Uses socket.select to make sure one can read/write on the socket.

        Returns:
           An integer labelling the status via bitwise or of the relevant members
           of Status.
        """

        if not self.waitstatus:
            try:
                # This can sometimes hang with no timeout.
                readable, writable, errored = select.select(
                    [], [self], [], SELECTTIMEOUT
                )
                if self in writable:
                    self.send_msg("status")
                    self.waitstatus = True
            except socket.error:
                return Status.Disconnected

        try:
            readable, writable, errored = select.select([self], [], [], SELECTTIMEOUT)
            if self in readable:
                reply = self.recv_msg(HDRLEN)
                self.waitstatus = False  # got some kind of reply
            else:
                # This is usually due to VERY slow clients.
                warning(
                    f" @SOCKET: Couldn't find readable socket in {SELECTTIMEOUT}s, will try again",
                    verbosity.low,
                )
                return Status.Busy
        except socket.timeout:
            warning(" @SOCKET:   Timeout in status recv!", verbosity.debug)
            return Status.Up | Status.Busy | Status.Timeout
        except:
            warning(
                " @SOCKET:   Other socket exception. Disconnecting client and trying to carry on.",
                verbosity.debug,
            )
            return Status.Disconnected

        if not len(reply) == HDRLEN:
            return Status.Disconnected
        elif reply == MESSAGE["ready"]:
            return Status.Up | Status.Ready
        elif reply == MESSAGE["needinit"]:
            return Status.Up | Status.NeedsInit
        elif reply == MESSAGE["havedata"]:
            return Status.Up | Status.HasData
        else:
            warning(" @SOCKET:    Unrecognized reply: " + str(reply), verbosity.low)
            return Status.Up

    def _getstatus_direct(self):
        """Gets driver status. Relies on blocking send/recv, which might lead to
        timeouts with slow networks.

        Returns:
           An integer labelling the status via bitwise or of the relevant members
           of Status.
        """

        if not self.waitstatus:
            try:
                self.send_msg("status")
                self.waitstatus = True
            except socket.error:
                return Status.Disconnected
        try:
            reply = self.recv_msg(HDRLEN)
            self.waitstatus = False  # got some kind of reply
        except socket.timeout:
            warning(" @SOCKET:   Timeout in status recv!", verbosity.debug)
            return Status.Up | Status.Busy | Status.Timeout
        except:
            warning(
                " @SOCKET:   Other socket exception. Disconnecting client and trying to carry on.",
                verbosity.debug,
            )
            return Status.Disconnected

        if not len(reply) == HDRLEN:
            return Status.Disconnected
        elif reply == MESSAGE["ready"]:
            return Status.Up | Status.Ready
        elif reply == MESSAGE["needinit"]:
            return Status.Up | Status.NeedsInit
        elif reply == MESSAGE["havedata"]:
            return Status.Up | Status.HasData
        else:
            warning(" @SOCKET:    Unrecognized reply: " + str(reply), verbosity.low)
            return Status.Up

    # depending on the system either _select or _direct can be slightly faster
    # if you're network limited it might be worth experimenting changing this
    _getstatus = _getstatus_select

    def get_status(self):
        """Sets (and returns) the client internal status. Wait for an answer if
        the client is busy."""
        status = self._getstatus()
        while status & Status.Busy:
            status = self._getstatus()
        self.status = status
        return status

    def initialize(self, rid, pars):
        """Sends the initialisation string to the driver.

        Args:
           rid: The index of the request, i.e. the replica that
              the force calculation is for.
           pars: The parameter string to be sent to the driver.

        Raises:
           InvalidStatus: Raised if the status is not NeedsInit.
        """

        if self.status & Status.NeedsInit:
            try:
                # combines all messages in one to reduce latency
                self.sendall(
                    MESSAGE["init"]
                    + np.int32(rid)
                    + np.int32(len(pars))
                    + pars.encode()
                    + np.int32(self.nat)
                    + np.int32(self.nbeads)
                    + Message(self.pos_bufname)
                    + Message(self.h_bufname)
                    + Message(self.ih_bufname)
                    + Message(self.pot_bufname)
                    + Message(self.f_bufname)
                    + Message(self.vir_bufname)
                )
            except:
                self.get_status()
                return
        else:
            raise InvalidStatus("Status in init was " + self.status)

    def sendpos(self, r):
        """Sends the position and cell data to the driver.

        Args:
           r: Request dict containing ``pos`` and ``cell`` data. The data are
              written to shared memory before sending the control header.

        Raises:
           InvalidStatus: Raised if the status is not Ready.
        """
        global TIMEOUT  # we need to update TIMEOUT in case of sendall failure
        # this is to handle the batch mode for multiple bead positions
        timers.start
        
        if self.status & Status.Ready:
            try:
                self.np_to_shm(r) 
                self.sendall(
                    MESSAGE["posdata"]  # header
                )
                self.status = Status.Up | Status.Busy
            except socket.timeout:
                warning(
                    f"Timeout in sendall after {TIMEOUT}s: resetting status and increasing timeout",
                    verbosity.quiet,
                )
                self.status = Status.Timeout
                TIMEOUT *= 2
                return
            except Exception as exc:
                warning(
                    f"Other exception during posdata receive: {exc}", verbosity.quiet
                )
                raise exc
        else:
            raise InvalidStatus("Status in sendpos was " + self.status)

    def getforce(self):
        """Gets the potential energy, force and virial from the driver.

        For SHM drivers the socket is only used to synchronize when the shared
        buffers are ready; data are read from shared memory.

        Raises:
           InvalidStatus: Raised if the status is not HasData.
           Disconnected: Raised if the driver has disconnected.

        Returns:
           A list of the form [potential, force, virial, extra].
        """

        if self.status & Status.HasData:
            self.send_msg("getforce")
            reply = ""
            while True:
                try:
                    reply = self.recv_msg()
                except socket.timeout:
                    warning(
                        " @SOCKET:   Timeout in getforce, trying again!", verbosity.low
                    )
                    continue
                except:
                    warning(
                        " @SOCKET:   Error while receiving message: %s" % (reply),
                        verbosity.low,
                    )
                    raise Disconnected()
                if reply == MESSAGE["forceready"]:
                    break
                else:
                    warning(
                        " @SOCKET:   Unexpected getforce reply: %s" % (reply),
                        verbosity.low,
                    )
                if reply == "":
                    raise Disconnected()
        else:
            raise InvalidStatus("Status in getforce was " + str(self.status))


        return self.shm_to_np()

    
    def alloc_shm(self, r):
        """Allocate and map shared-memory buffers for the first request.

        Initializes shared-memory blocks sized to match the current request
        data layout and creates NumPy views into those buffers.

        Args:
           r: Request dict with ``pos`` and ``cell`` data. ``pos`` may be
              shape (3*nat,) for a single bead or (nbeads, 3*nat).

        Raises:
           Exception: Propagates allocation/mapping errors.
        """
        try:
            if r["pos"].ndim == 1:
                self.nat = len(r["pos"]) // 3
                self.nbeads = 1
            else:
                self.nat = r["pos"].shape[1] // 3
                self.nbeads = r["pos"].shape[0]
            
            # assuming float64 in size*8
            self.pos_shm =  shared_memory.SharedMemory(create=True, size=r["pos"].size*8, name=self.pos_bufname)
            self.h_shm = shared_memory.SharedMemory(create=True, size=9*8, name=self.h_bufname)
            self.ih_shm = shared_memory.SharedMemory(create=True, size=9*8, name=self.ih_bufname)
            
            self.pot_shm = shared_memory.SharedMemory(create=True, size=self.nbeads*8, name=self.pot_bufname)
            self.f_shm = shared_memory.SharedMemory(create=True, size=r["pos"].size*8, name=self.f_bufname)
            self.vir_shm = shared_memory.SharedMemory(create=True, size=self.nbeads*9*8, name=self.vir_bufname)

            self.pos_snp = np.ndarray(r["pos"].shape, dtype=np.float64, buffer=self.pos_shm.buf)
            self.h_snp = np.ndarray((3,3), dtype=np.float64, buffer=self.h_shm.buf)
            self.ih_snp = np.ndarray((3,3), dtype=np.float64, buffer=self.ih_shm.buf)
            
            self.pot_snp =np.ndarray((self.nbeads), dtype=np.float64, buffer=self.pot_shm.buf)
            self.f_snp =np.ndarray(r["pos"].shape, dtype=np.float64, buffer=self.f_shm.buf)
            self.vir_snp =np.ndarray((self.nbeads,3,3), dtype=np.float64, buffer=self.vir_shm.buf)
        except Exception as exc:
            warning(
                f"Exception occured:: {exc}", verbosity.quiet
            )
            raise exc


    def np_to_shm(self, r):
        """Write request data into shared-memory buffers.

        Args:
           r: Request dict with ``pos`` and ``cell`` (h, ih) data.
        """
        self.pos_snp, self.h_snp, self.ih_snp = r["pos"], r["cell"][0], r["cell"][1]

    def shm_to_np(self):
        """Read potential, forces, and virial from shared memory.

        Returns:
           List ``[potential, force, virial, extra]`` matching Driver semantics.
           ``extra`` (``mxtradict``) is currently a placeholder and not yet implemented.
        """
        mxtradict = ""
        if self.nbeads == 1:
            return [float(self.pot_snp), self.f_snp.squeeze(), self.vir_snp.squeeze(), mxtradict]
        else:
            return [self.pot_snp, self.f_snp, self.vir_snp, mxtradict]
        

    def dispatch(self, r):
        """Dispatches a request r and looks after it setting results
        once it has been evaluated. This is meant to be launched as a
        separate thread, and takes care of all the communication related to
        the request.

        On the first dispatch, shared-memory buffers are allocated and their
        names are sent to the driver during initialization.
        """
        
        if self.first_dispatch:
            self.alloc_shm(r)
            self.first_dispatch = False

        #timers.start("[++++++]Get Stat1")
        if not self.status & Status.Up:
            warning(
                " @SOCKET:   Inconsistent client state in dispatch thread! (I)",
                verbosity.low,
            )
            return
        r["t_dispatched"] = time.time()

        self.get_status()
        if self.status & Status.NeedsInit:
            self.initialize(r["id"], r["pars"]) # also sending here buffer names to find shm on driver
            self.status = self.get_status()

        if not (self.status & Status.Ready):
            warning(
                " @SOCKET:   Inconsistent client state in dispatch thread! (II)",
                verbosity.low,
            )
            return

        r["start"] = time.time()

        self.sendpos(r)

        self.get_status()
        if not (self.status & Status.HasData):
            warning(
                " @SOCKET:   Inconsistent client state in dispatch thread! (III)",
                verbosity.low,
            )
            return
        
        try:
            r["result"] = self.getforce()
        except Disconnected:
            self.status = Status.Disconnected
            return
        except Exception as exc:
            warning(
                f"Other exception during force receive: {exc}", verbosity.quiet
            )
            raise exc


        r["result"][0] -= r["offset"]

        if len(r["result"][1]) != len(r["active"]):
            raise InvalidSize

        # If only a piece of the system is active, resize forces and reassign
        if len(r["active"]) != len(r["pos"]):
            rftemp = r["result"][1]
            r["result"][1] = np.zeros(len(r["pos"]), dtype=np.float64)
            r["result"][1] = r["result"][1].at[r["active"]].set(rftemp)
        r["t_finished"] = time.time()
        self.lastreq = r["id"]  #

        # updates the status of the client before leaving
        self.get_status()

        # marks the request as done as the very last thing
        r["status"] = "Done"


class Driver(DriverSocket):
    """Deals with communication between the client and driver code.

    Deals with sending and receiving the data from the driver code. Keeps track
    of the status of the driver. Initialises the driver forcefield, sends the
    position and cell data, and receives the force data.

    Attributes:
       waitstatus: Boolean giving whether the driver is waiting to get a status answer.
       status: Keeps track of the status of the driver.
       lastreq: The ID of the last request processed by the client.
       locked: Flag to mark if the client has been working consistently on one image.
    """

    def __init__(self, sock):
        """Initialises Driver.

        Args:
           socket: A socket through which the communication should be done.
        """

        super(Driver, self).__init__(sock)
        self.waitstatus = False
        self.status = Status.Up
        self.lastreq = None
        self.locked = False
        self.exit_on_disconnect = False

    def shutdown(self, how=socket.SHUT_RDWR):
        """Tries to send an exit message to clients to let them exit gracefully."""

        if self.exit_on_disconnect:
            trd = threading.Thread(
                target=softexit.trigger, kwargs={"message": "Client shutdown."}
            )
            trd.daemon = True
            trd.start()

        self.send_msg("exit")
        self.status = Status.Disconnected

        super(DriverSocket, self).shutdown(how)

    def _getstatus_select(self):
        """Gets driver status. Uses socket.select to make sure one can read/write on the socket.

        Returns:
           An integer labelling the status via bitwise or of the relevant members
           of Status.
        """

        if not self.waitstatus:
            try:
                # This can sometimes hang with no timeout.
                readable, writable, errored = select.select(
                    [], [self], [], SELECTTIMEOUT
                )
                if self in writable:
                    self.send_msg("status")
                    self.waitstatus = True
            except socket.error:
                return Status.Disconnected

        try:
            readable, writable, errored = select.select([self], [], [], SELECTTIMEOUT)
            if self in readable:
                reply = self.recv_msg(HDRLEN)
                self.waitstatus = False  # got some kind of reply
            else:
                # This is usually due to VERY slow clients.
                warning(
                    f" @SOCKET: Couldn't find readable socket in {SELECTTIMEOUT}s, will try again",
                    verbosity.low,
                )
                return Status.Busy
        except socket.timeout:
            warning(" @SOCKET:   Timeout in status recv!", verbosity.debug)
            return Status.Up | Status.Busy | Status.Timeout
        except:
            warning(
                " @SOCKET:   Other socket exception. Disconnecting client and trying to carry on.",
                verbosity.debug,
            )
            return Status.Disconnected

        if not len(reply) == HDRLEN:
            return Status.Disconnected
        elif reply == MESSAGE["ready"]:
            return Status.Up | Status.Ready
        elif reply == MESSAGE["needinit"]:
            return Status.Up | Status.NeedsInit
        elif reply == MESSAGE["havedata"]:
            return Status.Up | Status.HasData
        else:
            warning(" @SOCKET:    Unrecognized reply: " + str(reply), verbosity.low)
            return Status.Up

    def _getstatus_direct(self):
        """Gets driver status. Relies on blocking send/recv, which might lead to
        timeouts with slow networks.

        Returns:
           An integer labelling the status via bitwise or of the relevant members
           of Status.
        """

        if not self.waitstatus:
            try:
                self.send_msg("status")
                self.waitstatus = True
            except socket.error:
                return Status.Disconnected
        try:
            reply = self.recv_msg(HDRLEN)
            self.waitstatus = False  # got some kind of reply
        except socket.timeout:
            warning(" @SOCKET:   Timeout in status recv!", verbosity.debug)
            return Status.Up | Status.Busy | Status.Timeout
        except:
            warning(
                " @SOCKET:   Other socket exception. Disconnecting client and trying to carry on.",
                verbosity.debug,
            )
            return Status.Disconnected

        if not len(reply) == HDRLEN:
            return Status.Disconnected
        elif reply == MESSAGE["ready"]:
            return Status.Up | Status.Ready
        elif reply == MESSAGE["needinit"]:
            return Status.Up | Status.NeedsInit
        elif reply == MESSAGE["havedata"]:
            return Status.Up | Status.HasData
        else:
            warning(" @SOCKET:    Unrecognized reply: " + str(reply), verbosity.low)
            return Status.Up

    # depending on the system either _select or _direct can be slightly faster
    # if you're network limited it might be worth experimenting changing this
    _getstatus = _getstatus_select

    def get_status(self):
        """Sets (and returns) the client internal status. Wait for an answer if
        the client is busy."""
        status = self._getstatus()
        while status & Status.Busy:
            status = self._getstatus()
        self.status = status
        return status

    def initialize(self, rid, pars):
        """Sends the initialisation string to the driver.

        Args:
           rid: The index of the request, i.e. the replica that
              the force calculation is for.
           pars: The parameter string to be sent to the driver.

        Raises:
           InvalidStatus: Raised if the status is not NeedsInit.
        """

        if self.status & Status.NeedsInit:
            try:
                # combines all messages in one to reduce latency
                self.sendall(
                    MESSAGE["init"]
                    + np.int32(rid)
                    + np.int32(len(pars))
                    + pars.encode()
                )
            except:
                self.get_status()
                return
        else:
            raise InvalidStatus("Status in init was " + self.status)

    def sendpos(self, pos, h_ih):
        """Sends the position and cell data to the driver.

        Args:
           pos: An array containing the atom positions.
           cell: A cell object giving the system box.

        Raises:
           InvalidStatus: Raised if the status is not Ready.
        """
        global TIMEOUT  # we need to update TIMEOUT in case of sendall failure
        # this is to handle the batch mode for multiple bead positions
        timers.start
        if pos.ndim == 2:
            nat = pos.shape[1] // 3
            nbeads = pos.shape[0]
        else:
            nat = len(pos) // 3
            nbeads = 1
        if self.status & Status.Ready:
            try:
                # reduces latency by combining all messages in one
                self.sendall(
                    MESSAGE["posdata"]  # header
                    + h_ih[0].tobytes()  # cell
                    + h_ih[1].tobytes()  # inverse cell
                    + np.int32(len(pos) // 3).tobytes()  # length of position array
                    + pos.tobytes()  # positions
                )
                self.status = Status.Up | Status.Busy
            except socket.timeout:
                warning(
                    f"Timeout in sendall after {TIMEOUT}s: resetting status and increasing timeout",
                    verbosity.quiet,
                )
                self.status = Status.Timeout
                TIMEOUT *= 2
                return
            except Exception as exc:
                warning(
                    f"Other exception during posdata receive: {exc}", verbosity.quiet
                )
                raise exc
        else:
            raise InvalidStatus("Status in sendpos was " + self.status)

    def getforce(self):
        """Gets the potential energy, force and virial from the driver.

        Raises:
           InvalidStatus: Raised if the status is not HasData.
           Disconnected: Raised if the driver has disconnected.

        Returns:
           A list of the form [potential, force, virial, extra].
        """

        if self.status & Status.HasData:
            self.send_msg("getforce")
            reply = ""
            while True:
                try:
                    reply = self.recv_msg()
                except socket.timeout:
                    warning(
                        " @SOCKET:   Timeout in getforce, trying again!", verbosity.low
                    )
                    continue
                except:
                    warning(
                        " @SOCKET:   Error while receiving message: %s" % (reply),
                        verbosity.low,
                    )
                    raise Disconnected()
                if reply == MESSAGE["forceready"]:
                    break
                else:
                    warning(
                        " @SOCKET:   Unexpected getforce reply: %s" % (reply),
                        verbosity.low,
                    )
                if reply == "":
                    raise Disconnected()
        else:
            raise InvalidStatus("Status in getforce was " + str(self.status))

        mu = np.float64()
        mu = self.recvall(mu)

        mlen = np.int32()
        mlen = self.recvall(mlen)
        """
        if self.nbeads > 1:
            mf = np.zeros((mbeads, 3 * mlen), np.float64)
            mf = self.recvall(mf)

            mvir = np.zeros((mbeads, 3, 3), np.float64)
            mvir = self.recvall(mvir)
        else:
            mf = np.zeros(3 * mlen, np.float64)
            mf = self.recvall(mf)

            mvir = np.zeros((3, 3), np.float64)
            mvir = self.recvall(mvir)
        """
        mf = np.zeros(3 * mlen, np.float64)
        mf = self.recvall(mf)

        mvir = np.zeros((3, 3), np.float64)
        mvir = self.recvall(mvir)

        # Machinery to return a string as an "extra" field.
        # Comment if you are using a ancient patched driver that does not return anything!
        # Actually, you should really update your driver, you're like a decade behind.
        mlen = np.int32()
        mlen = self.recvall(mlen)
        if mlen > 0:
            mxtra = np.zeros(mlen, np.dtype("S1"))
            mxtra = self.recvall(mxtra)
            mxtra = bytearray(mxtra).decode("utf-8")
        else:
            mxtra = ""
        mxtradict = {}
        if mxtra:
            try:
                mxtradict = json.loads(mxtra)
                info(
                    "@driver.getforce: Extra string JSON has been loaded.",
                    verbosity.debug,
                )
            except:
                # if we can't parse it as a dict, issue a warning and carry on
                info(
                    "@driver.getforce: Extra string could not be loaded as a dictionary. Extra="
                    + mxtra,
                    verbosity.debug,
                )
                mxtradict = {}
                pass
            if "raw" in mxtradict:
                raise ValueError(
                    "'raw' cannot be used as a field in a JSON-formatted extra string"
                )

            mxtradict["raw"] = mxtra
        return [mu, mf, mvir, mxtradict]

    def dispatch(self, r):
        """Dispatches a request r and looks after it setting results
        once it has been evaluated. This is meant to be launched as a
        separate thread, and takes care of all the communication related to
        the request.
        """

        if not self.status & Status.Up:
            warning(
                " @SOCKET:   Inconsistent client state in dispatch thread! (I)",
                verbosity.low,
            )
            return

        r["t_dispatched"] = time.time()

        self.get_status()
        if self.status & Status.NeedsInit:
            self.initialize(r["id"], r["pars"])
            self.status = self.get_status()

        if not (self.status & Status.Ready):
            warning(
                " @SOCKET:   Inconsistent client state in dispatch thread! (II)",
                verbosity.low,
            )
            return

        r["start"] = time.time()
        self.sendpos(r["pos"][r["active"]], r["cell"])

        self.get_status()
        if not (self.status & Status.HasData):
            warning(
                " @SOCKET:   Inconsistent client state in dispatch thread! (III)",
                verbosity.low,
            )
            return

        try:
            r["result"] = self.getforce()
        except Disconnected:
            self.status = Status.Disconnected
            return

        r["result"][0] -= r["offset"]

        if len(r["result"][1]) != len(r["pos"][r["active"]]):
            raise InvalidSize

        # If only a piece of the system is active, resize forces and reassign
        if len(r["active"]) != len(r["pos"]):
            rftemp = r["result"][1]
            r["result"][1] = np.zeros(len(r["pos"]), dtype=np.float64)
            r["result"][1][r["active"]] = rftemp
        r["t_finished"] = time.time()
        self.lastreq = r["id"]  #

        # updates the status of the client before leaving
        self.get_status()

        # marks the request as done as the very last thing
        r["status"] = "Done"


class InterfaceSocket(object):
    """Host server class.

    Deals with distribution of all the jobs between the different client servers
    and both initially and as clients either finish or are disconnected.
    Deals with cleaning up after all calculations are done. Also deals with the
    threading mechanism, and cleaning up if the interface is killed.

    Attributes:
       address: A string giving the name of the host network.
       port: An integer giving the port the socket will be using.
       slots: An integer giving the maximum allowed backlog of queued clients.
       mode: A string giving the type of socket used.
       latency: A float giving the number of seconds the interface will wait
          before updating the client list.
       timeout: A float giving a timeout limit for considering a calculation dead
          and dropping the connection.
       server: The socket used for data transmission.
       clients: A list of the driver clients connected to the server.
       requests: A list of all the jobs required in the current PIMD step.
       jobs: A list of all the jobs currently running.
       _poll_thread: The thread the poll loop is running on.
       _prev_kill: Holds the signals to be sent to clean up the main thread
          when a kill signal is sent.
       _poll_true: A boolean giving whether the thread is alive.
       _poll_iter: An integer used to decide whether or not to check for
          client connections. It is used as a counter, once it becomes higher
          than the pre-defined number of steps between checks the socket will
          update the list of clients and then be reset to zero.
    """

    def __init__(
        self,
        address="localhost",
        port=31415,
        slots=4,
        mode="unix",
        shm=False,
        timeout=1.0,
        match_mode="auto",
        exit_on_disconnect=False,
        max_workers=128,
        sockets_prefix="/tmp/ipi_",
    ):
        """Initialises interface.

        Args:
           address: An optional string giving the name of the host server.
              Defaults to 'localhost'.
           port: An optional integer giving the port number. Defaults to 31415.
           slots: An optional integer giving the maximum allowed backlog of
              queueing clients. Defaults to 4.
           mode: An optional string giving the type of socket. Defaults to 'unix'.
           timeout: Length of time waiting for data from a client before we assume
              the connection is dead and disconnect the client.
            max_workers: Maximum number of threads launched concurrently

        Raises:
           NameError: Raised if mode is not 'unix' or 'inet'.
        """

        self.address = address
        self.port = port
        self.slots = slots
        self.mode = mode
        self.shm = shm
        self.timeout = timeout
        self.sockets_prefix = sockets_prefix
        self.poll_iter = UPDATEFREQ  # triggers pool_update at first poll
        self.prlist = []  # list of pending requests
        self.match_mode = match_mode  # heuristics to match jobs and active clients
        self.requests = None  # these will be linked to the request list of the FFSocket object using the interface
        self.exit_on_disconnect = exit_on_disconnect
        self.max_workers = max_workers
        self.offset = 0.0  # a constant energy offset added to the results returned by the driver (hacky but simple)

    def open(self):
        """Creates a new socket.

        Used so that we can create a interface object without having to also
        create the associated socket object.
        """

        if self.mode == "unix":
            self.server = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
            try:
                self.server.bind(self.sockets_prefix + self.address)
                info(
                    " @interfacesocket.open: Created unix socket with address "
                    + self.address,
                    verbosity.medium,
                )
            except socket.error:
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
                + " is not implemented (should be unix/inet)"
            )

        self.server.listen(self.slots)
        self.server.settimeout(SERVERTIMEOUT)

        # these are the two main objects the socket interface should worry about and manage
        self.clients = []  # list of active clients (working or ready to compute)
        self.jobs = []  # list of jobs
        self.executor = ThreadPoolExecutor(max_workers=self.max_workers)

    def close(self):
        """Closes down the socket."""

        info(
            " @interfacesocket.close: Shutting down the driver interface.",
            verbosity.low,
        )

        for c in self.clients:
            try:
                c.shutdown(socket.SHUT_RDWR)
                c.close()
            except:
                pass

        # flush it all down the drain
        self.clients = []
        self.jobs = []

        try:
            self.server.shutdown(socket.SHUT_RDWR)
            self.server.close()
        except:
            info(
                " @interfacesocket.close: Problem shutting down the server socket. Will just continue and hope for the best.",
                verbosity.low,
            )
        if self.mode == "unix":
            os.unlink(self.sockets_prefix + self.address)

    def poll(self):
        """Called in the main thread loop.

        Runs until either the program finishes or a kill call is sent. Updates
        the pool of clients every UPDATEFREQ loops and loops every latency seconds.
        The actual loop is in the associated forcefield class.
        """

        # makes sure to remove the last dead client as soon as possible -- and to get clients if we are dry
        if (
            self.poll_iter >= UPDATEFREQ
            or len(self.clients) == 0
            or (len(self.clients) > 0 and not (self.clients[0].status & Status.Up))
        ):
            self.poll_iter = 0
            self.pool_update()

        self.poll_iter += 1
        self.pool_distribute()

    def pool_update(self):
        """Deals with keeping the pool of client drivers up-to-date during a
        force calculation step.

        Deals with maintaining the client list. Clients that have
        disconnected are removed and their jobs removed from the list of
        running jobs and new clients are connected to the server.
        """

        # check for disconnected clients
        for c in self.clients[:]:
            if not (c.status & Status.Up):
                try:
                    warning(
                        " @SOCKET:   Client "
                        + str(c.peername)
                        + " died or got unresponsive(C). Removing from the list.",
                        verbosity.low,
                    )
                    c.shutdown(socket.SHUT_RDWR)
                    c.close()
                except socket.error:
                    pass
                c.status = Status.Disconnected
                self.clients.remove(c)
                # requeue jobs that have been left hanging
                for [k, j, tc] in self.jobs[:]:
                    tc.result()
                    if j is c:
                        self.jobs = [
                            w for w in self.jobs if not (w[0] is k and w[1] is j)
                        ]  # removes pair in a robust way

                        k["status"] = "Queued"
                        k["start"] = -1

        if len(self.clients) == 0:
            searchtimeout = SERVERTIMEOUT
        else:
            searchtimeout = 0.0

        keepsearch = True
        while keepsearch:
            readable, writable, errored = select.select(
                [self.server], [], [], searchtimeout
            )
            if self.server in readable:
                client, address = self.server.accept()
                client.settimeout(TIMEOUT)
               
                if self.shm:
                    driver = SHMDriver(client)
                    info(" @interfacesocket.pool_update: Using SHM communication", verbosity.low)
                else:
                    driver = Driver(client)
                
                info(
                    " @interfacesocket.pool_update:   Client asked for connection from "
                    + str(address)
                    + ". Now hand-shaking.",
                    verbosity.low,
                )
                driver.get_status()
                if driver.status | Status.Up:
                    driver.exit_on_disconnect = self.exit_on_disconnect
                    self.clients.append(driver)
                    info(
                        " @interfacesocket.pool_update:   Handshaking was successful. Added to the client list.",
                        verbosity.low,
                    )
                    self.poll_iter = UPDATEFREQ  # if a new client was found, will try again harder next time
                    searchtimeout = SERVERTIMEOUT
                else:
                    warning(
                        " @SOCKET:   Handshaking failed. Dropping connection.",
                        verbosity.low,
                    )
                    client.shutdown(socket.SHUT_RDWR)
                    client.close()
            else:
                keepsearch = False

    def pool_distribute(self):
        """Deals with keeping the list of jobs up-to-date during a force
        calculation step.

        Deals with maintaining the jobs list. Gets data from drivers that have
        finished their calculation and removes that job from the list of running
        jobs, adds jobs to free clients and initialises the forcefields of new
        clients.
        """

        ttotal = tdispatch = tcheck = 0
        ttotal -= time.time()

        # get clients that are still free
        freec = self.clients[:]
        for [r2, c, ct] in self.jobs:
            freec.remove(c)

        # fills up list of pending requests if empty, or if clients are abundant
        if len(self.prlist) == 0 or len(freec) > len(self.prlist):
            self.prlist = [r for r in self.requests if r["status"] == "Queued"]

        if self.match_mode == "auto":
            match_seq = ["match", "none", "free", "any"]
        elif self.match_mode == "any":
            match_seq = ["any"]
        elif self.match_mode == "lock":
            match_seq = ["match", "none"]

        # first: dispatches jobs to free clients (if any!)
        # tries first to match previous replica<>driver association, then to get new clients, and only finally send the a new replica to old drivers
        ndispatch = 0
        tdispatch -= time.time()
        while len(freec) > 0 and len(self.prlist) > 0:
            for match_ids in match_seq:
                for fc in freec[:]:
                    if self.dispatch_free_client(fc, match_ids):
                        freec.remove(fc)
                        ndispatch += 1

                    if len(self.prlist) == 0:
                        break
            # if using lock mode, check that there is a at least one client-replica match in the lists freec and prlist.
            # If not, we break out of the while loop
            if self.match_mode == "lock":
                break

            if len(freec) > 0:
                self.prlist = [r for r in self.requests if r["status"] == "Queued"]
        tdispatch += time.time()

        # now check for client status
        if len(self.jobs) == 0:
            for c in self.clients:
                if (
                    c.status == Status.Disconnected
                ):  # client disconnected. force a pool_update
                    self.poll_iter = UPDATEFREQ
                    return

        # check for finished jobs
        nchecked = 0
        nfinished = 0
        tcheck -= time.time()
        for [r, c, ct] in self.jobs[:]:
            chk = self.check_job_finished(r, c, ct)
            if chk == 1:
                nfinished += 1
            elif chk == 0:
                self.poll_iter = UPDATEFREQ  # client disconnected. force a pool_update
            nchecked += 1
        tcheck += time.time()

        ttotal += time.time()
        # info("POLL TOTAL: %10.4f  Dispatch(N,t):  %4i, %10.4f   Check(N,t):   %4i, %10.4f" % (ttotal, ndispatch, tdispatch, nchecked, tcheck), verbosity.debug)

        if nfinished > 0:
            # don't wait, just try again to distribute
            self.pool_distribute()

    def dispatch_free_client(self, fc, match_ids="any", send_threads=[]):
        """
        Tries to find a request to match a free client.
        """

        # first, makes sure that the client is REALLY free
        if not (fc.status & Status.Up):
            return False
        if fc.status & Status.HasData:
            return False
        if not (fc.status & (Status.Ready | Status.NeedsInit | Status.Busy)):
            warning(
                " @SOCKET: Client "
                + str(fc.peername)
                + " is in an unexpected status "
                + str(fc.status)
                + " at (1). Will try to keep calm and carry on.",
                verbosity.low,
            )
            return False

        for r in self.prlist[:]:
            if match_ids == "match" and fc.lastreq is not r["id"]:
                continue
            elif match_ids == "none" and fc.lastreq is not None:
                continue
            elif (
                self.match_mode == "lock"
                and match_ids == "none"
                and (r["id"] in [c.lastreq for c in self.clients])
            ):
                # if using lock mode and the user connects more clients than there are replicas, do not allow this client to
                # be matched with a pending request.
                continue

            elif match_ids == "free" and fc.locked:
                continue

            # makes sure the request is marked as running and the client included in the jobs list
            fc.locked = fc.lastreq is r["id"]

            r["offset"] = (
                self.offset
            )  # transmits with the request an offset value for the energy (typically zero)

            r["status"] = "Running"
            self.prlist.remove(r)
            info(
                " @interfacesocket.dispatch_free_client: %s Assigning [%5s] request id %4s to client with last-id %4s (% 3d/% 3d : %s)"
                % (
                    time.strftime("%y/%m/%d-%H:%M:%S"),
                    match_ids,
                    str(r["id"]),
                    str(fc.lastreq),
                    self.clients.index(fc),
                    len(self.clients),
                    str(fc.peername),
                ),
                verbosity.high,
            )
            # fc_thread = threading.Thread(
            #    target=fc.dispatch, name="DISPATCH", kwargs={"r": r}
            # )
            fc_thread = self.executor.submit(fc.dispatch, r=r)
            self.jobs.append([r, fc, fc_thread])
            # fc_thread.daemon = True
            # fc_thread.start()
            return True

        return False

    def check_job_finished(self, r, c, ct):
        """
        Checks if a job has been completed, and retrieves the results
        """

        if r["status"] == "Done":
            ct.result()
            self.jobs = [
                w for w in self.jobs if not (w[0] is r and w[1] is c)
            ]  # removes pair in a robust way
            return 1

        if (
            self.timeout > 0
            and r["start"] > 0
            and time.time() - r["start"] > self.timeout
        ):
            warning(
                " @SOCKET:  Timeout! request has been running for "
                + str(time.time() - r["start"])
                + " sec.",
                verbosity.low,
            )
            warning(
                " @SOCKET:   Client "
                + str(c.peername)
                + " died or got unresponsive(A). Disconnecting.",
                verbosity.low,
            )
            try:
                c.shutdown(socket.SHUT_RDWR)
            except socket.error:
                pass
            c.close()
            c.status = Status.Disconnected
            return 0  # client will be cleared and request resuscitated in poll_update

        return -1

"""Deals with the socket communication between the i-PI and drivers.

Deals with creating the socket, transmitting and receiving data, accepting and
removing different driver routines and the parallelization of the force
calculation.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import sys
import os
import socket
import select
import string
import time
import threading

import numpy as np

from ipi.utils.depend import dstrip
from ipi.utils.messages import verbosity, warning, info


__all__ = ['InterfaceSocket']


HDRLEN = 12
UPDATEFREQ = 10
TIMEOUT = 0.02
SERVERTIMEOUT = 5.0 * TIMEOUT
NTIMEOUT = 20


def Message(mystr):
    """Returns a header of standard length HDRLEN."""

    return string.ljust(string.upper(mystr), HDRLEN)


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

    Attributes:
       _buf: A string buffer to hold the reply from the other connection.
    """

    def __init__(self, socket):
        """Initialises DriverSocket.

        Args:
           socket: A socket through which the communication should be done.
        """

        super(DriverSocket, self).__init__(_sock=socket)
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
        return self.sendall(Message(msg))

    def recv_msg(self, l=HDRLEN):
        """Get the next message send through the socket.

        Args:
           l: Length of the accepted message. Defaults to HDRLEN.
        """
        return self.recv(l)

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
        if (blen > len(self._buf)):
            self._buf.resize(blen)
        bpos = 0
        ntimeout = 0

        while bpos < blen:
            timeout = False

            # pre-2.5 version.
            try:
                bpart = ""
                bpart = self.recv(blen - bpos)
                if len(bpart) == 0:
                    raise socket.timeoout    # if this keeps returning no data, we are in trouble....
                self._buf[bpos:bpos + len(bpart)] = np.fromstring(bpart, np.byte)
            except socket.timeout:
                #warning(" @SOCKET:   Timeout in recvall, trying again!", verbosity.low)
                timeout = True
                ntimeout += 1
                if ntimeout > NTIMEOUT:
                    warning(" @SOCKET:  Couldn't receive within %5d attempts. Time to give up!" % (NTIMEOUT), verbosity.low)
                    raise Disconnected()
                pass
            if not timeout and len(bpart) == 0:
                raise Disconnected()
            bpos += len(bpart)

            # post-2.5 version: slightly more compact for modern python versions
            # try:
            #   bpart = 1
            #   bpart = self.recv_into(self._buf[bpos:], blen-bpos)
            # except socket.timeout:
            #   print " @SOCKET:   Timeout in status recvall, trying again!"
            #   timeout = True
            #   pass
            # if (not timeout and bpart == 0):
            #   raise Disconnected()
            # bpos += bpart
            # TODO this Disconnected() exception currently just causes the program to hang.
            # This should do something more graceful

        if np.isscalar(dest):
            return np.fromstring(self._buf[0:blen], dest.dtype)[0]
        else:
            return np.fromstring(self._buf[0:blen], dest.dtype).reshape(dest.shape)


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

    def __init__(self, socket):
        """Initialises Driver.

        Args:
           socket: A socket through which the communication should be done.
        """

        super(Driver, self).__init__(socket=socket)
        self.waitstatus = False
        self.status = Status.Up
        self.lastreq = None
        self.locked = False
        self.threadlock = threading.Lock()

    def shutdown(self, how=socket.SHUT_RDWR):
        """Tries to send an exit message to clients to let them exit gracefully."""

        self.sendall(Message("exit"))
        self.status = Status.Disconnected
        super(DriverSocket, self).shutdown(how)

    def poll(self):
        """Waits for driver status."""

        with self.threadlock:
            self.status = Status.Disconnected  # sets disconnected as failsafe status, in case _getstatus fails and exceptions are ignored upstream
            self.status = self._getstatus()

    def _getstatus(self):
        """Gets driver status.

        Returns:
           An integer labelling the status via bitwise or of the relevant members
           of Status.
        """

        if not self.waitstatus:
            try:
                # This can sometimes hang with no timeout.
                # Using the recommended 60 s.
                readable, writable, errored = select.select([], [self], [], 60)
                if self in writable:
                    self.sendall(Message("status"))
                    self.waitstatus = True
            except socket.error:
                return Status.Disconnected

        try:
            reply = self.recv(HDRLEN)
            self.waitstatus = False  # got some kind of reply
        except socket.timeout:
            warning(" @SOCKET:   Timeout in status recv!", verbosity.trace)
            return Status.Up | Status.Busy | Status.Timeout
        except:
            return Status.Disconnected

        if not len(reply) == HDRLEN:
            return Status.Disconnected
        elif reply == Message("ready"):
            return Status.Up | Status.Ready
        elif reply == Message("needinit"):
            return Status.Up | Status.NeedsInit
        elif reply == Message("havedata"):
            return Status.Up | Status.HasData
        else:
            warning(" @SOCKET:    Unrecognized reply: " + str(reply), verbosity.low)
            return Status.Up

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
                self.sendall(Message("init"))
                self.sendall(np.int32(rid))
                self.sendall(np.int32(len(pars)))
                self.sendall(pars)
            except:
                self.poll()
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

        if (self.status & Status.Ready):
            try:
                self.sendall(Message("posdata"))
                self.sendall(h_ih[0])
                self.sendall(h_ih[1])
                self.sendall(np.int32(len(pos) / 3))
                self.sendall(pos)
                self.status = Status.Up | Status.Busy
            except:
                self.poll()
                return
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

        with self.threadlock:
            if (self.status & Status.HasData):
                self.sendall(Message("getforce"));
                reply = ""
                while True:
                    try:
                        reply = self.recv_msg()
                    except socket.timeout:
                        warning(" @SOCKET:   Timeout in getforce, trying again!", verbosity.low)
                        continue
                    if reply == Message("forceready"):
                        break
                    else:
                        warning(" @SOCKET:   Unexpected getforce reply: %s" % (reply), verbosity.low)
                    if reply == "":
                        raise Disconnected()
            else:
                raise InvalidStatus("Status in getforce was " + str(self.status))

            mu = np.float64()
            mu = self.recvall(mu)

            mlen = np.int32()
            mlen = self.recvall(mlen)
            mf = np.zeros(3 * mlen, np.float64)
            mf = self.recvall(mf)

            mvir = np.zeros((3, 3), np.float64)
            mvir = self.recvall(mvir)

            #! Machinery to return a string as an "extra" field. Comment if you are using a old patched driver that does not return anything!
            mlen = np.int32()
            mlen = self.recvall(mlen)
            if mlen > 0:
                mxtra = np.zeros(mlen, np.character)
                mxtra = self.recvall(mxtra)
                mxtra = "".join(mxtra)
            else:
                mxtra = ""

        return [mu, mf, mvir, mxtra]


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
       server: The socket used for data transmition.
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

    def __init__(self, address="localhost", port=31415, slots=4, mode="unix", timeout=1.0, match_mode="auto"):
        """Initialises interface.

        Args:
           address: An optional string giving the name of the host server.
              Defaults to 'localhost'.
           port: An optional integer giving the port number. Defaults to 31415.
           slots: An optional integer giving the maximum allowed backlog of
              queueing clients. Defaults to 4.
           mode: An optional string giving the type of socket. Defaults to 'unix'.
           latency: An optional float giving the time in seconds the socket will
              wait before updating the client list. Defaults to 1e-3.
           timeout: Length of time waiting for data from a client before we assume
              the connection is dead and disconnect the client.

        Raises:
           NameError: Raised if mode is not 'unix' or 'inet'.
        """

        self.address = address
        self.port = port
        self.slots = slots
        self.mode = mode
        self.timeout = timeout
        self.poll_iter = UPDATEFREQ  # triggers pool_update at first poll
        self.prlist = []
        self.match_mode = match_mode
        self.threadlock = threading.Lock()

    def open(self):
        """Creates a new socket.

        Used so that we can create a interface object without having to also
        create the associated socket object.
        """

        if self.mode == "unix":
            self.server = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
            try:
                self.server.bind("/tmp/ipi_" + self.address)
                info("Created unix socket with address " + self.address, verbosity.medium)
            except socket.error:
                raise RuntimeError("Error opening unix socket. Check if a file " + ("/tmp/ipi_" + self.address) + " exists, and remove it if unused.")

        elif self.mode == "inet":
            self.server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            self.server.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            self.server.bind((self.address, self.port))
            info("Created inet socket with address " + self.address + " and port number " + str(self.port), verbosity.medium)
        else:
            raise NameError("InterfaceSocket mode " + self.mode + " is not implemented (should be unix/inet)")

        self.server.listen(self.slots)
        self.server.settimeout(SERVERTIMEOUT)
        self.clients = []
        self.jobs = []

    def close(self):
        """Closes down the socket."""

        info(" @SOCKET: Shutting down the driver interface.", verbosity.low)

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
            info(" @SOCKET: Problem shutting down the server socket. Will just continue and hope for the best.", verbosity.low)
        if self.mode == "unix":
            os.unlink("/tmp/ipi_" + self.address)

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
                    warning(" @SOCKET:   Client " + str(c.peername) + " died or got unresponsive(C). Removing from the list.", verbosity.low)
                    c.shutdown(socket.SHUT_RDWR)
                    c.close()
                except socket.error:
                    pass
                c.status = Status.Disconnected
                self.clients.remove(c)
                # requeue jobs that have been left hanging
                with self.threadlock:
                    for [k, j] in self.jobs[:]:
                        if j is c:
                            self.jobs = [w for w in self.jobs if not (w[0] is k and w[1] is j)]  # removes pair in a robust way

                            k["status"] = "Queued"
                            k["start"] = -1

        if len(self.clients) == 0:
            searchtimeout = SERVERTIMEOUT
        else:
            searchtimeout = 0.0

        keepsearch = True
        while keepsearch:
            readable, writable, errored = select.select([self.server], [], [], searchtimeout)
            if self.server in readable:
                client, address = self.server.accept()
                client.settimeout(TIMEOUT)
                driver = Driver(client)
                info(" @SOCKET:   Client asked for connection from " + str(address) + ". Now hand-shaking.", verbosity.low)
                driver.poll()
                if (driver.status | Status.Up):
                    self.clients.append(driver)
                    info(" @SOCKET:   Handshaking was successful. Added to the client list.", verbosity.low)
                    self.poll_iter = UPDATEFREQ   # if a new client was found, will try again harder next time
                    searchtimeout = SERVERTIMEOUT
                else:
                    warning(" @SOCKET:   Handshaking failed. Dropping connection.", verbosity.low)
                    client.shutdown(socket.SHUT_RDWR)
                    client.close()
            else:
                keepsearch = False

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
            warning(" @SOCKET: Client " + str(fc.peername) + " is in an unexpected status " + str(fc.status) + " at (1). Will try to keep calm and carry on.", verbosity.low)
            return False

        for r in self.prlist[:]:
            if match_ids == "match" and not fc.lastreq is r["id"]:
                continue
            elif match_ids == "none" and not fc.lastreq is None:
                continue
            elif match_ids == "free" and fc.locked:
                continue
            info(" @SOCKET: %s Assigning [%5s] request id %4s to client with last-id %4s (% 3d/% 3d : %s)" % (time.strftime("%y/%m/%d-%H:%M:%S"), match_ids, str(r["id"]), str(fc.lastreq), self.clients.index(fc), len(self.clients), str(fc.peername)), verbosity.high)

            while fc.status & Status.Busy:
                fc.poll()

            if fc.status & Status.NeedsInit:
                fc.initialize(r["id"], r["pars"])
                fc.poll()
                while fc.status & Status.Busy:  # waits for initialization to finish. hopefully this is fast
                    fc.poll()

            if fc.status & Status.Ready:
                r["status"] = "Running"
                #fc.sendpos(r["pos"][r["active"]], r["cell"])
                nt = threading.Thread(target=fc.sendpos, name="sendpos", args = (r["pos"][r["active"]], r["cell"]))
                nt.daemon = True
                nt.start()
                send_threads.append(nt)

                r["t_dispatched"] = time.time()
                r["start"] = time.time()  # sets start time for the request
                # fc.poll()
                with self.threadlock:
                    self.jobs.append([r, fc])
                fc.locked = (fc.lastreq is r["id"])

                # removes r from the list of pending jobs
                self.prlist.remove(r)

                return True
            else:
                warning(" @SOCKET: Client " + str(fc.peername) + " is in an unexpected status " + str(fc.status) + " at (2). Will try to keep calm and carry on.", verbosity.low)

    def check_job_finished(self, r, c):
        """
            Checks if a job has been completed, and retrieves the results
        """

        if not c.status & (Status.Ready | Status.NeedsInit):
            c.poll()
        if not (c.status & Status.HasData):
            return False
        try:
            r["result"] = c.getforce()
            if len(r["result"][1]) != len(r["pos"][r["active"]]):
                raise InvalidSize
            # If only a piece of the system is active, resize forces and reassign
            rftemp = r["result"][1]
            r["result"][1] = np.zeros(len(r["pos"]), dtype=np.float64)
            r["result"][1][r["active"]] = rftemp
        except Disconnected:
            c.status = Status.Disconnected
            return False
        except InvalidSize:
            warning(" @SOCKET:   Client returned an inconsistent number of forces. Will mark as disconnected and try to carry on.", verbosity.low)
            c.status = Status.Disconnected
            return False
        except:
            warning(" @SOCKET:   Client got in a awkward state during getforce. Will mark as disconnected and try to carry on.", verbosity.low)
            c.status = Status.Disconnected
            return False

        c.poll()
        while c.status & Status.Busy:  # waits, but check if we got stuck.
            if self.timeout > 0 and r["start"] > 0 and time.time() - r["start"] > self.timeout:
                warning(" @SOCKET:  Timeout! HASDATA for bead " + str(r["id"]) + " has been running for " + str(time.time() - r["start"]) + " sec.", verbosity.low)
                warning(" @SOCKET:   Client " + str(c.peername) + " died or got unresponsive(A). Disconnecting.", verbosity.low)
                try:
                    c.shutdown(socket.SHUT_RDWR)
                except socket.error:
                    pass
                c.close()
                c.status = Status.Disconnected
                return False
            c.poll()
        if not (c.status & Status.Up):
            warning(" @SOCKET:   Client died a horrible death while getting forces. Will try to cleanup.", verbosity.low)
            return False
        r["status"] = "Done"
        r["t_finished"] = time.time()
        c.lastreq = r["id"]  # saves the ID of the request that the client has just processed
        return True

    def join_threads(self, sthreads, cthreads):
        for st in sthreads:
            while st.isAlive():
                st.join(1.0)

        for ct, r, c in cthreads:
            while ct is not None and ct.isAlive():
                ct.join(1.0)

            if r["status"] == "Done":
                with self.threadlock:
                    self.jobs = [w for w in self.jobs if not (w[0] is r and w[1] is c)]  # removes pair in a robust way

            if self.timeout > 0 and c.status != Status.Disconnected and r["start"] > 0 and time.time() - r["start"] > self.timeout:
                warning(" @SOCKET:  Timeout! Request for bead " + str(r["id"]) + " has been running for " + str(time.time() - r["start"]) + " sec.", verbosity.low)
                warning(" @SOCKET:   Client " + str(c.peername) + " died or got unresponsive(B). Disconnecting.", verbosity.low)
                try:
                    c.shutdown(socket.SHUT_RDWR)
                except socket.error:
                    e = sys.exc_info()
                    warning(" @SOCKET:  could not shut down cleanly the socket. %s: %s in file '%s' on line %d" % (e[0].__name__, e[1], os.path.basename(e[2].tb_frame.f_code.co_filename), e[2].tb_lineno), verbosity.low)
                c.close()
                c.poll()
                c.status = Status.Disconnected

    def pool_distribute(self):
        """Deals with keeping the list of jobs up-to-date during a force
        calculation step.

        Deals with maintaining the jobs list. Gets data from drivers that have
        finished their calculation and removes that job from the list of running
        jobs, adds jobs to free clients and initialises the forcefields of new
        clients.
        """

        ttotal = tdispatch = tcheck = tjoin = 0
        ttotal -= time.time()

        # get clients that are still free
        with self.threadlock:
            freec = self.clients[:]
            for [r2, c] in self.jobs:
                freec.remove(c)

        # fills up list of pending requests if empty, or if clients are abundant
        if len(self.prlist) == 0 or len(freec)>len(self.prlist):
            self.prlist = [r for r in self.requests if r["status"] == "Queued"]

        if self.match_mode == "auto":
            match_seq = ["match", "none", "free", "any"]
        elif self.match_mode == "any":
            match_seq = ["any"]

        # first: dispatches jobs to free clients (if any!)
        # tries first to match previous replica<>driver association, then to get new clients, and only finally send the a new replica to old drivers
        sthreads = []
        ndispatch = 0
        tdispatch -= time.time()
        while len(freec) > 0 and len(self.prlist) > 0:
            for match_ids in match_seq:
                for fc in freec[:]:
                    if self.dispatch_free_client(fc, match_ids, sthreads):
                        freec.remove(fc)
                        ndispatch += 1
                    if len(self.prlist)==0:
                        break
            if len(freec)>0:
                self.prlist = [r for r in self.requests if r["status"] == "Queued"]
        tdispatch += time.time()

        # force a pool_update if there are requests pending
        # now check for client status
        if len(self.jobs)==0:
            for c in self.clients:
                if c.status == Status.Disconnected:  # client disconnected. force a pool_update
                    self.poll_iter = UPDATEFREQ
                    return
                if not c.status & (Status.Ready | Status.NeedsInit):
                    c.poll()

        # check for finished jobs
        cthreads = []
        nchecked = 0
        tcheck -= time.time()
        for [r, c] in self.jobs[:]:
            #self.check_job_finished(r,c)
            #ct = None
            ct = threading.Thread(target=self.check_job_finished, name="", args=(r,c))
            ct.daemon = True
            ct.start()
            cthreads.append((ct, r, c))
            nchecked += 1
        tcheck += time.time()

        tjoin -= time.time()
        self.join_threads(sthreads,cthreads)
        tjoin += time.time()

        ttotal += time.time()
        #print "POLL TOTAL: ", ttotal, "JOIN: ", tjoin, "DISPATCH:", ndispatch, tdispatch, "CHECK:", nchecked, tcheck

    def poll(self):
        """The main thread loop.

        Runs until either the program finishes or a kill call is sent. Updates
        the pool of clients every UPDATEFREQ loops and loops every latency
        seconds until _poll_true becomes false.
        """

        # makes sure to remove the last dead client as soon as possible -- and to get clients if we are dry
        if (self.poll_iter >= UPDATEFREQ or len(self.clients) == 0 or
             (len(self.clients) > 0 and not(self.clients[0].status & Status.Up))):
            self.poll_iter = 0
            self.pool_update()

        self.poll_iter += 1
        self.pool_distribute()

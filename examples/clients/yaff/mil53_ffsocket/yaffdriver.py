#!/usr/bin/env python

import socket
import numpy as np
from datetime import datetime

from ipi.interfaces.sockets import Message

# from yaff import ForceField
from yaff.log import log, timer

HDRLEN = 12
L_INT = 4
L_FLOAT = 8
L_CHAR = 1


class MySocket(object):
    def __init__(self, servername, mode="unix", port=31415, verbose=False):
        if mode == "unix":
            # Create a client socket
            self.s = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
            # Connect to i-PI server
            self.s.connect("/tmp/ipi_%s" % servername)
        elif mode == "inet":
            # Create a client socket
            self.s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            # Connect to i-PI server
            self.s.connect((servername, port))
        else:
            raise NotImplementedError
        self.s.setblocking(1)
        self.verbose = verbose

    def await_header(self):
        data = self.s.recv(HDRLEN * L_CHAR)
        if self.verbose:
            print("RECV", (data,))
        with log.section("DRIVER"):
            log("RECV %s %s" % (data, datetime.now()))
        return data

    def await_data(self, size):
        data = ""
        while len(data) < size:
            stub = self.s.recv(size - len(data))
            data += stub
            if self.verbose:
                print("RECV %i/%i" % (len(data), size))  # , hexlify(stub)
        return data

    def send_header(self, header):
        assert len(header) <= 12
        header = header.ljust(12, " ")
        if self.verbose:
            print("SEND", (header,))
        with log.section("DRIVER"):
            log("SEND %s %s" % (header, datetime.now()))
        self.s.send(header)

    def send_data(self, data):
        if self.verbose:
            print("SEND", len(data))  # , hexlify(data)
        self.s.send(data)


class YAFFDriver(object):
    """
    Use Yaff as a driver that calculates forces for i-PI
    """

    log_name = "DRIVER"

    def __init__(self, socket, ff):
        self.ff = ff
        self.s = socket
        self.isinit = False
        self.hasdata = False
        self.natom = -1
        self.rvecs = np.zeros((3, 3))
        self.gvecs = np.zeros((3, 3))

    def receive_pos(self):
        self.rvecs[:] = (
            np.fromstring(self.s.await_data(L_FLOAT * 9), np.float64).reshape((3, 3)).T
        )
        self.gvecs[:] = np.fromstring(
            self.s.await_data(L_FLOAT * 9), np.float64
        ).reshape((3, 3))
        natom = np.fromstring(self.s.await_data(L_INT), np.int32)[0]
        if self.natom < 0:
            self.natom = natom
            self.pos = np.zeros((self.natom, 3))
            self.gpos = np.zeros((self.natom, 3))
            self.vtens = np.zeros((3, 3))
        self.pos[:] = np.fromstring(
            self.s.await_data(L_FLOAT * 3 * self.natom), np.float64
        ).reshape((natom, 3))

    def compute(self):
        self.ff.update_rvecs(self.rvecs)
        self.ff.update_pos(self.pos)
        self.ff.nlist.update()
        self.gpos[:] = 0.0
        self.vtens[:] = 0.0
        self.e = self.ff.compute(self.gpos, self.vtens)
        # Sign conventions for i-PI
        self.gpos[:] *= -1
        self.vtens[:] *= -1
        self.hasdata = True

    def send_forces(self):
        self.s.send_header(Message("forceready"))
        self.s.send_data(np.array([self.e]))
        self.s.send_data(np.array([self.natom]))
        self.s.send_data(self.gpos.ravel())
        self.s.send_data(self.vtens.ravel())
        self.s.send_data(np.array([np.int32(9)]))
        self.s.send_data("YAFF Done")
        self.hasdata = False

    def run(self):
        # Run indefinitely (until i-PI sends exit message)
        while True:
            with timer.section("WAITHEADER"):
                # log("WAITHEADER %s" % datetime.now())
                header = self.s.await_header()
                # log("GOTHEADER %s" % datetime.now())
            if header == Message("status"):
                if not self.isinit:
                    self.s.send_header(Message("needinit"))
                elif self.hasdata:
                    self.s.send_header(Message("havedata"))
                else:
                    self.s.send_header(Message("ready"))
            elif header == Message("init"):
                # ibead = np.fromstring(self.s.await_data(L_INT), np.int32)[0]
                # len_init = np.fromstring(self.s.await_data(L_INT), np.int32)[0]
                # init = self.s.await_data(len_init * L_CHAR)
                if log.do_high:
                    with log.section(self.log_name):
                        # log( "YAFF driver initialized for pid %s (bead %d - %s)" % (os.getpid(), ibead, init) )
                        log("INIT %s" % datetime.now())
                self.isinit = True
            elif header == Message("posdata"):
                with timer.section("RECVPOS"), log.section(self.log_name):
                    self.receive_pos()
                    log("RECV %-12s %s" % ("POS", datetime.now()))
                with timer.section("COMPUTE"):
                    self.compute()
            elif header == Message("getforce"):
                with timer.section("SENDFORCE"), log.section(self.log_name):
                    self.send_forces()
                    log("SEND %-12s %s" % ("FORCE", datetime.now()))
            elif header == Message("exit"):
                if log.do_high:
                    with log.section(self.log_name):
                        # log('i-PI finished, stopping Yaff driver')
                        log("EXIT %s" % datetime.now())
                break
            else:
                raise NotImplementedError("Received unknown message %s" % header)

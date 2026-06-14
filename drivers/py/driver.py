#!/usr/bin/env python3
import socket
import argparse
import numpy as np
from multiprocessing import shared_memory, resource_tracker
from ipi.pes import Dummy_driver, load_pes, __drivers__
from ipi.utils.io.inputs import read_args_kwargs

description = """
Minimal example of a Python driver connecting to i-PI and exchanging energy, forces, etc.
"""


def recv_data(sock, data):
    """Fetches binary data from i-PI socket."""
    blen = data.itemsize * data.size
    buf = np.zeros(blen, np.byte)

    bpos = 0
    while bpos < blen:
        timeout = False
        try:
            bpart = 1
            bpart = sock.recv_into(buf[bpos:], blen - bpos)
        except socket.timeout:
            print(" @SOCKET:   Timeout in status recvall, trying again!")
            timeout = True
            pass
        if not timeout and bpart == 0:
            raise RuntimeError("Socket disconnected!")
        bpos += bpart
    if np.isscalar(data):
        return np.frombuffer(buf[0:blen], data.dtype)[0]
    else:
        return np.frombuffer(buf[0:blen], data.dtype).reshape(data.shape)


def send_data(sock, data):
    """Sends binary data to i-PI socket."""

    if np.isscalar(data):
        data = np.array([data], data.dtype)
    buf = data.tobytes()
    sock.send(buf)


HDRLEN = 12  # number of characters of the default message strings


def Message(mystr):
    """Returns a header of standard length HDRLEN."""

    # convert to bytestream since we'll be sending this over a socket
    return str.ljust(str.upper(mystr), HDRLEN).encode()


def run_driver(
    unix=False,
    address="",
    port=12345,
    driver=Dummy_driver(),
    f_verbose=False,
    sockets_prefix="/tmp/ipi_",
    shm=False,
):
    """Minimal socket client for i-PI."""

    # Opens a socket to i-PI. shm uses a unix socket for the control handshake;
    # only the bulk payload travels through shared memory (same node only).
    if unix or shm:
        sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        sock.connect(sockets_prefix + address)
    else:
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        # this reduces latency for the small messages passed by the i-PI protocol
        sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
        sock.connect((address, port))

    f_init = False
    f_data = False

    # batched evaluation: batch_n>1 is announced by i-PI in the INIT string
    batch_n = 1
    results_batch = None

    # shared-memory transport: segments are attached lazily on the first POSDATA
    shm_kinds = ("pos", "h", "ih", "pot", "force", "vir")
    shm_handles = {}
    shm_views = {}

    # initializes structure arrays
    cell = np.zeros((3, 3), float)
    icell = np.zeros((3, 3), float)
    pos = np.zeros(0, float)

    # initializes return arrays
    pot = 0.0
    force = np.zeros(0, float)
    vir = np.zeros((3, 3), float)
    while True:  # ah the infinite loop!
        header = sock.recv(HDRLEN)
        if f_verbose:
            print("Received ", header)
        if header == Message("STATUS"):
            # responds to a status request
            if not f_init:
                sock.sendall(Message("NEEDINIT"))
            elif f_data:
                sock.sendall(Message("HAVEDATA"))
            else:
                sock.sendall(Message("READY"))
        elif header == Message("INIT"):
            # initialization
            rid = recv_data(sock, np.int32())
            initlen = recv_data(sock, np.int32())
            initstr = recv_data(sock, np.empty(initlen, dtype="S1"))
            # i-PI announces request batching via a reserved 'batch_size' token
            init_text = initstr.tobytes().decode("utf-8", errors="replace")
            for token in init_text.split(","):
                k, sep, v = token.partition(":")
                if sep and k.strip() == "batch_size":
                    try:
                        batch_n = int(v.strip())
                    except ValueError:
                        pass
            if f_verbose:
                print(rid, initstr, "batch_size:", batch_n)
            f_init = True  # we are initialized now
        elif header == Message("POSDATA") and shm:
            # shared-memory transport: only natoms (and, on the first call, the
            # segment names) arrive on the socket; the batch_n structures are
            # read straight from shared memory. Handles single (batch_n=1) and
            # batched runs with one code path.
            nat = recv_data(sock, np.int32())
            if not shm_views:
                for kind in shm_kinds:
                    nmlen = recv_data(sock, np.int32())
                    name = (
                        recv_data(sock, np.empty(nmlen, dtype="S1"))
                        .tobytes()
                        .decode("utf-8")
                    )
                    handle = shared_memory.SharedMemory(name=name)
                    # only attaching: let i-PI own the segment lifetime
                    resource_tracker.unregister(handle._name, "shared_memory")
                    shm_handles[kind] = handle
                shm_views["pos"] = np.ndarray(
                    (batch_n, nat, 3), np.float64, buffer=shm_handles["pos"].buf
                )
                shm_views["h"] = np.ndarray(
                    (batch_n, 3, 3), np.float64, buffer=shm_handles["h"].buf
                )
                shm_views["ih"] = np.ndarray(
                    (batch_n, 3, 3), np.float64, buffer=shm_handles["ih"].buf
                )
                shm_views["pot"] = np.ndarray(
                    (batch_n,), np.float64, buffer=shm_handles["pot"].buf
                )
                shm_views["force"] = np.ndarray(
                    (batch_n, nat, 3), np.float64, buffer=shm_handles["force"].buf
                )
                shm_views["vir"] = np.ndarray(
                    (batch_n, 3, 3), np.float64, buffer=shm_handles["vir"].buf
                )
            cell_list = [shm_views["h"][i] for i in range(batch_n)]
            pos_list = [shm_views["pos"][i] for i in range(batch_n)]

            ##### THIS IS THE TIME TO DO SOMETHING WITH THE POSITIONS!
            results_batch = driver(cell_list, pos_list)
            f_data = True
        elif header == Message("POSDATA") and batch_n > 1:
            # batched structural information: a single atom count, then batch_n
            # cell blocks (h, ih) followed by a (batch_n, 3*nat) position array.
            nat = recv_data(sock, np.int32())
            cellbuf = recv_data(sock, np.zeros(batch_n * 18, np.float64))
            posbuf = recv_data(sock, np.zeros((batch_n, nat * 3), np.float64))
            cell_list = []
            pos_list = []
            for i in range(batch_n):
                cell_list.append(cellbuf[i * 18 : i * 18 + 9].reshape(3, 3))
                pos_list.append(posbuf[i].reshape(nat, 3))

            ##### THIS IS THE TIME TO DO SOMETHING WITH THE POSITIONS!
            results_batch = driver(cell_list, pos_list)
            f_data = True
        elif header == Message("POSDATA"):
            # receives structural information
            cell = recv_data(sock, cell)
            icell = recv_data(
                sock, icell
            )  # inverse of the cell. mostly useless legacy stuff
            nat = recv_data(sock, np.int32())
            if len(pos) == 0:
                # shapes up the position array
                pos.resize((nat, 3), refcheck=False)
                force.resize((nat, 3), refcheck=False)
            else:
                if len(pos) != nat:
                    raise RuntimeError("Atom number changed during i-PI run")
            pos = recv_data(sock, pos)

            ##### THIS IS THE TIME TO DO SOMETHING WITH THE POSITIONS!
            pot, force, vir, extras = driver(cell, pos)
            f_data = True
        elif header == Message("GETFORCE") and shm:
            # write the batch_n results into shared memory BEFORE the FORCEREADY
            # ack, so i-PI never reads a stale buffer; only the per-structure
            # extra strings travel on the socket.
            for i in range(batch_n):
                res = results_batch[i]
                shm_views["pot"][i] = res[0]
                shm_views["force"][i] = np.asarray(res[1]).reshape(nat, 3)
                shm_views["vir"][i] = np.asarray(res[2]).reshape(3, 3)
            sock.sendall(Message("FORCEREADY"))
            for res in results_batch:
                extra = res[3]
                send_data(sock, np.int32(len(extra)))
                sock.sendall(extra.encode("utf-8"))
            f_data = False
        elif header == Message("GETFORCE") and batch_n > 1:
            sock.sendall(Message("FORCEREADY"))
            # batched reply: potentials, atom count, forces, virials, then one
            # length-prefixed extra string per structure.
            pots = np.array([res[0] for res in results_batch], np.float64)
            send_data(sock, pots)
            send_data(sock, np.int32(nat))
            forces = np.array(
                [np.asarray(res[1]).reshape(nat * 3) for res in results_batch],
                np.float64,
            )
            send_data(sock, forces)
            virs = np.array(
                [np.asarray(res[2]).reshape(3, 3) for res in results_batch],
                np.float64,
            )
            send_data(sock, virs)
            for res in results_batch:
                extra = res[3]
                send_data(sock, np.int32(len(extra)))
                sock.sendall(extra.encode("utf-8"))
            f_data = False
        elif header == Message("GETFORCE"):
            sock.sendall(Message("FORCEREADY"))

            # sanity check in the returned values (catches bugs and inconsistencies in the implementation)
            if not isinstance(force, np.ndarray) and force.dtype == np.float64:
                raise ValueError(
                    "driver returned forces with the wrong type: we need a "
                    "numpy.ndarray containing 64-bit floating points values"
                )

            if not isinstance(vir, np.ndarray) and vir.dtype == np.float64:
                raise ValueError(
                    "driver returned virial with the wrong type: we need a "
                    "numpy.ndarray containing 64-bit floating points values"
                )

            if len(force.flatten()) != len(pos.flatten()):
                raise ValueError(
                    "driver returned forces with the wrong size: number of "
                    "atoms and dimensions must match positions"
                )

            if len(vir.flatten()) != 9:
                raise ValueError(
                    "driver returned a virial tensor which does not have 9 components"
                )

            send_data(sock, np.float64(pot))
            send_data(sock, np.int32(nat))
            send_data(sock, force)
            send_data(sock, vir)
            send_data(sock, np.int32(len(extras)))
            sock.sendall(extras.encode("utf-8"))

            f_data = False
        elif header == Message("EXIT"):
            print("Received exit message from i-PI. Bye bye!")
            if shm:
                for handle in shm_handles.values():
                    try:
                        handle.close()
                    except Exception:
                        pass
            return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "-u",
        "--unix",
        action="store_true",
        default=False,
        help="Use a UNIX domain socket.",
    )
    parser.add_argument(
        "-a",
        "--address",
        type=str,
        default="localhost",
        help="Host name (for INET sockets) or name of the UNIX domain socket to connect to.",
    )
    parser.add_argument(
        "-S",
        "--sockets_prefix",
        type=str,
        default="/tmp/ipi_",
        help="Prefix used for the unix domain sockets. Ignored when using TCP/IP sockets.",
    )
    parser.add_argument(
        "-p",
        "--port",
        type=int,
        default=12345,
        help="TCP/IP port number. Ignored when using UNIX domain sockets.",
    )
    parser.add_argument(
        "-m",
        "--mode",
        type=str,
        default="dummy",
        choices=list(__drivers__.keys()) + ["custom"],
        help="""Type of potential to be used to compute the potential and its derivatives.
        """,
    )
    parser.add_argument(
        "-P",
        "--pes_path",
        type=str,
        default=None,
        help="""File path for 'custom' PES (it should end with .py).
        """,
    )
    parser.add_argument(
        "-o",
        "--param",
        type=str,
        default="",
        help="""Parameters required to run the driver. Comma-separated list of values
        """,
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="Verbose output.",
    )
    parser.add_argument(
        "--shm",
        action="store_true",
        default=False,
        help="Exchange the bulk position/force payload through shared memory "
        "(implies a UNIX domain control socket; same-node only). The batch size, "
        "if any, is taken from the INIT string.",
    )

    args = parser.parse_args()

    driver_args, driver_kwargs = read_args_kwargs(args.param)

    # import the driver class
    cls = load_pes(args.mode, args.pes_path)

    d_f = cls(*driver_args, **driver_kwargs)

    run_driver(
        unix=args.unix,
        address=args.address,
        port=args.port,
        driver=d_f,
        f_verbose=args.verbose,
        sockets_prefix=args.sockets_prefix,
        shm=args.shm,
    )

#!/usr/bin/env python3
import socket
import argparse
import numpy as np
from ipi.pes import Dummy_driver, load_pes, __drivers__
from ipi.utils.io.inputs import read_args_kwargs
from multiprocessing import resource_tracker, shared_memory


description = """
Minimal example of a Python driver connecting to i-PI and exchanging energy, forces, etc.
"""

HDRLEN = 12  # number of characters of the default message strings


def Message(mystr):
    """Returns a header of standard length HDRLEN."""

    # convert to bytestream since we'll be sending this over a socket
    return str.ljust(str.upper(mystr), HDRLEN).encode()



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


def open_driver_socket(unix, address, port, sockets_prefix):
    """Opens the control socket to i-PI."""

    if unix:
        sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        sock.connect(sockets_prefix + address)
    else:
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
        sock.connect((address, port))
    return sock


def validate_force_payload(force, vir, pos):
    """Checks returned force/virial arrays before sending them back."""

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

    if force.size != pos.size:
        raise ValueError(
            "driver returned forces with the wrong size: number of "
            "atoms and dimensions must match positions"
        )

    if vir.size != 9:
        raise ValueError(
            "driver returned a virial tensor which does not have 9 components"
        )


def recv_fixed(sock, length):
    """Receives an exact number of bytes."""

    return recv_data(sock, np.empty(length, dtype=np.uint8)).tobytes()


def attach_shared_memory(name):
    """Attach to server-owned shared memory without tracking it for unlink."""

    shm = shared_memory.SharedMemory(name=name)
    resource_tracker.unregister(shm._name, "shared_memory")
    return shm


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

    sock = open_driver_socket(unix, address, port, sockets_prefix)

    f_init = False
    f_data = False

    # initializes structure arrays
    cell = np.zeros((3, 3), float)
    icell = np.zeros((3, 3), float)
    pos = np.zeros(0, float)

    # initializes return arrays
    pot = 0.0
    force = np.zeros(0, float)
    vir = np.zeros((3, 3), float)
    extras = ""

    cell_shm = None
    icell_shm = None
    pot_shm = None
    pos_shm = None
    f_shm = None
    vir_shm = None
    pot_snp = None
    pos_snp = None
    f_snp = None
    vir_snp = None

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
            initstr = recv_data(sock, np.chararray(initlen))

            if f_verbose:
                print(rid, "Initing...")

            if shm:
                nat = recv_data(sock, np.int32())
                nbeads = recv_data(sock, np.int32())
                if nbeads != 1:
                    raise RuntimeError(
                        "run_driver(shm=True) only supports single-bead shm mode; "
                        "use the MPI shm driver for batched mode"
                    )

                pos_bufname = recv_fixed(sock, HDRLEN).decode("utf-8").strip()
                h_bufname = recv_fixed(sock, HDRLEN).decode("utf-8").strip()
                ih_bufname = recv_fixed(sock, HDRLEN).decode("utf-8").strip()
                pot_bufname = recv_fixed(sock, HDRLEN).decode("utf-8").strip()
                f_bufname = recv_fixed(sock, HDRLEN).decode("utf-8").strip()
                vir_bufname = recv_fixed(sock, HDRLEN).decode("utf-8").strip()

                if f_verbose:
                    print(rid, initstr)
                    print(f"natoms, nbeads: {nat}, {nbeads}")
                    print("Driver initing SHM with buffer names:")
                    print(pos_bufname, h_bufname, ih_bufname, pot_bufname, f_bufname, vir_bufname)

                cell_shm = attach_shared_memory(h_bufname)
                icell_shm = attach_shared_memory(ih_bufname)
                pot_shm = attach_shared_memory(pot_bufname)
                pos_shm = attach_shared_memory(pos_bufname)
                f_shm = attach_shared_memory(f_bufname)
                vir_shm = attach_shared_memory(vir_bufname)

                cell_snp = np.ndarray((3, 3), dtype=np.float64, buffer=cell_shm.buf)
                icell_snp = np.ndarray((3, 3), dtype=np.float64, buffer=icell_shm.buf)
                pot_snp = np.ndarray((1), dtype=np.float64, buffer=pot_shm.buf)
                pos_snp = np.ndarray((nat * 3), dtype=np.float64, buffer=pos_shm.buf)
                f_snp = np.ndarray((nat * 3), dtype=np.float64, buffer=f_shm.buf)
                vir_snp = np.ndarray((3, 3), dtype=np.float64, buffer=vir_shm.buf)

            f_init = True  # we are initialized now
        elif header == Message("POSDATA"):
            if shm:
                pos = pos_snp
                cell = cell_snp
                icell = icell_snp
            else:
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
        elif header == Message("GETFORCE"):
            validate_force_payload(force, vir, pos)

            if shm:
                pot_snp[:] = pot
                f_snp[:] = force
                vir_snp[:] = vir
                sock.sendall(Message("FORCEREADY"))
            else:
                sock.sendall(Message("FORCEREADY"))
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
                cell_shm.close()
                icell_shm.close()
                pot_shm.close()
                pos_shm.close()
                f_shm.close()
                vir_shm.close()
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
        "-s",
        "--shm",
        action="store_true",
        default=False,
        help="Use shared memory communication",
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

    args = parser.parse_args()

    driver_args, driver_kwargs = read_args_kwargs(args.param)

    # import the driver class
    cls = load_pes(args.mode, args.pes_path)

    d_f = cls(*driver_args, **driver_kwargs)

    if args.shm:
        run_driver(
            unix=True,
            address=args.address,
            port=args.port,
            driver=d_f,
            f_verbose=args.verbose,
            sockets_prefix=args.sockets_prefix,
            shm=True,
        )
    else:
        run_driver(
            unix=args.unix,
            address=args.address,
            port=args.port,
            driver=d_f,
            f_verbose=args.verbose,
            sockets_prefix=args.sockets_prefix,
        )

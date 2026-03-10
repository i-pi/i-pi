#!/usr/bin/env python3
import socket
import argparse
import numpy as np
from ipi.pes import Dummy_driver, load_pes, __drivers__
from ipi.utils.io.inputs import read_args_kwargs
from multiprocessing import shared_memory


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


def run_driver(
    unix=False,
    address="",
    port=12345,
    driver=Dummy_driver(),
    f_verbose=False,
    sockets_prefix="/tmp/ipi_",
):
    """Minimal socket client for i-PI."""

    # Opens a socket to i-PI
    if unix:
        sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        sock.connect(sockets_prefix + address)
    else:
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        # this reduces latency for the small messages passed by the i-PI protocol
        sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
        sock.connect((address, port))

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
            #initstr = recv_data(sock, np.empty(initlen, dtype=np.uint8))

            if f_verbose:
                print(rid, "Initing...")
            f_init = True  # we are initialized now
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

            if force.size != pos.size:
                raise ValueError(
                    "driver returned forces with the wrong size: number of "
                    "atoms and dimensions must match positions"
                )

            if vir.size != 9:
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
            return

def run_shmdriver(
    unix=True,
    address="",
    port=12345,
    driver=Dummy_driver(),
    f_verbose=False,
    sockets_prefix="/tmp/ipi_",
):
    """Minimal socket client for i-PI."""

    # Opens a socket to i-PI
    if unix:
        sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        sock.connect(sockets_prefix + address)
    else:
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        # this reduces latency for the small messages passed by the i-PI protocol
        sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
        sock.connect((address, port))

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
            print("start initing")
            rid = recv_data(sock, np.int32())
            initlen = recv_data(sock, np.int32())
            initstr = recv_data(sock, np.chararray(initlen))
            print("initstr")
                  
            # receiving nat and nbeads at init, necessary to be able to allocate the correct numpy shapes in SHM
            nat = recv_data(sock, np.int32())
            nbeads = recv_data(sock, np.int32())
          
            if f_verbose:
                print(rid, initstr)
                print(f"natoms, nbeads: {nat}, {nbeads}")
            
            # reading buffer names for shaerd memory access
            pos_bufname = sock.recv(HDRLEN).decode("utf-8").strip()
            h_bufname = sock.recv(HDRLEN).decode("utf-8").strip()
            ih_bufname = sock.recv(HDRLEN).decode("utf-8").strip()
            pot_bufname = sock.recv(HDRLEN).decode("utf-8").strip()
            f_bufname = sock.recv(HDRLEN).decode("utf-8").strip()
            vir_bufname = sock.recv(HDRLEN).decode("utf-8").strip()

            print("Driver initing SHM with buffer names:")
            print(pos_bufname, h_bufname, ih_bufname, pot_bufname, f_bufname, vir_bufname)

            #shm objects
            cell_shm = shared_memory.SharedMemory(name=h_bufname)
            icell_shm = shared_memory.SharedMemory(name=ih_bufname)
            pot_shm = shared_memory.SharedMemory(name=pot_bufname)
            pos_shm = shared_memory.SharedMemory(name=pos_bufname)
            f_shm = shared_memory.SharedMemory(name=f_bufname)
            vir_shm = shared_memory.SharedMemory(name=vir_bufname)


            # allocating numpy arrays in shared memory on the same buffer as in the server
            cell_snp = np.ndarray((3,3), dtype=np.float64, buffer=cell_shm.buf)
            icell_snp = np.ndarray((3,3), dtype=np.float64, buffer=icell_shm.buf)
            pot_snp = np.ndarray((1), dtype=np.float64, buffer=pot_shm.buf)
            pos_snp = np.ndarray((nat*3), dtype=np.float64, buffer=pos_shm.buf)
            f_snp = np.ndarray((nat*3), dtype=np.float64, buffer=f_shm.buf)
            vir_snp = np.ndarray((3,3), dtype=np.float64, buffer=vir_shm.buf)

            
            
            f_init = True  # we are initialized now
        
        elif header == Message("POSDATA"):
            # receives structural information
            pos = pos_snp[:]
            cell = cell_snp[:]
            icell = icell_snp[:]
            
            ''' 
            if len(pos) == 0:
                # shapes up the position array
                pos.resize((nat, 3), refcheck=False)
                force.resize((nat, 3), refcheck=False)
            else:
                if len(pos) != nat:
                    raise RuntimeError("Atom number changed during i-PI run")
            '''

            ##### THIS IS THE TIME TO DO SOMETHING WITH THE POSITIONS!
            pot, force, vir, extras = driver(cell, pos)
            f_data = True
        elif header == Message("GETFORCE"):

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

            if force.size != pos.size:
                raise ValueError(
                    "driver returned forces with the wrong size: number of "
                    "atoms and dimensions must match positions"
                )

            if len(vir.flatten()) != 9:
                raise ValueError(
                    "driver returned a virial tensor which does not have 9 components"
                )

            pot_snp = pot
            f_snp = force[:]
            vir_snp = vir[:]

            sock.sendall(Message("FORCEREADY"))
            f_data = False
        elif header == Message("EXIT"):
            print("Received exit message from i-PI. Bye bye!")
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
        run_shmdriver(
            unix=True,
            address=args.address,
            port=args.port,
            driver=d_f,
            f_verbose=args.verbose,
            sockets_prefix=args.sockets_prefix,
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

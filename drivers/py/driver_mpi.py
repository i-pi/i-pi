#!/usr/bin/env python3
'''
import builtins, traceback

_real_import = builtins.__import__

def tracing_import(name, *args, **kwargs):
    if name.startswith("jax"):
        print(f"\n*** JAX IMPORT TRIGGERED: {name} ***")
        traceback.print_stack()
    return _real_import(name, *args, **kwargs)

builtins.__import__ = tracing_import
'''


import socket
import argparse
import numpy as np
import time


from ipi.pes import load_pes, __drivers__

from ipi.utils.io.inputs import read_args_kwargs
from mpi4py import MPI
from ipi.utils.timing_manager import timers
from multiprocessing import resource_tracker, shared_memory

description = """
Minimal MPI Python driver for i-PI batch mode.

Rank 0 handles the socket control protocol. In plain unix/inet mode it
receives one batched position payload and distributes work across MPI ranks.
In SHM mode it can either use Scatter/Gather or, with ``--shm-local-slices``,
let one rank per bead read and write its own shared-memory slice directly.
"""
def recv_data(sock, data):
    """Fetches binary data from i-PI socket."""
    blen = data.itemsize * data.size
    buf = np.zeros(blen, np.byte)

    bpos = 0
    while bpos < blen:
        timeout = False
        try:
            bpart = sock.recv_into(buf[bpos:], blen - bpos)
        except socket.timeout:
            print(" @SOCKET:   Timeout in status recvall, trying again!")
            timeout = True
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
    sock.sendall(buf)

def send_datas(sock, *args, extras=None):
    """Send multiple scalars/arrays and an optional extras string in one go."""
    buffers = []
    for data in args:
        if np.isscalar(data):
            data = np.array([data], dtype=type(data))
        buffers.append(data.tobytes())

    sock.send(b"".join(buffers))

    if extras is not None:
        sock.sendall(extras.encode("utf-8"))

HDRLEN = 12  # number of characters of the default message strings

def Message(mystr):
    """Returns a header of standard length HDRLEN."""
    return str.ljust(str.upper(mystr), HDRLEN).encode()


def attach_shared_memory(name):
    """Attach to server-owned shared memory without tracking it for unlink."""

    shm = shared_memory.SharedMemory(name=name)
    resource_tracker.unregister(shm._name, "shared_memory")
    return shm

def run_driver(
    unix=False,
    address="",
    port=12345,
    driver=None,
    f_verbose=False,
    f_debug=False,
    sockets_prefix="/tmp/ipi_",
    shm=True,
    shm_local_slices=False,
):
    """Batch socket client for i-PI with optional SHM data exchange.

    Rank 0 handles the socket protocol. In plain unix/inet mode it receives
    one batched request and distributes work across MPI ranks. In SHM mode it
    maps the server-owned buffers and either uses Scatter/Gather or direct
    per-rank bead slices before signaling readiness to i-PI.
    """

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    t_start = time.perf_counter()
    
    def dbg_print(msg):
        if f_debug and rank == 0:
            elapsed = time.perf_counter() - t_start
            print(f"[{elapsed:8.4f} s] Rank {rank}: {msg}", flush=True)

    # Only rank 0 opens and communicates on the socket
    if rank == 0:
        dbg_print("Connecting to socket...")
        if unix:
            sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
            sock.connect(sockets_prefix + address)
        else:
            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
            sock.connect((address, port))
        dbg_print("Socket connected.")
    else:
        sock = None  # Other ranks don't have socket

    f_init = False
    f_data = False
    request_count = 0

    # Initialize arrays (all ranks keep them)
    cell = np.zeros((3, 3), float)
    icell = np.zeros((3, 3), float)
    pos = np.zeros(0, float)

    pot = 0.0
    force = np.zeros(0, float)
    vir = np.zeros((3, 3), float)

    cell_shm = None
    icell_shm = None
    pot_shm = None
    pos_shm = None
    f_shm = None
    vir_shm = None
    cell_snp = None
    icell_snp = None
    pot_snp = None
    pos_snp = None
    f_snp = None
    vir_snp = None

    while True:
        timers.set_step(request_count)
        if rank == 0:
            timers.start("Recv Header")
            header = sock.recv(HDRLEN)
            if f_verbose:
                print(f"Rank {rank} received ", header)
            timers.stop("Recv Header")
        else:
            header = None

        header = comm.bcast(header, root=0)

        if header == Message("STATUS"):
            if rank == 0:
                if not f_init:
                    sock.sendall(Message("NEEDINIT"))
                elif f_data:
                    sock.sendall(Message("HAVEDATA"))
                else:
                    sock.sendall(Message("READY"))
                    if f_verbose:
                        print(f"Driver ready")

        elif header == Message("INIT"):
            if rank == 0:
                rid = recv_data(sock, np.int32())
                initlen = recv_data(sock, np.int32())
                initstr = recv_data(sock, np.zeros(initlen, np.dtype("S1")))
                init_payload = None
                if shm:
                    # receiving nat and nbeads at init, necessary to be able to allocate the correct numpy shapes in SHM
                    nat = recv_data(sock, np.int32())
                    nbeads = recv_data(sock, np.int32())

                    # reading buffer names for shared memory access
                    pos_bufname = sock.recv(HDRLEN).decode("utf-8").strip()
                    h_bufname = sock.recv(HDRLEN).decode("utf-8").strip()
                    ih_bufname = sock.recv(HDRLEN).decode("utf-8").strip()
                    pot_bufname = sock.recv(HDRLEN).decode("utf-8").strip()
                    f_bufname = sock.recv(HDRLEN).decode("utf-8").strip()
                    vir_bufname = sock.recv(HDRLEN).decode("utf-8").strip()

                    print("Driver initing SHM with buffer names:")
                    print(pos_bufname, h_bufname, ih_bufname, pot_bufname, f_bufname, vir_bufname)

                    init_payload = (
                        int(nat),
                        int(nbeads),
                        pos_bufname,
                        h_bufname,
                        ih_bufname,
                        pot_bufname,
                        f_bufname,
                        vir_bufname,
                    )
                
                if f_verbose:
                    print(f"Rank {rank} initing with rid, initstr: {rid}, {initstr}")
                f_init = True
            else:
                init_payload = None

            if shm:
                init_payload = comm.bcast(init_payload, root=0)
                (
                    nat,
                    nbeads,
                    pos_bufname,
                    h_bufname,
                    ih_bufname,
                    pot_bufname,
                    f_bufname,
                    vir_bufname,
                ) = init_payload

                if shm_local_slices and size != int(nbeads):
                    raise RuntimeError(
                        "SHM local-slices mode requires one MPI rank per bead "
                        f"(got size={size}, nbeads={nbeads})."
                    )

                if rank == 0 or shm_local_slices:
                    cell_shm = attach_shared_memory(h_bufname)
                    icell_shm = attach_shared_memory(ih_bufname)
                    pot_shm = attach_shared_memory(pot_bufname)
                    pos_shm = attach_shared_memory(pos_bufname)
                    f_shm = attach_shared_memory(f_bufname)
                    vir_shm = attach_shared_memory(vir_bufname)

                    # allocating numpy arrays in shared memory on the same buffer as in the server
                    cell_snp = np.ndarray((3,3), dtype=np.float64, buffer=cell_shm.buf)
                    icell_snp = np.ndarray((3,3), dtype=np.float64, buffer=icell_shm.buf)
                    pot_snp = np.ndarray((nbeads), dtype=np.float64, buffer=pot_shm.buf)
                    pos_snp = np.ndarray((nbeads, nat*3), dtype=np.float64, buffer=pos_shm.buf)
                    f_snp = np.ndarray((nbeads, nat*3), dtype=np.float64, buffer=f_shm.buf)
                    vir_snp = np.ndarray((nbeads, 3,3), dtype=np.float64, buffer=vir_shm.buf)
            f_init = comm.bcast(f_init, root=0)

        elif header == Message("POSDATA"):
            request_count += 1
            timers.set_step(request_count)
            if rank == 0:
                timers.start("Read Posdata")
                if shm:
                    pos = pos_snp
                    cell = cell_snp
                    icell = icell_snp
                else:
                    cell = recv_data(sock, cell)
                    icell = recv_data(sock, icell)
                    nat_token = int(recv_data(sock, np.int32()))
                    if nat_token >= 0:
                        raise RuntimeError(
                            "MPI driver plain-socket mode expects batched positions"
                        )
                    nat = -nat_token
                    nbeads = int(recv_data(sock, np.int32()))
                    pos = recv_data(
                        sock, np.zeros((nbeads, nat * 3), dtype=np.float64)
                    )
                
                timers.stop("Read Posdata")
            else:
                if not (shm and shm_local_slices):
                    nbeads = None
                    nat = None

            if shm and shm_local_slices:
                pos_local = pos_snp[rank].reshape((nat, 3))
                cell = cell_snp
                icell = icell_snp

                timers.start("Run Driver")
                results = driver(cell, pos_local)
                timers.stop("Run Driver")

                pot_snp[rank] = results[0]
                f_snp[rank, :] = np.asarray(results[1], dtype=np.float64).reshape(nat * 3)
                vir_snp[rank, :, :] = np.asarray(results[2], dtype=np.float64).reshape(3, 3)

                comm.Barrier()

                if rank == 0:
                    pot = pot_snp
                    force = f_snp
                    vir = vir_snp
            else:
                timers.start("MPI Bcast Metadata")

                # fist, broadcast nat, nbeads, and cell data to other ranks
                if rank == 0:
                        buf = np.concatenate([
                            np.atleast_1d(np.int32(nbeads)),
                            np.atleast_1d(np.int32(nat)),
                            cell.ravel(),
                            icell.ravel(),
                        ])
                else:
                    buf = np.empty(20, dtype=np.float64)
                
                comm.Bcast([buf, MPI.DOUBLE], root=0)

                # unpack after broadcast
                nbeads = np.int32(buf[0])
                nat = np.int32(buf[1])
                cell = buf[2:11].reshape(3,3)
                icell = buf[11:].reshape(3,3)
                
                timers.stop("MPI Bcast Metadata")
                
                if rank == 0:
                    # pos is shape (nbeads, nat, 3)
                    sendbuf = pos
                else:
                    sendbuf = None

                # Allocate local slice on each worker
                recvbuf = np.empty((nat, 3), dtype=np.float64)

                timers.start("MPI Scatter Pos")
                comm.Scatter(sendbuf, recvbuf, root=0)
                pos_local = recvbuf
                timers.stop("MPI Scatter Pos")

                timers.start("Run Driver")
                results = driver(cell, pos_local)
                
                timers.stop("Run Driver")

                local_pot = np.asarray([results[0]], dtype=np.float64)
                local_force = np.asarray(results[1], dtype=np.float64).reshape(nat * 3)
                local_vir = np.asarray(results[2], dtype=np.float64).reshape(9)

                if rank == 0:
                    if shm:
                        # Gather directly into the SHM-backed result arrays so
                        # GETFORCE only needs to notify i-PI that results are ready.
                        pot = pot_snp
                        force = f_snp.reshape((nbeads, nat * 3))
                        vir = vir_snp.reshape((nbeads, 9))
                    else:
                        pot = np.empty(size, dtype=np.float64)
                        force = np.empty((size, nat * 3), dtype=np.float64)
                        vir = np.empty((size, 9), dtype=np.float64)
                else:
                    pot = None
                    force = None
                    vir = None

                timers.start("MPI Gather Results")
                comm.Gather([local_pot, MPI.DOUBLE], [pot, MPI.DOUBLE], root=0)
                comm.Gather([local_force, MPI.DOUBLE], [force, MPI.DOUBLE], root=0)
                comm.Gather([local_vir, MPI.DOUBLE], [vir, MPI.DOUBLE], root=0)
                extras = comm.gather(results[3], root=0)
                timers.stop("MPI Gather Results")

                if rank == 0:
                    pot = np.asarray(pot, dtype=np.float64)
                    force = np.asarray(force, dtype=np.float64).reshape((nbeads, nat * 3))
                    vir = np.asarray(vir, dtype=np.float64).reshape((nbeads, 3, 3))

            f_data = True

        elif header == Message("GETFORCE"):
            if rank == 0:
                # we update SHM first, than send message via socket
                '''
                if not isinstance(force, np.ndarray) or force.dtype != np.float64:
                    raise ValueError("Forces must be float64 numpy array")
                if not isinstance(vir, np.ndarray) or vir.dtype != np.float64:
                    raise ValueError("Virial must be float64 numpy array")
                if len(force.flatten()) != len(pos.flatten()):
                    raise ValueError("Force array shape mismatch")
                if len(vir.flatten()) != nbeads * 9:
                    raise ValueError("Virial shape mismatch")
                '''
                if shm:
                    timers.start("Notify Force Ready")
                    sock.sendall(Message("FORCEREADY"))
                    timers.stop("Notify Force Ready")
                else:
                    timers.start("Write Force Results")
                    sock.sendall(Message("FORCEREADY"))
                    send_data(sock, pot)
                    send_data(sock, np.int32(nat))
                    send_data(sock, force)
                    send_data(sock, vir)
                    send_data(sock, np.int32(0))
                    timers.stop("Write Force Results")
                
                f_data = False
        
        elif header == Message("EXIT"):
            if rank == 0:
                print("Received exit message from i-PI. Bye bye!")
                if shm and cell_shm is not None:
                    cell_shm.close()
                    icell_shm.close()
                    pot_shm.close()
                    pos_shm.close()
                    f_shm.close()
                    vir_shm.close()
                timers.summary("driver_timing_breakdown")
            elif shm and shm_local_slices and cell_shm is not None:
                cell_shm.close()
                icell_shm.close()
                pot_shm.close()
                pos_shm.close()
                f_shm.close()
                vir_shm.close()
            return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=description)

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    parser.add_argument(
        "--debug", action="store_true", default=False,
        help="Detailed timing output for debugging MPI and socket performance."
    )
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
        "--shm-local-slices",
        action="store_true",
        default=False,
        help="In SHM batch mode on a single node, let each MPI rank read/write its own bead slices directly from shared memory instead of using Scatter/Gather.",
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

    run_driver(
        unix=True if args.shm else args.unix,
        address=args.address,
        port=args.port,
        driver=d_f,
        f_verbose=args.verbose,
        f_debug=args.debug,
        sockets_prefix=args.sockets_prefix,
        shm=args.shm,
        shm_local_slices=args.shm_local_slices,
    )

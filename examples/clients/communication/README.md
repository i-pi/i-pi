Communication transports
=========================

i-PI delegates energy/force/virial evaluation to a *driver* (a client of the
i-PI server). Several transports move the data between server and driver; they
provide identical information and differ only in *how* it is moved. These
examples are organised by transport so each one is easy to find. See the
`_transports` section of `docs/src/distributed.rst` for the full description.

| Subdir    | `<forcefield>` block        | Where the driver runs        | Typical use |
|-----------|-----------------------------|------------------------------|-------------|
| `direct/` | `<ffdirect>`                | in-process (same interpreter)| Python PES from `ipi.pes`; no inter-process communication, lowest overhead |
| `unix/`   | `<ffsocket mode='unix'>`    | same machine                 | fast/empirical potentials run locally; low-latency UNIX-domain socket |
| `inet/`   | `<ffsocket mode='inet'>`    | any machine on the network   | *ab initio* clients, possibly on another host or HPC cluster |
| `shm/`    | `<ffsocket mode='shm'>`     | same machine                 | as `unix`, but the bulk payload travels through shared memory — preferable for large systems |
| `mpi/`    | `<ffmpi>`                   | co-launched MPI ranks        | HPC runs that launch i-PI and the drivers in one MPI job |

`direct/` is special: the potential runs inside i-PI's own Python process, so
there is no driver process and no wire protocol. It is therefore **limited to the
PES modules shipped in `ipi.pes`** (and custom PES files loaded the same way) —
external or compiled codes must use one of the socket/MPI transports. More
in-process / custom-PES examples live in `../pes/`.

Both drivers, one protocol
--------------------------

Every socket example here runs the bundled **Python** driver `i-pi-py_driver`,
but the compiled **Fortran** driver `i-pi-driver` speaks the identical wire
protocol — swap one command for the other and the run is unchanged. The relevant
flags:

- `unix` : `-u -a <name>`            (filesystem rendezvous name)
- `inet` : `-a <host> -p <port>`     (no `-u`)
- `shm`  : `-u -a <name> --shm`      (unix handshake, payload in shared memory)

The Fortran driver is built with `make -C ../../../drivers/f90` (add `MPI=1` for
the `mpi/` examples). Each example ships a `Makefile` (`make` runs it) with the
exact launch command; the Fortran alternative is noted in the Makefile comments.

Examples
--------

- `direct/harmonic/`        — harmonic PES evaluated in-process via `<ffdirect>`
  (no driver process; PES modules only)
- `unix/harmonic/`          — harmonic PES over a UNIX socket (the minimal case)
- `unix/harmonic-batched/`  — same, with `<batch_size>` so one request carries
  several queued structures (here the whole ring polymer)
- `inet/harmonic/`          — harmonic PES over a TCP socket (`address` + `port`)
- `shm/dummy-batched/`      — batched requests over a shared-memory socket, with a
  vectorized no-op PES as a minimal-cost target to gauge communication overhead
- `mpi/`                    — `<ffmpi>` examples (see `mpi/README.md`): a minimal
  harmonic run, two tagged drivers in one launch, and q-TIP4P/f water with a
  dipole `extras` field

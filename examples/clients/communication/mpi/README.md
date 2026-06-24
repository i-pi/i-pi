MPI transport examples (`ffmpi`)
================================

These examples use the `<ffmpi>` forcefield, where i-PI and its drivers exchange
positions/forces over MPI instead of a socket. There is no address or port: i-PI
and the drivers are launched together in a single MPMD job, i-PI being rank 0 of
`MPI_COMM_WORLD` and each driver group a block of the remaining ranks. This avoids
the socket layer entirely and runs over the machine's native interconnect. See
the `_ffmpi-transport` section of `docs/src/distributed.rst` for the full
description.

Each example ships a `run.sh` with the exact launch command and a `Makefile`
(`make` runs `run.sh`).

harmonic/
    The minimal case: i-PI plus four driver ranks evaluating a harmonic PES, in
    one MPMD launch

        mpirun -n 1 i-pi input.xml : -n 4 i-pi-py_driver --mpi -m harmonic -o 1

two-drivers/
    Two independent `<ffmpi>` forcefields (`ffA`, `ffB`) in one launch. Each
    driver group is tagged with `--mpi-name` matching its `<ffmpi name='...'>`,
    and i-PI routes the ranks accordingly; `ffB` is run as several groups with
    different rank counts to show the mechanism. Also combines MPI with
    `<batch_size>`.

water-dipole/
    q-TIP4P/f water PIMD that also retrieves the molecular dipole (an `extras`
    string) from the driver, with all 8 ring-polymer beads bundled into one
    batched request. Uses the compiled Fortran driver. If it was built with MPI
    the run goes over MPI; otherwise `run.sh` falls back to a shared-memory
    socket (`input-shm.xml`), which produces identical output.

Drivers
-------

The Python driver `i-pi-py_driver --mpi` works out of the box once `mpi4py` is
installed. The Fortran driver needs an MPI build:

    make -C ../../../../drivers/f90 MPI=1    # builds drivers/f90/driver.x with --mpi

A driver that itself runs on several MPI ranks (an MPI-parallel code) is launched
with `--group-size N`; i-PI then talks to one root rank per group.

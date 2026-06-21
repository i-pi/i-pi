#!/bin/bash
# q-TIP4P/f water PIMD that also retrieves the molecular dipole (an 'extras'
# string) from the driver.
#
# Uses MPI if the Fortran driver was built with it (make -C drivers/f90 MPI=1),
# otherwise it falls back to a shared-memory socket. Either way all 8 ring-polymer
# beads are bundled into one batched request, and the dipole is written to
# simulation.dipole_* identically for both transports.

DRIVER=${DRIVER:-../../../../drivers/f90/driver.x}   # compiled Fortran driver (q-TIP4P/f lives here)
ADDR=qtip4pf-water

if [ ! -x "$DRIVER" ]; then
    echo "Build the Fortran driver first:"
    echo "    make -C ../../../../drivers/f90          # socket build"
    echo "    make -C ../../../../drivers/f90 MPI=1    # adds the --mpi transport"
    exit 1
fi

if command -v mpirun >/dev/null 2>&1 && nm "$DRIVER" 2>/dev/null | grep -qi mpi_init; then
    echo "# driver built with MPI -> running over MPI (1 batched rank)"
    mpirun -n 1 i-pi input.xml : -n 1 "$DRIVER" --mpi -m qtip4pf
else
    echo "# MPI unavailable -> falling back to a shared-memory socket"
    i-pi input-shm.xml &
    sleep 3
    "$DRIVER" -u -a "$ADDR" -m qtip4pf --shm
    wait
fi

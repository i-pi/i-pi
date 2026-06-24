#!/bin/bash
# Runs i-PI and four MPI driver ranks in a single MPMD launch. i-PI is rank 0
# of MPI_COMM_WORLD; the four driver ranks evaluate a harmonic PES (k=1) and
# exchange positions/forces with i-PI over MPI.
#
# Swap `i-pi-py_driver --mpi -m harmonic -o 1` for the compiled Fortran driver
# `i-pi-driver --mpi -m harm3d -o 1` (build it with `make -C drivers/f90 MPI=1`)
# to check that the Fortran driver interoperates over the identical wire format.

set -e

mpirun -n 1 i-pi input.xml : -n 4 i-pi-py_driver --mpi -m harmonic -o 1

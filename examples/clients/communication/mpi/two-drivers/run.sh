#!/bin/bash
# Two independent MPI driver forcefields (ffA, ffB) in one MPMD launch. Each
# driver group is tagged with --mpi-name matching its <ffmpi name='...'>; i-PI
# routes the ranks to the right forcefield. Both run the harmonic PES here (k=1
# and k=2), but they could be different codes.
# We also run multiple instances of ffB with different number of ranks to 
# provide a full demonstration of the mechanism

mpirun -n 1 i-pi input.xml \
  : -n 1 i-pi-py_driver --mpi --mpi-name ffA -m harmonic -o 1 \
  : -n 1 i-pi-py_driver --mpi --mpi-name ffB -m harmonic -o 2 \
  : -n 2 i-pi-py_driver --mpi --mpi-name ffB -m harmonic -o 2

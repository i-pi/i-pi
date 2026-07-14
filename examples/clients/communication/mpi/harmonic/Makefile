# Makefile for the MPI harmonic oscillator example
#
# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

.PHONY: all clean harmonic
all: harmonic

# i-PI (rank 0) and the driver ranks are launched together in one MPMD job.
harmonic:
	mpirun -n 1 i-pi input.xml : -n 4 i-pi-py_driver --mpi -m harmonic -o 1

clean:
	rm -f *simulation.* RESTART EXIT

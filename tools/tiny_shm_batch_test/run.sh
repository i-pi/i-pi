#!/usr/bin/env bash

# Keep this close to the original profiling run, but tiny.
module load intel/mpi

export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export I_MPI_FABRICS=shm:ofi

source /scratch/c_ccmd/GroupBin/i-pi-stable/env.sh
export PYTHONUNBUFFERED=1

rm -f simulation* RESTART /tmp/ipi_mydriver 2>/dev/null || true

export IPI_SHM_DEBUG=1 
export IPI_SHM_DEBUG_FILE=shmdebug.log
i-pi input.xml | tee ipi.log 2>&1 &

sleep 5

mpirun -np 2 i-pi-py_mpidriver -m dummy -u -a mydriver --shm > driver.log 2>&1

wait 

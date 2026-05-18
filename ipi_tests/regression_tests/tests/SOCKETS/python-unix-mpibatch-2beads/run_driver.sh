#!/usr/bin/env bash
set -e
shift
source /etc/profile >/dev/null 2>&1 || true
module load intel/mpi
export PYTHONPATH=/scratch/c_ccmd/GroupBin/i-pi-stable:
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
mpirun -np 2 /scratch/c_ccmd/GroupBin/i-pi-stable/bin/i-pi-py_mpidriver "$@" 

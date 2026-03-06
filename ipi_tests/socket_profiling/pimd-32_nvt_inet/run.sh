#!/bin/bash 
#SBATCH --job-name=test
#SBATCH --time=1:00:00

#SBATCH --partition=cpu
#SBATCH --mem-per-cpu=2000
#SBATCH -N 1  --ntasks-per-node=34



# Source I-PI env here
source /scratch/c_ccmd/GroupBin/i-pi-stable/env.sh
export PYTHONUNBUFFERED=1



NCORE=1
export MKL_NUM_THREADS=$NCORE
export OPENBLAS_NUM_THREADS=$NCORE
export NUMEXPR_NUM_THREADS=$NCORE
export OMP_NUM_THREADS=$NCORE

# For a detailed timing breakdown
export ENABLE_TIMING_MANAGER=1

NBEADS=32

for i in 8 100 1000 10000
do
    ./run_profiling.sh $i $NBEADS
    wait
    rm simulation* RESTART
    mv timing_breakdown.dat timing_breakdown_${i}.dat
    mv timing_breakdown.log timing_breakdown_${i}.log
done


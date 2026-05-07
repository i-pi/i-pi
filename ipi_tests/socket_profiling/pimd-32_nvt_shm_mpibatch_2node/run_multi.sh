#!/bin/bash
#SBATCH --job-name=test
#SBATCH --time=02:00:00
#SBATCH --partition=cpu
#SBATCH --mem-per-cpu=2000
#SBATCH -N 2  --ntasks-per-node=32
#SBATCH --cpus-per-task=1

source /scratch/c_ccmd/GroupBin/i-pi-stable/env.sh

module load intel/mpi
module load intel/tbb
module load intel/compiler-rt
module load intel/mkl


# For detailed timing analysis
export ENABLE_TIMING_MANAGER=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

export PYTHONUNBUFFERED=1

IPI_INPUT=test-10000.xml
IPI_PORT=50000
NUM_BEADS=32

# -------------------------------
# Determine server and driver nodes
# -------------------------------
ALL_NODES=($(scontrol show hostnames $SLURM_NODELIST))
SERVER_NODE=${ALL_NODES[0]}
DRIVER_NODE=${ALL_NODES[@]:1}

echo "Server node: $SERVER_NODE"
echo "Driver node: $DRIVER_NODE"

# -------------------------------
# Launch i-PI server on server node
# -------------------------------
# explicit binding is important in multinode setup!
srun -N1 -n1 -w $SERVER_NODE --exclusive --cpu-bind=cores  bash  -c "
    source /scratch/c_ccmd/GroupBin/i-pi-stable/env.sh
    hostname > server.host
    i-pi ${IPI_INPUT}
" &> log.ipi &

sleep 50
SERVER_HOST=$(cat server.host)
echo "i-PI server running on $SERVER_HOST"

# -------------------------------
# Launch MPI driver
# -------------------------------
# Source once on driver node
srun -N1 -n1 -w $DRIVER_NODE bash -c "
    source /scratch/c_ccmd/GroupBin/i-pi-stable/env.sh
    sleep 1000
    " &
echo "Sourcing done"
sleep 10


# -------------------------------
# Launch MPI: rank 0 on NODE0, others on NODE1
# -------------------------------
mpiexec \
    -hosts $SERVER_NODE -n 1  i-pi-py_mpidriver -m dummy -u -a mydriver --shm : \
    -hosts $DRIVER_NODE -n 31  i-pi-py_mpidriver -m dummy -u -a mydriver --shm


wait
echo "All i-PI processes finished."

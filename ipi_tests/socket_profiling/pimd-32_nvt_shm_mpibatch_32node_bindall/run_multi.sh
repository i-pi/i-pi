#!/bin/bash
#SBATCH --job-name=test
#SBATCH --time=2:00:00
#SBATCH --partition=cpu --mem-per-cpu=2000
#SBATCH --nodes=32
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1

source /scratch/c_ccmd/GroupBin/i-pi-stable/env.sh

module load intel/mpi
module load intel/tbb
module load intel/compiler-rt
module load intel/mkl

export PYTHONUNBUFFERED=1
export ENABLE_TIMING_MANAGER=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

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
srun -N1 -n1 -w $SERVER_NODE --exclusive --cpu-bind=cores  bash  -c "
    source /scratch/c_ccmd/GroupBin/i-pi-stable/env.sh
    rm /tmp/ipi_*
    hostname > server.host
    i-pi ${IPI_INPUT}
" &> log.ipi &

sleep 50
SERVER_HOST=$(cat server.host)
echo "i-PI server running on $SERVER_HOST"

export I_MPI_PIN=1
export I_MPI_PIN_DOMAIN=core
MPI_CMD="mpiexec"

MPI_CMD+=" -hosts ${ALL_NODES[0]} -n 1 i-pi-py_mpidriver -m dummy -u -a mydriver --shm"

# Remaining ranks: distribute evenly across nodes 1..N
RANK=1
for NODE in "${ALL_NODES[@]:1}"; do
    # this is to keep something alive on each node (not necessary for sourcing env, but it is if you need to mount sth...), your program should run until sleep 1000
    srun -N1 -n1 -w $NODE bash -c "
    source /scratch/c_ccmd/GroupBin/i-pi-stable/env.sh
    sleep 1000
    " &
    # Decide how many ranks per node (here, 1 per node)
    MPI_CMD+=" : -hosts ${NODE} -n 1 i-pi-py_mpidriver -m dummy -u -a mydriver --shm"
    RANK=$((RANK+1))
    # Stop if we reach NUM_BEADS
    if [ $RANK -ge $NUM_BEADS ]; then
        break
    fi
done

sleep 10

echo "$MPI_CMD"

$MPI_CMD
wait

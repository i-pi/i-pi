#!/bin/bash
#SBATCH --job-name=test
#SBATCH --time=02:00:00

#SBATCH --partition=cpu
#SBATCH --mem-per-cpu=2000
#SBATCH -N 2  --ntasks-per-node=32
#SBATCH --cpus-per-task=1

# Load environment
source /scratch/c_ccmd/GroupBin/i-pi-stable/env.sh
export PYTHONUNBUFFERED=1

# For detailed timing analysis
export ENABLE_TIMING_MANAGER=1


NCORE=1
export MKL_NUM_THREADS=$NCORE
export OPENBLAS_NUM_THREADS=$NCORE
export NUMEXPR_NUM_THREADS=$NCORE
export OMP_NUM_THREADS=$NCORE

# Input and port
IPI_INPUT=test-10000.xml
IPI_PORT=50000

# Identify server node (first node in SLURM allocation)
SERVER_NODE=$(scontrol show hostnames $SLURM_NODELIST | head -n1)
echo "Server will run on: $SERVER_NODE"

# -------------------------------
# Launch i-PI server on first node
# -------------------------------
echo "Launching i-PI server..."
srun -N1 -n1 -w $SERVER_NODE --exclusive --cpu-bind=cores bash -c "
    source /scratch/c_ccmd/GroupBin/i-pi-stable/env.sh
    i-pi ${IPI_INPUT}
" &> log.ipi &

# Give server enough time to start and bind the port
sleep 50

NUM_BEADS=32
DRIVERS_PER_NODE=$((NUM_BEADS / 1))

# Get the two nodes for drivers (excluding the server node)
ALL_NODES=($(scontrol show hostnames $SLURM_NODELIST))
SERVER_NODE=${ALL_NODES[0]}
DRIVER_NODES=("${ALL_NODES[@]:1}")  # pick next nodes

echo "Launching $NUM_BEADS driver processes..."
echo "Server node: $SERVER_NODE"
echo "Driver nodes: ${DRIVER_NODES[@]}"

NODE=$DRIVER_NODES

# -------------------------------
# Launch i-PI server on first two nodes
# -------------------------------

echo "Launching drivers on node $NODE..."
srun -N1 -n1  -w $NODE bash -c "
    # Mount once per node
    source /scratch/c_ccmd/GroupBin/i-pi-stable/env.sh

    # Launch multiple driver tasks in background
    for ((i=1;i<=${DRIVERS_PER_NODE};i++)); do
        i-pi-py_driver -m dummy -a ${SERVER_NODE} -p ${IPI_PORT} &
    done

    wait
" &> log.driver.$NODE &

# Wait for all processes to finish
wait
echo "All i-PI processes finished."

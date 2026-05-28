#!/bin/bash
#SBATCH -N 32                
#SBATCH --ntasks-per-node=2  
#SBATCH -t 01:59:00
#SBATCH --mem-per-cpu=2000
#SBATCH --cpus-per-task=1
#

# Load environment
source /scratch/c_ccmd/GroupBin/i-pi-stable/env.sh
export PYTHONUNBUFFERED=1

# For detailed timing breakdown
export ENABLE_TIMING_MANAGER=1

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
sleep 10

# -------------------------------
# Launch i-PI drivers on remaining nodes
# -------------------------------
NUM_BEADS=32
echo "Launching $NUM_BEADS driver processes..."

for nbead in $(seq 1 $NUM_BEADS); do
    sleep 0.2
    srun -N1 -n1  bash -c "
        source /scratch/c_ccmd/GroupBin/i-pi-stable/env.sh
        i-pi-py_driver -m dummy -a ${SERVER_NODE} -p ${IPI_PORT}
    " &> log.driver.$nbead &
done

# Wait for all processes to finish
wait
echo "All i-PI processes finished."

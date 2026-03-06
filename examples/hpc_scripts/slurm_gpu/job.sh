#!/bin/bash -l
#SBATCH -D ./                   # working directory
#SBATCH -J example              # job name
#SBATCH --nodes=1               # number of nodes
#SBATCH --partition=p.ada       # GPU partition
#SBATCH --gres=gpu:a100:4       # 4 GPUs per node
#SBATCH --ntasks-per-node=32    # 32 tasks (clients) per node
#SBATCH --nvmps                 # NVIDIA MPS for multi-process GPU jobs
#SBATCH --time=01:00:00         # max walltime
 
module purge
module load my_modules

# this is optional
# Start nvidia-smi logging every 10 seconds in the background
(while true; do nvidia-smi > gpu.out; sleep 10; done) &

# launch i-PI as usual
python -u i-pi input.xml &

# this will run 32 times launch.sh, each one with a different value of SLURM_LOCALID, ranging from 0 to 31
srun launch.sh &
 
wait
#!/bin/bash
# Gather SLURM info
local_id="${SLURM_LOCALID:-0}"
num_gpus=$(nvidia-smi -L | wc -l)
num_tasks="${SLURM_NTASKS_PER_NODE:-1}"

# Compute tasks per GPU (min 1)
tasks_per_gpu=$(( num_gpus > 0 ? num_tasks / num_gpus : 1 ))
(( tasks_per_gpu < 1 )) && tasks_per_gpu=1

# Compute GPU ID and clamp
gpu_id=$(( local_id / tasks_per_gpu ))
(( gpu_id >= num_gpus )) && gpu_id=$(( num_gpus - 1 ))

# Apply assignment
export CUDA_VISIBLE_DEVICES="$gpu_id"

# run driver
python i-pi-py_driver ...
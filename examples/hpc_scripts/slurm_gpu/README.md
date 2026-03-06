# Torch-Based Client with GPU Oversubscription

Many modern Machine Learning Interatomic Potentials (MLIPs) are implemented in PyTorch, making them capable of running on both CPUs and GPUs. On HPC clusters, however, the hardware configuration often looks like this:

- **HPC clusters**
  - Many slow CPU cores  
  - Few fast GPU devices
- **Simulations**
  - PIMD with many replicas  
  - Systems containing ~1000 atoms  

Running on CPUs typically allows better parallelization, whereas single-point evaluations are intrinsically faster on GPUs. Unfortunately, small systems do not fully utilize GPU capacity, and the limited number of GPUs restricts the number of clients you can run.

For example, if you run 32 replicas, you can usually run at most as many clients as available CPU cores (e.g., 72+ on modern nodes). But if the node has only 1–4 GPUs, you cannot simply assign one GPU per replica.

## Solution: GPU Oversubscription

**GPU oversubscription** allows more than one client to run on the same GPU. This repository includes:

- A sample submission script: `job.sh`
- A launcher script used to configure GPU oversubscription: `launch.sh`

The core idea is to assign `CUDA_VISIBLE_DEVICES` in a round-robin fashion across the available GPUs so that PyTorch runs on the designated device.

## Practical Example

Consider an **i-PI simulation** with:

- 32 beads  
- one PyTorch model per bead  
- a node with **4 GPUs**

In `job.sh`, the line `srun launch.sh` 
is executed 32 times (as defined by `#SBATCH --ntasks-per-node=32`).  
Each invocation of `launch.sh` receives a distinct `SLURM_LOCALID` ranging from 0 to 31.

`launch.sh` uses this ID to:

- assign `CUDA_VISIBLE_DEVICES`  
- optionally modify socket addresses, driver input files, or any other per-client parameters

### Device Mapping Example

| Model (Replica) | `CUDA_VISIBLE_DEVICES` |
|-----------------|------------------------|
| 0               | 0 |
| 1               | 1 |
| 2               | 2 |
| 3               | 3 |
| 4               | 0 |
| 5               | 1 |
| 6               | 2 |
| …               | … |
| 31              | 3 |

This round-robin mapping allows all 32 models to run efficiently across just 4 GPUs.

---

## Optional: GPU Logging
The line `(while true; do nvidia-smi > gpu.out; sleep 10; done) &` in `job.sh` starts a lightweight background loop that records GPU usage every 10 seconds. This is optional and simply helps monitor oversubscription behavior during the run.


Have fun :)

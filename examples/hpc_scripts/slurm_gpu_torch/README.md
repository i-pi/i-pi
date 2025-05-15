Torch-based client with GPU oversubscription
=======================

Nowadays many Machine Learning Interatomic Potentials (MLIP) are based on torch, which can be run on both cpu and gpu.

In cases, the usual situation is the following:
- HPC clusters with:
   - many slow CPU cores
   - few fast GPU cores
- simulations with:
   - PIMD with many replcas
   - too large systems

Usually running on CPU allows a better parallelization, but "single-point" calculation are in principle much faster on GPU, but the limited size of the system does not allow to fully exploit the GPU.

Usually if you are running with 32 replicas, you can instantiate at most a number of clients equal to the number of cores.
If you are running on CPU it is easy to have 72 cores or more, but the usual case is that only few GPU cores (1,2, or 4) are available.

Here is a possible solution: oversubscribe the GPU.
This means allowing more clients to run on the same GPU.

This example provides an example submission `slurm.sh` and a `python` function that you **have to** call before importing `torch`.

Have fun :)
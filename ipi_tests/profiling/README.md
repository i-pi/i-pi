Standardized profiling setups for i-PI
======================================

This folder contains a few inputs designed to stress-test the i-PI implementation, and determine internal bottlenecks.
To emphasize the overhead associated with the Python implementation and the socket communication, these examples 
use dummy drivers that do not compute actual forces, or skip altogether the force evaluation.

The inputs provided are for very small systems with 8 atoms (so the focus is on the constant overhead of the simulation)
but it's possible to perform runs with an arbitrary number of atoms using `run_profiling.sh natoms`.
The -p option activates wall-clock-time profiling using `yappi`.

Examples include:

`classical_md_gas` runs a simple NVE dynamics using an ideal-gas driver

`classical_md_noff` runs a simple NVE dynamics. sets force weights to zero so there is no forcefield evaluation

`remd_scpimd_gas` runs a super-duper calculation with replica exchange, constant-pressure NpT ensemble and Suzuki-Chin high-order path integrals

`remd_scpimd_noff` as above, but without forcefield 

`classical_npt_noff` intermediate setup with classical NPT sampling

`pimd_npt_noff`, `pimd_npt_gas`, `pimd_npt_gas_inet` simple NPT path integral calculations, with different tyles of no-op forcefields

`pimd-32_npt_noff`, `pimd-32_npt_noff-threads` as above, but with 32 beads, with and without multithreading for FFT and PRNG. 
           run setting `OMP_NUM_THREADS` to a value consistent with the thread information in the input.xml file (2 by default)

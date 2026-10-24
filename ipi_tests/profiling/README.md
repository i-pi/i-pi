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

Benchmarks of the internal cost per step
----------------------------------------

`benchmark.py` measures the time i-PI spends per MD step when force evaluation is
negligible. It runs in-process and calls the step function directly, so outputs and
checkpointing are excluded. By default it uses a direct forcefield with the `dummy` PES,
so the full force request/gather pipeline is exercised without sockets; with
`--forcefield none` the forcefield has zero weight and is never evaluated, as in the
`*_noff` inputs, so only the integration of the equations of motion is timed. The default cases cover
NVE, NVT (Langevin and SVR), NPT (isotropic and flexible cell) and path-integral NVT/NPT
with 8 and 32 beads, for 100 to 10000 atoms. Each combination runs in a separate process,
and results are written to a JSON file:

```
python benchmark.py -o new.json                       # all cases and sizes
python benchmark.py -c pimd_npt -n 100,1000 -o a.json  # a subset
python benchmark.py --backend torch --device cuda -o torch.json
python benchmark.py --ipi-root /path/to/other/i-pi -o ref.json
python benchmark.py --compare ref.json new.json        # ratios of the timings
python benchmark.py --compare r1.json,r2.json n1.json,n2.json
```

By default the `ipi` package is taken from the source tree containing the script;
`--ipi-root` points to a different tree, e.g. to compare against another branch.
Timings of large systems are sensitive to other activity on the machine: to detect
small regressions, alternate several runs of the reference and of the new code, and
compare the lists of files, which uses the fastest timing of each combination.

With `--forcefield metatomic` the forces come from a metatomic model, evaluated in
process by a direct forcefield that computes all the beads in one batch. The system is
then a supercell (`-s`) of the 32-water box (`--structure water`, e.g. for PET-MAD) or
of the 32-atom nickel box (`--structure nickel`, for the Lennard-Jones model in
`examples/clients/metatomic/lj`). The model device is independent of the i-PI array
backend, and the time taken by the driver alone on the same structures is also reported,
so that the i-PI overhead is the difference between the two timings:

```
python -c "import upet; upet.save_upet(model='pet-mad', size='xs', version='1.6.0', \
    output='pet-mad-xs-v1.6.0.pt')"
python benchmark.py --forcefield metatomic --model pet-mad-xs-v1.6.0.pt --structure water \
    --model-device cuda -c nve,pimd4_nvt -s 1,2,3 -o numpy.json
python benchmark.py --forcefield metatomic --model pet-mad-xs-v1.6.0.pt --structure water \
    --model-device cuda --backend torch --device cuda -c nve,pimd4_nvt -s 1,2,3 -o torch.json
```

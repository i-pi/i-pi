# MACE execution-path comparison

This folder compares the same ten-step NVE trajectory using:

1. `ffsocket` with ASE's `MACECalculator` in `run_mace.py`;
2. `ffsocket` with i-PI's `bin/i-pi-py_driver`; and
3. the in-process `ffdirect` MACE driver.

From this folder, run:

```bash
./run.sh
```

The script uses the active Python environment for every process. It downloads
`../mace.model` through `../getmodel.sh` if the model is missing. MACE, ASE,
and NumPy must be installed in that environment.

Results are written to `results/`:

- `run.log`, `preflight.log`, and the per-process logs retain diagnostic output.

If a run fails, `results/FAILED` records the failed stage and points to its log.
An existing `results/` directory is preserved with a timestamp before a new run.

The regression test compares coordinates, energies, and conserved quantities
from the three runs with an absolute tolerance of `1e-9`.

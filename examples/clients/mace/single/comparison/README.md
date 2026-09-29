# MACE dynamics and timing comparison

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
NumPy, and Matplotlib must be installed in that environment.

Results are written to `results/`:

- `timings.csv` contains end-to-end wall times, including process startup and
  model loading;
- `trajectory_comparison.txt` reports the numerical trajectory checks;
- `timing.png` compares the three wall times;
- `run.log`, `preflight.log`, and the per-process logs retain diagnostic output.

If a run fails, `results/FAILED` records the failed stage and points to its log.
An existing `results/` directory is preserved with a timestamp before a new run.

Coordinates, energies, and conserved quantities are compared with an absolute
tolerance of `1e-9`. The three paths are expected to be identical at the
precision written to the output files.

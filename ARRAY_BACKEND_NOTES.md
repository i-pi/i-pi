# Array backend (numpy / torch): notes

Notes for the array-backend PR: the changes of behaviour it introduces, and the work left for later PRs.
See also `CONSTRAINTS_REVIEW.md` for constrained dynamics.

## Behaviour changes in this PR (numpy backend)

Bugs found while porting, fixed because they were renames that had not been carried over, obvious typos,
or because they were approved explicitly:

- Normal modes: `nm_rescale.b2tob1` with open paths discarded the open-path transformation (for-else bug,
  since ece43207). Affects ring-polymer contraction combined with open paths.
- Flexible/MTK barostats: when the cell velocity has an all-zero diagonal (e.g. `hfix` on xx, yy, zz) the
  Taylor branch of dt sinh(v dt)/(v dt) missed the factor dt (since 33acc29a). Now a shared helper
  `barostats.dt_sinch`, with a unit test.
- Constrained dynamics could not run from an input file (fixatoms/fixatoms_dof mismatch, removed
  `Forces.forces_mts`), plus small issues in constrtools; see `CONSTRAINTS_REVIEW.md`. New unit test.
- Planetary and Al6xxx KMC could not bind (`Forces.copy` is now `clone`; `Motion.__init__` was not called).
- L-BFGS (`mintools.L_BFGS`), checked against Nocedal & Wright (Alg. 7.4/7.5) and Numerical Recipes
  (dfpmin/lnsrch): the "moving uphill" reset used +g/|g| (uphill) and did not reach the caller; pairs failing
  the curvature condition were stored anyway (indefinite inverse Hessian in later steps) and the direction of
  that step was set to -dx. Now: reset to -g/|g|; such pairs are not stored; the direction is always -H g
  from the stored history. Identical results when all pairs pass the curvature test.
- NEB/string: climb-stage Hessian check (`if` -> `elif`), full-Hessian trimming by `3*natoms` (was `natoms`);
  string mode renames (`fixatoms_dof`, `fixatoms_mask`) and default `tr_trm` (string mode is still disabled
  at input level).
- Instanton: `Fix.get_active_array` condition that was always true; `Fix.get_full_vector` t=3 crash
  (`self.fix.get_mask`); `get_imvector` error message raising TypeError instead of ValueError.
- Phonons: NMFD/eNMFD refinement was overwritten by the unrefined matrix (water/qtip4pf: NMFD eigenvalues now
  match FD to 5e-7, were 75% off).
- IMF restart: `v0` was not set and the basis-size schedule differed from a fresh run.
- Ensemble bias/Hamiltonian weights, PILE time scales, boson indices and cycle energies are stored on the host
  (`host=True`): same values, no device syncs.

## Work after this PR

### Pre-existing bugs found and not fixed (need a decision)
- `GeopMotion.reset` uses the undefined `self.initial_values`.
- NEB/string climb writes the gradient (minus the climbing force) into the force output; maybe intended.
- Instanton friction: without a spectral density, `FrictionMapper.set_fric_spec_dens` crashes
  (`np.arange(omegak.shape)`); `FrictionMapper.evaluate` adds the total friction energy to every bead (counted
  nbeads times in line-search energies). Dead code: `DummyOptimizer.opt_coef`, the `banded = False` branch of
  `LanczosOptimizer`; `Mapper.__call__` passes `new_disc` positionally into `SpringMapper`'s `ret`.
- VSCF with `grid=True` crashes: `scipy.interpolate.interp2d` was removed in scipy 1.14 (the examples run only
  one step).
- SCP/VSCF with `nparallel > 1` clone a 1-bead system to more beads, which raises.
- SCP: `sd` displacement mode unreachable while `precheck` is on; `get_KnD` accumulates `sqtdm` into `sqtDm`
  using `sqtd` (looks like a typo; fixing it changes results).
- `apply_asr` (phonontools, hesstools) uses `eig` on the symmetric inertia tensor; `eigh` would be cleaner.
- Planetary output: `netstring_encoded_savez` writes the bytes of `np.savez` through `io.StringIO` into a
  text-mode file (Python 2 leftover), so saving the frequency matrix fails.
- `Properties.get_vcom(bead=...)` uses an `Atoms` object as momenta, and the bead index arrives as a string.
- Constrained dynamics: see `CONSTRAINTS_REVIEW.md` (Eckart index check, angle and linear-group conditioning,
  tolerance units, solver defaults, batching), and add a regression test with reference outputs.
- Flaky test: `test_softexit.py::test_signal_during_outputs_completes_the_step` failed once in six full runs
  ("release unlocked lock": two signal threads racing in `softexit.cleanup`).
- Failing also on main: regression test `COMMITTEE/XTRA-JSON`; `test_setup.py::test_driver_forcebuild` outside
  tox (driver symlink).

### GPU performance and `torch.compile`
- The cell and barostat state stay on the device, so that whole steps can be compiled. The remaining host
  syncs in step code are graph breaks: the `dt_sinch` branch on the diagonal and `eig` (flexible/MTK
  barostats), `float(K)` in the SVR thermostat, the scalar acceptance logic of stochastic cell rescaling, the
  convergence loop of the constraint solver.
- Fixed per-step cost on CUDA (Quadro P400): ~0.35 ms NVE, ~0.9 ms Langevin, ~2.7 ms NPT, 4-6.5 ms PIMD,
  independent of size (many small kernel launches). NPT on CUDA is slower than torch on CPU: fusing the 3x3
  cell/barostat operations, CUDA graphs or `torch.compile` should remove most of this.
- Forcefields always exchange numpy arrays, so positions and forces are copied to the host at each force
  evaluation. A device-resident path for Python potentials (e.g. ML potentials through `ffdirect`) would avoid
  the copies. `forces.transfer_forces_manual` could also convert its arguments with `to_numpy` itself.
- Host-side code that could move to the device if it becomes a bottleneck: properties estimators that loop over
  atoms (they run on numpy copies), the bosonic exchange recursions (cycle energies and prefix/suffix
  potentials), the constraint solver (batch groups of the same shape).
- Benchmarks to repeat on a compute GPU (the numbers above come from a 2 GB Quadro P400 with slow float64).

### Other backends
- cupy and jax: `xp` comes from `array_api_compat`, the PRNG has one implementation per backend, and all
  writes to the state go through the depend objects (`set`, `__setitem__`, in-place operators), which is where
  jax would need functional updates.

### Testing
- CI has no GPU: the torch job runs on the CPU, with implicit tensor->numpy conversions made to fail (as on
  GPUs). The sync-free CUDA test and CUDA runs must be done by hand on a GPU machine
  (`pytest --ipi-backend torch --ipi-device cuda`).
- The examples are not run with torch in CI (only the unit and regression tests); the benchmark is not run
  in CI.

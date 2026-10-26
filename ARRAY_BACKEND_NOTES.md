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
- SCP and VSCF with `nparallel > 1` failed at bind (a 1-bead system cloned to more beads). SCP now evaluates
  `nparallel` configurations per step on `nparallel` copies of the system; VSCF, which evaluates one
  configuration at a time, ignores `nparallel` with a warning.
- Hamiltonian replica exchange: the Hamiltonian weights of the ensemble were a plain array, separate from the
  depend array that the scaled force components depend on (since 5e863291), so exchanging them did not change
  the bias: the acceptance used the weights before the exchange, and each replica kept its Hamiltonian. They
  are now a depend array, exposed as `hweights`. New unit test.
- Force copies (`dump_state`/`load_state`, `transfer_forces`): forces that were not up to date when the state
  was saved, in practice those of components with zero weight, which are not computed, were loaded as up to
  date: reading such a component afterwards (e.g. `forces_component(k, weighted=False)`) failed, or returned
  the forces of other positions. They are now recomputed when read. Unit test extended.
- Ensemble bias/Hamiltonian weights, PILE time scales, boson indices and cycle energies are stored on the host
  (`host=True`): same values, no device syncs.
- Direct forcefields whose driver sets `backend_arrays` (the metatomic driver) receive positions and cell as
  arrays of the backend, and return energy, forces and virial as arrays of the backend, so with the torch
  backend they stay on the device. Results are unchanged (bit-identical on the CPU for the LJ example model).
  The `dummy` PES does the same, so that direct forcefields with `<pes>dummy</pes>` (e.g. in
  `ipi_tests/profiling/benchmark.py`) exercise this path.
- The state stored before the steps for soft-exit checkpoints (`chk.store()`, every `safe_stride` steps)
  keeps the arrays of the backend on their device, and converts them to numpy only when a checkpoint is
  written, so storing no longer waits for the GPU. Stored scalars of the backend are copies, as the live
  ones can be updated in place.
- Arrays that were built by filling an array in place are computed with whole-array operations, which also
  work with arrays that cannot be modified (JAX): `mtom3`, `f_gather`, `invert_ut3x3`, `get_prop_pq_ms`,
  `kstress_mts`/`stress_mts`, the 3x3 <-> 6 conversions of the flexible barostats (a product with a constant
  selection matrix, which also avoids copying index arrays to the device at each call), the centroid-virial
  kinetic stress estimators. Results change at the level of rounding errors. `Cell.array_pbc` returns the
  wrapped positions instead of modifying its argument.
- Depend objects that become computed after they are created (thermostats of the normal modes, constraints,
  contracted beads, bias scaling) are wired with `dcompute` or `dpipe` instead of setting `_func` by hand.
  `dcompute` raises an `AssertionError` if the object is already computed, synchronized, or a view. `dpipe`
  still replaces an earlier pipe (the time step of the normal modes follows the motion, then each dynamics
  bound to them), but now raises for synchronized objects and views, whose function it used to overwrite.
  Dependencies added after an object is created go through `ddeps`, which raises if the object is not
  computed; the `add_dependency`/`add_dependant` methods are removed.
- Depend objects are stored in attributes with a trailing underscore, `obj.x_` instead of `obj._x`, and
  `dproperties` exposes their value as `obj.x`. Entries of `dproperties` with no depend object behind them
  (they raised a `KeyError` when read) are removed.

The unit tests marked as strict xfail on main for the `b2tob1`, NMFD/eNMFD and `Fix.get_full_vector` (t=3)
bugs pass, and are no longer marked.

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
- SCP: `sd` displacement mode unreachable while `precheck` is on; `get_KnD` accumulates `sqtdm` into `sqtDm`
  using `sqtd` (looks like a typo; fixing it changes results). `ik` with classical statistics gives a NaN
  Hessian when the first displacement is large compared to the sampled width (all reweighting factors vanish).
- `apply_asr` (phonontools, hesstools) uses `eig` on the symmetric inertia tensor; `eigh` would be cleaner.
- Planetary output: `netstring_encoded_savez` writes the bytes of `np.savez` through `io.StringIO` into a
  text-mode file (Python 2 leftover), so saving the frequency matrix fails.
- `Properties.get_vcom(bead=...)` uses an `Atoms` object as momenta, and the bead index arrives as a string.
- Constrained dynamics: see `CONSTRAINTS_REVIEW.md` (Eckart index check, angle and linear-group conditioning,
  tolerance units, solver defaults, batching), and add a regression test with reference outputs.
- Flaky test: `test_softexit.py::test_signal_during_outputs_completes_the_step` failed once in six full runs
  ("release unlocked lock": two signal threads racing in `softexit.cleanup`).
- Failing also on main: `test_metatomic.py` and `test_metatomic_reference.py` error at setup with Python 3.14
  (`AtomisticModel.__annotations__`, TorchScript).
- The other known bugs marked as strict xfail in the unit tests fail in the same way as on main.

### GPU performance and `torch.compile`
- The cell and barostat state stay on the device, so that whole steps can be compiled. The remaining host
  syncs in step code are graph breaks: the `dt_sinch` branch on the diagonal and `eig` (flexible/MTK
  barostats), `float(K)` in the SVR thermostat, the scalar acceptance logic of stochastic cell rescaling, the
  convergence loop of the constraint solver.
- Fixed per-step cost on CUDA (Quadro P400): ~0.35 ms NVE, ~0.9 ms Langevin, ~2.7 ms NPT, 4-6.5 ms PIMD,
  independent of size (many small kernel launches). NPT on CUDA is slower than torch on CPU: fusing the 3x3
  cell/barostat operations, CUDA graphs or `torch.compile` should remove most of this.
- Forcefields exchange numpy arrays, so positions and forces are copied to the host at each force evaluation,
  except for direct forcefields whose driver sets `backend_arrays` (the metatomic driver), which receive and
  return arrays of the backend. Other torch drivers (e.g. `_mace.py`) could opt in the same way.
  `FFRotations` and `FFCommittee` convert the results of such forcefields to numpy, as they combine them with
  numpy; `FFRotations` could be ported to `xp` (rotation matrices on the device), although rotational
  averaging is also being moved into the models. `forces.transfer_forces_manual` could also convert its
  arguments with `to_numpy` itself.
- The metatomic driver synchronizes the device before each batched model call (`_compute_chunk`).
- Storing the state for soft exits walks the whole input tree: ~0.3 ms per step with numpy, ~0.5-1.4 ms
  with torch on CUDA (P400, 100-10000 atoms), more than an NVE step. Runs with cheap forcefields should use
  `safe_stride` > 1.
- Host-side code that could move to the device if it becomes a bottleneck: properties estimators that loop over
  atoms (they run on numpy copies), the bosonic exchange recursions (cycle energies and prefix/suffix
  potentials), the constraint solver (batch groups of the same shape).
- Benchmarks to repeat on a compute GPU (the numbers above come from a 2 GB Quadro P400 with slow float64).
- H100 with metatomic models through `ffdirect` (in-process steps without outputs, ms/step, i-PI overhead
  = step - model): with the torch backend the overhead is ~0.44 ms (NVE) and ~2.2 ms (4-bead PIMD), independent
  of size from 32 to 11000 atoms, i.e. kernel launches. Keeping the forcefield on the device lowered it by
  0.1-0.3 ms at small sizes and by 0.8 ms at 11000 atoms (LJ PIMD 20.3 -> 18.1 ms/step); PET-MAD-xs 4-bead
  PIMD of 96 water atoms went from 20.1 to 17.8 ms/step (numpy backend: 18.3). With the numpy backend and the
  model on the GPU, the overhead of 4-bead PIMD grows to 15 ms (1 BLAS thread) - 30 ms (16 threads) at 11000
  atoms, and was several times larger on a node shared with other jobs.

### Other backends
- cupy and jax: `xp` comes from `array_api_compat`, the PRNG has one implementation per backend, and all
  writes to the state go through the depend objects (`set`, `__setitem__`, in-place operators), which is where
  jax would need functional updates.
- JAX, eager (experiment on `jax-backend`): `set_backend("jax")` enables 64-bit types; depend arrays holding
  JAX arrays become `depend_immutable`, whose writes replace the stored array and whose views (per-bead
  Atoms, thermostat `bind`s, `obj.x_[i] op= v`) read the parent at each access and write through it; masks and
  arrays of indices select copies, as with the other backends.
  `depend_array` writes go through two hooks (`_assign`, `_apply`), the only change for numpy/torch
  (~0.05 us per write, within noise in the benchmarks). Indexing views uses jitted take/put with the index as
  static argument (9-15 us), as eager JAX indexing costs 60-400 us. The PRNG splits a stored key per draw,
  in a single jitted call with the draw (40 us instead of 150 us). Only these array operations are jitted,
  not the algorithms. The MD path builds no arrays in place (see the behaviour changes above).
  Equivalence with numpy (same Gaussian draws, 1e-10): NVE, NVT Langevin/SVR, NPT isotropic/flexible, PIMD
  NVT/NPT; GLE xfails (`stab_cholesky` and the auxiliary momenta `s` are modified in place). Regression tests
  with `--ipi-backend jax`: 47/69 pass (all MD), and the MD restart tests; the others stop at plain-array writes in geometry
  optimizers, NEB, phonons, instantons, exchange, open paths and `stab_cholesky`.
  Eager JAX is slow: ms/step without forces, CPU, numpy / torch / jax, 100 atoms: NVE 0.049 / 0.12 / 0.78,
  NVT 0.11 / 0.29 / 1.8, PIMD-4 0.56 / 1.4 / 11.8, NPT 0.29 / 0.98 / 6.3, NPT flexible 0.93 / 2.6 / 13.4; at
  10000 atoms PIMD-4 28.9 / 26.5 / 34.8. Each eager JAX operation costs 11-16 us (arithmetic), 26 us (reshape), 40-85 us (indexing), 70 us
  (`zeros_like`), against <1 us with numpy; the depend helpers are <10% of the step. Performance needs
  jit of whole steps, which would need the depend state passed in and out of the traced function (as the
  CUDA-graph replay does).

### Testing
- CI has no GPU: the torch job runs on the CPU, with implicit tensor->numpy conversions made to fail (as on
  GPUs). The sync-free CUDA test and CUDA runs must be done by hand on a GPU machine
  (`pytest --ipi-backend torch --ipi-device cuda`).
- The torch CI jobs run the unit tests, and the regression tests and examples with the torch backend
  (`tox -e torch`, `tox -e torch-examples`).
- On CUDA, `test_restart.py::test_restart_matches_regression_reference[instanton]` fails about one time in
  three, with a 3e-10 difference in one element of the final Hessian (non-deterministic GPU arithmetic).

# Array backend (numpy / torch): notes

Notes for the array-backend PR: the changes of behaviour it introduces, and the work left for later PRs.
See also `CONSTRAINTS_REVIEW.md` for constrained dynamics.

## Behaviour changes in this PR (numpy backend)

Bugs found while porting, fixed because they were renames that had not been carried over, obvious typos,
or because they were approved explicitly:

- Normal modes: `nm_rescale.b2tob1` with open paths discarded the open-path transformation (for-else bug,
  since ece43207). Affects ring-polymer contraction combined with open paths.
- Flexible/MTK barostats: when the cell velocity has an all-zero diagonal (e.g. `hfix` on xx, yy, zz) the
  Taylor branch of dt sinh(v dt)/(v dt) missed the factor dt (since 33acc29a). Now computed with
  exp(+-v dt) by `barostats.cell_propagators`, from one exponential of 6x6 block matrices (no branch, no
  `eig`), with a unit test.
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
- The normal-mode transformation is done by default with a matrix product (`transform="matrix"`) instead of
  an FFT: it is faster up to ~128 beads on CPUs (numpy with or without PyFFTW, torch, jax; 1 or 12 threads), and
  on GPUs with fast double precision (H100: 3-6x up to 256 beads x 10000 atoms); the FFT wins from 256-512
  beads on CPUs, and from 32-64 beads on GPUs with slow double precision (Quadro P400). Results change at the
  level of rounding errors (the two agree to ~4e-13).
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
  syncs in step code are graph breaks: `float(K)` in the SVR thermostat, the scalar acceptance logic of stochastic cell rescaling, the
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

### CUDA graphs (`--jit` with torch on CUDA, branch feats/cuda-graphs)
- `ipi/utils/cudagraphs.py` has the generic machinery: `GraphReplay` runs a function eagerly for the first
  calls (two, and at least one in the calling thread, as cuBLAS creates handles per thread), records the next
  one in segments separated by eager functions, and replays them. Code in the recorded function marks the
  points where the recording may be interrupted with `jit.split(obj)`, and calls functions eagerly
  with `jit.eager(func, *args)` (plain calls when not recording); these markers, the recording of the
  operations on the depend objects and `JitError` are in `ipi/utils/jit.py`, shared with jax.jit.
- `ipi/engine/jit.py` decides which steps are replayed (`unsupported`), and `GraphStep` replaces the
  step of the integrator: the recording is interrupted at each force evaluation (`ForceBead.queue` calls
  `split`), with the forces of the component computed eagerly between the segments.
- For non-threaded `FFDirect` forcefields whose driver takes backend arrays (metatomic, dummy), the requests
  and the collection of the results are recorded too, and the segments end at the calls of the driver
  (`eager(FFDirect._results, ...)`, which also applies the offset and parses the extras): on replay the driver is called
  with the request tensors written by the graph, and its tensor results are copied (`torch._foreach_copy_`)
  into those of the recorded call, which the next graph reads; the other items (extras) replace those of the
  recorded result lists, which the depend objects hold. Drivers must not keep their input arrays.
- The depend objects notify a recorder (`depend.set_recorder`) of reads, writes and taints while recording.
  Before a replay, the objects that a segment reads before modifying them are refreshed, tensor scalars that
  were replaced by new objects are copied into the ones read by the graph, and host values must be
  unchanged (`JitError` otherwise, e.g. if the time step changes; the coefficients of the Langevin
  thermostats are tensor scalars, so the temperature can change). After a replay, the taint status
  of all the objects downstream of those modified, the `synchro.manual` of synchronized pairs, and the
  objects held by tensor scalars are restored as at the end of the recorded segment.
- Replayed: `Dynamics` with NVE/NVT/NPT/NST integrators (BAOAB and OBABO, MTS, PIMD) with no, Langevin or
  PILE_L thermostat, and the isotropic (BZP), flexible (MTK) or anisotropic (RGB) barostat with no or a
  Langevin thermostat; not replayed (eager, with a warning): stochastic cell rescaling, SVR, PILE_G, GLE,
  fixed atoms, open paths, bosons, multi-system motions, systems stepped in threads. The barostats read the
  temperature and the pressure as host values, so changing them stops a replayed run with `JitError`;
  the external stress of the anisotropic barostat is a device array, and its changes are replayed.
- NPT without forces (ms/step at 100 / 1000 / 10000 atoms): H100 (job 4665979), eager -> replayed (numpy):
  isotropic 1.42 -> 0.30 / 0.32 / 0.35 (0.17 / 0.34 / 2.35), flexible 3.58 -> 0.81 / 0.75 / 0.75 (0.38 / 0.59
  / 2.69), PIMD-8 isotropic 3.21 -> 0.55 / 0.60 / 0.83 (0.62 / 2.71 / 23.9). The replays are bound by the
  number of kernels (~185 per step for the isotropic barostat, ~400 for the flexible one, of which roughly 160
  are the two matrix exponentials of `cell_propagators`, with 20 Taylor terms and 10 squarings for arguments
  of ~1e-4). On a P400 the 3xN @ Nx3 products of the kinetic stress take ~0.9 ms each (cuBLAS FP64 on a single
  block), and dominate the replays above ~1000 atoms; they take ~10 us on the H100. Compiling (H100, job 4666020):
  `cell_propagators` alone, recorded, takes 85 kernels / 160 us per call eagerly, 65 / 127 us with `matrix_exp`
  compiled, 58 / 117 us compiled whole (the 6x6 matrix products stay cuBLAS calls), 43 / 78 us with 8 Taylor
  terms and 4 squarings (9e-14 from the default series), 28 / 54 us with both; the replayed flexible step goes
  from 0.76 to 0.69 ms with `matrix_exp` compiled, 0.67 ms with `cell_propagators` compiled (bit-identical).
  `cell_propagators` now uses 8 Taylor terms and 4 squarings, and is compiled: replayed flexible step 0.81 / 0.75
  / 0.75 -> 0.60 / 0.53 / 0.54 ms at 100 / 1000 / 10000 atoms, ~400 -> 287 kernels per step (job 4666196). Fusing the whole exponential
  would need a custom kernel, or the series written on the six entries of the upper-triangular matrices
  as scalars (element-wise, so Inductor could fuse it, but hundreds of operations when eager).
- Replayed trajectories and outputs are bit-identical to eager ones (`test_cuda_graphs.py`). The regression
  tests run with graphs (scratch pytest plugin) pass, except the stochastic ones under `--ipi-backend torch`:
  the conftest patch that draws Gaussian numbers with numpy and copies them to the device cannot be recorded.
- P400, 64 atoms, eager -> replayed ms/step: NVE 0.29 -> 0.06, Langevin 0.75 -> 0.19, PIMD-4 2.7 -> 0.66 without
  forces; with the dummy forcefield NVE 0.67 -> 0.44 and PIMD-4 4.0 -> 1.5: the eager force evaluation
  (requests, driver call, gathers) is then the main cost, ~0.4 ms per evaluation for the dummy driver.
- H100 (job 4659635), ms/step without forces, 100-10000 atoms, eager -> replayed (numpy): NVE 0.14 -> 0.05
  (0.03-0.47), Langevin 0.37 -> 0.13 (0.07-1.7), PIMD-4 1.5 -> 0.37-0.46 (0.31-11.2), PIMD-8 2.1 -> 0.57-0.76
  (0.46-20.4): replayed steps are faster than numpy above ~100-1000 atoms. With the nickel LJ metatomic model
  the i-PI part of the step (total - driver alone) goes from 0.40 to 0.25 ms (NVE) and from 2.1 to 0.6 ms
  (PIMD-4) at all sizes (32-10976 atoms); PET-MAD-xs PIMD-4 with 96 atoms 16.3 -> 14.7 ms, the larger PET
  cases are dominated by the model.
- Recording the force pipeline of direct forcefields (P400, dummy driver, 100 atoms, replayed ms/step): NVE
  0.43 -> 0.16 (0.064 without forces, 0.032 for the driver alone), PIMD-4 1.37 -> 0.87, PIMD-8 2.19 -> 1.39;
  a replayed step with the dummy driver makes no syncs, and launches 5 kernels outside the graphs.
  H100 (job 4662134), replayed ms/step before -> after: dummy NVE 0.238 -> 0.085 (0.048 without forces,
  0.017 driver), PIMD-4 (4 driver calls) 0.69 -> 0.51; nickel LJ, i-PI part (total - driver alone, eager /
  before / after): NVE 0.45 / 0.31 / 0.12 ms at all sizes (32-10976 atoms), PIMD-4 2.23 / 0.73 / 0.47 (32
  atoms) to 2.4 / 1.09 / 0.79 (10976 atoms); PET-MAD-xs PIMD-4 with 96 atoms 16.4 / 14.75 / 14.50 ms (driver
  14.1), dominated by the model.
- Kernel fusion: at ~100 atoms the replays are bound by the number of kernels (~1.9 us each on the H100).
  `jit.compiled` compiles a function with torch.compile only when replaying CUDA graphs (plain call
  otherwise, or if torch.compile is unavailable, e.g. torch 2.9 on Python 3.14); its arguments should be
  arrays, as Python scalars are compiled in as constants (hence `ThermoLangevin.T`, `T2`, `S` are backend
  scalars, also making temperature changes replayable). Compiled: the COM removal (`_com_motion`), the
  Langevin step (`_langevin_step`, also used per mode by PILE_L) and the cell propagators of the flexible
  barostats (`cell_propagators`). H100 (jobs 4662227, 4662508),
  classical_nvt_noff 100 atoms: 67 -> 39 kernels, 0.127 -> 0.059 ms/step (numpy 0.074); 1000 atoms 0.121 ->
  0.072 (numpy 0.217); PIMD-4/8 without forces, 100 atoms, 0.374 -> 0.234 / 0.296 ms (numpy 0.311 / 0.464).
  Replays then match eager steps to rounding errors (fused reductions), so the tests use a tolerance.
- Normal-mode transforms on the H100 (job 4666371), round trip `b2nm` + `nm2b` recorded as a CUDA graph, us
  (kernels): FFT 30-38 (15) at 100-1000 atoms and 4-256 beads, 36-718 at 10000 atoms; compiled FFT ~25% faster
  (10 kernels); matrix 6-24 (2) at 100-1000 atoms, 9-185 at 10000 atoms, unchanged by compiling. The matrix
  product is 3-6x faster than the FFT at all the sizes tested (up to 256 beads x 10000 atoms), and agrees with it
  to 4e-13. With numpy (no PyFFTW, 16 cores) the matrix product was also faster at all these sizes.
  A PIMD NVT step with PILE_L does 7 transforms (6 without forces): the momenta go from beads to normal modes
  for each thermostat step and for the free ring-polymer step, and back for each COM removal and force kick.
  `transform='matrix'` in PIMD NVT without forces (job 4666942), ms/step fft -> matrix at 100/10000 atoms:
  replayed 4 beads 0.23 -> 0.16 / 0.31 -> 0.23, 8 beads 0.29 -> 0.22 / 0.47 -> 0.35, 32 beads 0.63 -> 0.55 /
  1.53 -> 1.28; eager 1.48 -> 1.07 (4 beads), 2.01 -> 1.61 (8), 5.26 -> 4.94 (32); numpy 8 beads 0.47 -> 0.33 /
  20.2 -> 15.2, 32 beads 1.35 -> 1.14 / 84.5 -> 65.8. With PET-MAD-xs (water, 384 and 768 atoms, 4 and 8
  beads, 25-74 ms/step) no difference within noise.
- Possible extensions: stochastic cell rescaling (host logic and `float` syncs), SVR without `float(K)`, PILE_G, evaluating the forces of all the components of a level at once
  (now each component is computed when its first bead is queued), cheaper staging of the inputs (host
  arrays such as `nmts` and `pdt` are compared with `np.array_equal` at each replay).

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
- JAX, jit of the MD steps (`--jit` with JAX, branch `feat/jax-jit`, CUDA graphs rebased on `feat/jax-backend`):
  `ipi.utils.jaxjit.JitReplay` uses the recording of the CUDA graphs (`_Recorder`, the depend hooks, the
  restoration of the taint status, the `split`/`eager` markers, `unsupported` and the warmup steps). The
  second eager step is observed to find the depend objects read before being written; the third is traced
  with these (and the PRNG key) as arguments of `jax.jit`, returning the values of the objects it writes,
  compiled ahead of time (`lower().compile()`, so it never retraces), and replayed: inputs brought up to
  date with `get()`, host inputs (`dt`, `nmts`, `splitting`, ...) compared with their traced values
  (`JitError` if they change), outputs stored, statuses restored. depend changes: the root `_value` of
  `depend_immutable` notifies reads (in-place operators and writes through views read without `get`), and
  full assignments replace the stored array without reading it. Unlike CUDA graphs, a trace cannot be
  split: the driver calls marked with `eager` are traced, so only direct forcefields whose potential works
  on JAX arrays (`dummy`) are compiled; other forcefields leave the steps eager, with a warning (they would
  need `jax.pure_callback`, i.e. a host round trip). Same trajectories as eager JAX to ~1e-15, same keys:
  NVE, NVT Langevin, PIMD PILE_L, NPT isotropic and flexible. CPU, 100 atoms, dummy forcefield, ms/step
  numpy / jax eager / jax jit: NVE 0.14 / 1.8 / 0.078, NVT 0.22 / 3.0 / 0.19, PIMD-4 0.93 / 16 / 0.78; at
  1000 atoms jit is slower than numpy with thermostats (NVT 0.43 vs 0.86, PIMD-4 2.6 vs 4.9), partly as
  XLA draws float64 Gaussian numbers 4x slower than numpy on the CPU (160 vs 40 us for 3000).
  H100 (kuma jobs 4679488, 4679541; Python 3.11, torch 2.12.1, JAX 0.10.2), ms/step without forces / with the
  dummy forcefield, 100 -> 10000 atoms: NVE graphs 0.031-0.036 / 0.078-0.084, jit 0.085-0.082 / 0.102-0.104;
  PIMD-8 (PILE_L) graphs 0.30-0.48 / 0.51-0.85, jit 0.30-0.35 / 0.43-0.48 (numpy 0.47-19.5 / 0.76-22.6, eager
  JAX 19-27). jit is 13-78x faster than eager JAX; CUDA graphs have the lower fixed cost (classical MD), jit
  scales better with beads and atoms (one XLA program, forces included). Jit PIMD is slower at 3162 atoms than
  at 10000 (not investigated). Jit vs eager JAX on the GPU: 2e-16 to 3e-15 after 10 steps, same keys.
  Full tables: report.md (untracked, repository root).
  torch.compile and jax.jit do not overlap: with torch the steps are CUDA graphs and only the functions
  marked `@compiled` are compiled (whole-step torch.compile breaks on the depend locks and the lazy forces);
  with JAX the whole step is one XLA program and `@compiled` functions are traced inline.
- One option, `--jit` (`set_backend(..., jit=True)`, `InteractiveSimulation(jit=True)`), compiles the steps
  with either backend: CUDA graphs with torch on a CUDA device, jax.jit with JAX (it was `--cuda-graphs` for
  torch); with numpy, or torch on the CPU, it is an error.

### Testing
- CI has no GPU: the torch job runs on the CPU, with implicit tensor->numpy conversions made to fail (as on
  GPUs). The sync-free CUDA test and CUDA runs must be done by hand on a GPU machine
  (`pytest --ipi-backend torch --ipi-device cuda`).
- The torch CI jobs run the unit tests, and the regression tests and examples with the torch backend
  (`tox -e torch`, `tox -e torch-examples`).
- On CUDA, `test_restart.py::test_restart_matches_regression_reference[instanton]` fails about one time in
  three, with a 3e-10 difference in one element of the final Hessian (non-deterministic GPU arithmetic).

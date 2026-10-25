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

### CUDA graphs (`--cuda-graphs`, branch feats/cuda-graphs)
- `ipi/utils/cudagraphs.py` has the generic machinery: `GraphReplay` runs a function eagerly for the first
  calls (two, and at least one in the calling thread, as cuBLAS creates handles per thread), records the next
  one in segments separated by eager functions, and replays them. Code in the recorded function marks the
  points where the recording may be interrupted with `cudagraphs.split(obj)`, and calls functions eagerly
  with `cudagraphs.eager(func, *args)` (plain calls when not recording).
- `ipi/engine/cudagraphs.py` decides which steps are replayed (`unsupported`), and `GraphStep` replaces the
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
  unchanged (`CUDAGraphError` otherwise, e.g. if the temperature changes). After a replay, the taint status
  of all the objects downstream of those modified, the `synchro.manual` of synchronized pairs, and the
  objects held by tensor scalars are restored as at the end of the recorded segment.
- Replayed: `Dynamics` with NVE/NVT integrators (BAOAB and OBABO, MTS, PIMD) with no, Langevin or PILE_L
  thermostat; not replayed (eager, with a warning): barostats, SVR, PILE_G, GLE, fixed atoms, open paths,
  bosons, multi-system motions, systems stepped in threads.
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
- Possible extensions: isotropic barostat (needs the tensor scalars staged as above, verified to replay
  exactly), SVR without `float(K)`, PILE_G, evaluating the forces of all the components of a level at once
  (now each component is computed when its first bead is queued), cheaper staging of the inputs (host
  arrays such as `nmts` and `pdt` are compared with `np.array_equal` at each replay).

### Other backends
- cupy and jax: `xp` comes from `array_api_compat`, the PRNG has one implementation per backend, and all
  writes to the state go through the depend objects (`set`, `__setitem__`, in-place operators), which is where
  jax would need functional updates.

### Testing
- CI has no GPU: the torch job runs on the CPU, with implicit tensor->numpy conversions made to fail (as on
  GPUs). The sync-free CUDA test and CUDA runs must be done by hand on a GPU machine
  (`pytest --ipi-backend torch --ipi-device cuda`).
- The torch CI jobs run the unit tests, and the regression tests and examples with the torch backend
  (`tox -e torch`, `tox -e torch-examples`).
- On CUDA, `test_restart.py::test_restart_matches_regression_reference[instanton]` fails about one time in
  three, with a 3e-10 difference in one element of the final Hessian (non-deterministic GPU arithmetic).

# i-PI next: development notes

Development notes for `next`, the branch that will become i-PI 4.0: the changes of behaviour relative to `main`
(i-PI 3.x), the state of the array backends (numpy, torch, JAX) and of the compiled MD steps (`--jit`), and the
work left. Update them with every change that affects results, performance or the developer API.

## Behaviour changes relative to main (numpy backend)

Bugs found while porting, fixed because they were renames that had not been carried over, obvious typos,
or because they were approved explicitly:

- Normal modes: `nm_rescale.b2tob1` with open paths discarded the open-path transformation (for-else bug,
  since ece43207). Affects ring-polymer contraction combined with open paths.
- Normal modes: the dynamical mass factors (`manual`, `pa-cmd` and `wmax-cmd` frequencies) did not depend
  on the frequencies of the free ring polymer, so they kept the values of the initial temperature when it
  changed, e.g. with replica exchange. Also on main. New unit test.
- Instantons: the symmetric band given to `scipy.linalg.solve_banded` (`instools.sym_band`) left the corners
  outside the band uninitialized, so the solver failed at random when they held infinities or NaNs, which
  it checks for. Also on main. New unit test.
- NEB and string method: a Hessian given for all the coordinates of the path (or of the climbing image) raised
  an AttributeError, as the mask of the moving coordinates was read before the gradient mappers were bound.
  Also on main. New unit test (NEB; the string method is disabled at input level).
- Flexible/MTK barostats: when the cell velocity has an all-zero diagonal (e.g. `hfix` on xx, yy, zz) the
  Taylor branch of dt sinh(v dt)/(v dt) missed the factor dt (since 33acc29a). Now computed with
  exp(+-v dt) by `barostats.cell_propagators`, from one exponential of 6x6 block matrices (no branch, no
  `eig`), with a unit test.
- Constrained dynamics could not run from an input file (fixatoms/fixatoms_dof mismatch, removed
  `Forces.forces_mts`), plus small issues in constrtools; see [Constrained dynamics](#constrained-dynamics). New unit test.
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
  (`self.fix.get_mask`), and t=4/5 (friction tensors with fixed atoms) crash, as the indices of the inner loops
  were not reset; `get_imvector` error message raising TypeError instead of ValueError.
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
  Also, the integrators skipped the bias (`ensemble.has_bias` was false) whenever the Hamiltonian weights
  differed from 1 (since 5e863291), so the replicas evolved with the unscaled Hamiltonian, and with no explicit
  bias, while the conserved quantity, the exchange acceptance and the barostat included them. This changes
  the results of Hamiltonian replica exchange runs. New unit test.
- Combining the extras of several force components appended the raw strings of the other components to
  those of the first component (since cf4419b5), so `extras_component_raw(0)` printed the strings of all the
  components once the combined extras had been computed. Unit test extended.
- PLUMED (`ffplumed`): when PLUMED biases the potential energy (metadynamics with `use_energy`), the virial
  it returned contained the physical virial twice instead of only the bias contribution, and `mtd_update`
  failed when the positions had moved since the last evaluation (a plain dict instead of a `ForceRequest`).
  Both were strict xfails on main. The initial structure, the physical forces and the positions given to
  `mtd_update` are converted to numpy at this boundary, so PLUMED also works with the torch backend (and
  again with numpy, which failed on the read-only masses of the depend properties).
- Missing dependencies of the forces: the fourth-order Suzuki-Chin term on alpha, the total SC virial on the
  per-bead ones (it depended on the SC potential), and the combined extras and the components added to a bias
  (`add_component`) on the component weights; they gave stale values when these changed (e.g. the extras of a
  bias whose weights are exchanged by replica exchange). New unit tests.
- A forcefield that returns a non-symmetric virial gets a warning, once per forcefield: as on main, only the
  upper triangle of the virial is used.
- Force copies (`dump_state`/`load_state`, `transfer_forces`): forces that were not up to date when the state
  was saved, in practice those of components with zero weight, which are not computed, were loaded as up to
  date: reading such a component afterwards (e.g. `forces_component(k, weighted=False)`) failed, or returned
  the forces of other positions. They are now recomputed when read. Unit test extended.
- SVR thermostat (also the centroid of PILE_G): the sign of the rescaling factor used twice the kinetic energy
  relative to the construction of Bussi, Donadio and Parrinello (JCP 126, 014101, App. A), from which its
  square is computed (since 2.0), so some of the momentum reversals were missed. The sampled distribution was
  correct; the dynamics change only with few degrees of freedom and tau not much larger than the time step
  (with 30 degrees of freedom and dt = tau, 1e-4 of the steps reverse the momenta). With one degree of
  freedom, the step is now the Langevin propagator. New unit test.
- Socket interface: a client that disconnected during the handshake was added to the pool (`status | Up`);
  request ids were compared by identity, so in `lock` mode a client could not get its replica back when the id
  was a large or recreated integer; with `consolidate_messages=False`, the request of a client that timed out
  was never queued again; a driver closing the connection before the force header made `getforce` loop
  forever; `InvalidStatus` raised a `TypeError` instead; `close` did not close the sockets whose `shutdown`
  failed. The batched `gas` and `harmonic` drivers of `FFDirect` evaluate structures with different numbers
  of atoms one at a time (they assumed equal sizes). All were strict xfails of the unit tests.
- `logsumlog` (sums of signed numbers stored as logarithms) gave the wrong result when the larger term was
  negative (-3 + 2 gave -5), which affected the isotope fractionation estimators (`isotope_scfep`,
  `isotope_tdfep`) when the centroid-virial kinetic term is negative for some configurations. `apply_asr`
  (phonons, SCP, `hesstools`) diagonalized the inertia tensor with `eig`, which gives non-orthogonal
  eigenvectors when two moments of inertia are equal (symmetric and spherical tops), so the removal of the
  rotations (`asr="poly"`) was not a projection; it uses `eigh`. `Cqp` (GLE tools) rescaled the drift matrix
  the wrong way, giving wrong <q^2> and <qs> for frequencies other than 1 (`gle_frequency_kernel` only uses the
  momentum column, which was right). The unused, and incorrect, `exp_ut3x3` and `eigensystem_ut3x3` are
  removed. Strict xfails of the unit tests, plus a new test for the symmetric top.
- Input and output: the planetary model could not save its frequency matrix (`netstring_encoded_savez` wrote
  the bytes of `np.savez` through a `StringIO`); the NPZ data is now stored as base64, which is safe in the
  text file. Property strings with an unmatched unit brace (`potential{electronvolt`, `potential}`) raise a
  `ValueError` instead of giving a wrong name or unit; an empty initialization file raises a `ValueError`
  (it was an `UnboundLocalError`). The dipoles of driven dynamics read from the legacy "Total dipole moment"
  strings were NaN (strings multiplied by a float), and only the first bead was read. The interpolation of
  instantons left `new_hessian.dat` open, so it could be incomplete when the function returned. Strict xfails
  of the unit tests.
- Engine and drivers: `GeopMotion.reset` with `bfgstrm` failed (it read the undefined `initial_values`), which
  stopped the Al6xxx KMC with that optimizer; it restores the initial trust radius, whose default (0.4) is set
  in `GeopMotion`. `reset` also left a converged optimizer converged, so that the KMC skipped every relaxation
  after the first converged one (also on main); it now clears the convergence, resets the Hessian of
  `damped_bfgs` and the adaptive step of the line search of `sd` and `cg`, and no longer zeroes the
  `corrections` and `scale` settings of L-BFGS (which then went to the restart). The spring energy and forces
  of the beads (`vpath`, `fpath`), which were only used in tests, are removed (`nm.vspring` and `nm.fspring`
  hold those of the path integral): they were wrong for open paths and bosons, and `fpath` did not depend on
  the masses. `beads.kin` failed (`kins` was not a property); `vcom(bead=...)` failed. The GLE state given
  in `<initialize><gle>` was looked for in the ensemble instead of the motion, so it could never be read. The
  temperature of a single normal mode (`temperature(nm=...)`) with a fixed centroid counted the
  fixed coordinates as nbeads-1 degrees of freedom in every mode. The cavity forcefield (`FFCavPhSocket`)
  passed its error messages as the soft-exit status, so it raised a `ValueError` instead of stopping, and
  read the dipole data before checking it was there. The ASE driver without stress rejected its own zero
  stress. The double-well drivers with an explicit bath (`DW_bath`, `DDW_with_explicit_bath`) could not run
  since the bath was converted to `compute_structure`, and returned their energy as an array, which numpy 2.5
  does not convert with `float()`; `DW_bath` built its bath again at every call. The force of the spline
  driver (`pes/spline.py`) was the derivative of the tabulated energy with respect to the scaled coordinate,
  missing the scaling factor (`factor_coord = 5`), so it was 5 times smaller than the derivative of the energy
  it reports: results with that driver change. A stuck FIRE minimization passed its message as the soft-exit
  status, so it raised "Unknown option for softexit status" instead of stopping; the status and message of
  `softexit.trigger` and `cleanup` are keyword-only, and the errors that a restart cannot resolve (also a
  request cancelled while a bead force waits for it) have status `bad`. Strict xfails of the unit tests, plus
  a test of `reset` and of the DDW driver.
- Tools: `neb_interpolate` called `Atoms.copy` (it is `clone`) and divided the read-only masses in place;
  `mux-positions` with wrapping or unwrapping set a float shape; `install_driver` ignored the exit status of
  git and make, so a failed compilation was reported as a failed copy; the help of `i-pi -S` described
  profiler files instead of the prefix of the sockets.
- Atom swaps (`atomswap`) with more than one bead: the swapped atoms keep their masses, so the spring energy
  of the ring polymers changes when their masses differ; the acceptance and the conserved quantity left this
  change out (since the motion was written), so path-integral simulations did not sample the right
  distribution. It is computed from the normal modes, as in the alchemical exchange. Results of PIMD with
  atom swaps change; classical ones do not. Strict xfail of the unit tests, extended to an accepted swap.
  Atom swaps and alchemical exchanges stop at bind when an exchanged atom has an open path or is a boson:
  their change of the spring energy is that of closed paths, while open paths use `o_omegak`, and the spring
  energy of bosons is the exchange potential, which does not split over atoms (the conserved quantity drifted,
  and alchemical exchanges of bosons failed).
- Hamiltonian replica exchange: the scaled force components of the bias scale the extras that are interpolated
  as physical quantities (`interpolate_extras`), as the energy; they entered the bias unscaled.
- Geometry optimization with `damped_bfgs` and fixed coordinates (`fixatoms_dof`): the step was limited by the
  displacements of triplets of free coordinates, which failed when their number is not a multiple of 3, and
  mixed the coordinates of different atoms otherwise. The displacements of the atoms are now computed from
  the free coordinates. Strict xfail of the unit tests, plus a test of the step measure. NEB and the string
  method (also when climbing) did not give the free coordinates to the step measure, so they still had the
  problem; `Damped_BFGS` now requires them. Also with fixed coordinates: conjugate gradient built its
  directions from all the forces, also those on the fixed atoms, which turned it into steepest descent when
  they were large (64 steps instead of 3 in the test); with more than one bead, BFGS, BFGSTRM, L-BFGS and
  damped BFGS failed, as their masks of the Hessian and of the L-BFGS history covered one bead. The
  optimizers of `GeopMotion` see only the free coordinates, through one code path for all systems (the same
  results without fixed atoms, to rounding). New unit tests of NEB with a fixed coordinate, of conjugate
  gradient with a fixed atom, and of optimizations of two beads. At an exact stationary point (zero
  gradient) the damped BFGS update made the Hessian NaN; it now takes no step. New unit test. In geometry
  optimization, `damped_bfgs` kept its Hessian, a direct one, in the inverse Hessian of BFGS (stored as
  `invhessian_bfgs`); it is now kept and stored as the Hessian of BFGSTRM (`hessian_trm`), as NEB and the
  string method keep theirs in `hessian`, so restarts of `damped_bfgs` optimizations written before start
  from the identity. New unit test of the restart of the Hessians.
- IMF and VSCF (`normalmodes`), checked against Kapil, Engel, Rossi, Ceriotti, JCTC 15, 5845 (2019), which
  describes this implementation:
  - VSCF with `grid=True` failed with scipy >= 1.14 (`interp2d` was removed); the pair couplings are
    interpolated with `RectBivariateSpline`, the same bicubic interpolation.
  - In the SCF, the thermal density of each mode was normalised after adding each state instead of once,
    which weights the excited states wrongly (it matters for modes with hbar w below ~2 kT); mixing the mean
    field (`alpha` < 1) converged to v_indep/alpha + mean field; the one-mode potentials used in the pair maps
    were taken by the index of the mode minus the number of zero modes, wrong when a list of modes is given
    (`modes.dat`) or when low-frequency modes are left out.
  - The VSCF free energy is the variational (Gibbs-Bogoliubov) free energy of the product of the thermal
    densities, sum_i A_i - sum_{i<j} <<V2_ij>>. The sum of the free energies of the modes in their mean
    fields (Eq. 19 of the paper, used so far) counts each pair coupling twice: for weak couplings the VSCF
    correction to the IMF free energy is about half of what it was. VSCF results change.
  - The free energy of each mode is converged to `athresh`, as documented (the change was divided by the
    number of modes), and the internal energies use Boltzmann weights relative to the ground state (they were
    NaN below a few K). The IMF free energy printed by VSCF includes the potential offset, as in IMF (once
    per primitive cell: with `nprim > 1` VSCF divided it by `nprim` twice).
  - The IMF and VSCF free energies are stored as `total_imf_free_energy` and `total_vscf_free_energy`.
  New tests on separable and coupled Morse potentials against exact references, and of the pair coupling of
  a polynomial.
- Instantons with fixed atoms (`fixatoms`): `hessian_update=recompute` and the final Hessian with friction
  expanded the Hessians of `get_hessian`, which are already full-size, once more, and failed; the NR optimizer
  combined the masses and spring forces of all the atoms with the Hessian of the active ones, and failed;
  L-BFGS added the active forces to the full search direction and passed the active positions to the mapper,
  which took all of them, and failed. The instanton `Mapper` now takes the positions of the active atoms, as
  the mappers of the other optimizers do, and NR uses its forces, as Lanczos does. `fixatoms_dof` that are
  not whole atoms, which were fixed as whole atoms when their number was a multiple of 3, are rejected. New
  unit tests of all the optimizers with a fixed atom.
- Instanton rates with friction, checked against Litman et al., JCP 156, 194106 (2022): without a tabulated
  Laplace transform of the friction kernel (`fric_spec_dens`), i.e. with frequency-independent friction, the
  setup failed; the friction energy of the ring polymer (Eq. 20) was added to the energy of every bead, so the
  energies summed by the optimizers counted it nbeads times while the gradient counted it once; the Laplace
  transform was not divided by its value at `fric_spec_dens_ener` (compared and evaluated in cm^-1 against
  frequencies in atomic units since d4dd9d70), which changes the results of the runs that set it. Removed
  `DummyOptimizer.opt_coef` and the unreachable full-matrix branch of the Lanczos step; `Mapper` passes
  `new_disc` to `SpringMapper` by keyword (it went into `ret`). New unit tests.
  With `frictionSD=False` (a constant friction tensor) the first step failed when writing the friction
  Hessian, which only exists with `frictionSD`; with `frictionSD` and fixed atoms, the friction Hessian
  integrated the friction along the path of all the atoms with the friction tensor of the active ones, and
  failed. New unit test, comparing one atom with and without a second, fixed, atom.
- Self-consistent phonons: with the default `precheck`, the `sd` displacement mode never moved the mean position
  nor updated the Hessian (its branch was an `elif` of the check); the reweighting factors of the samples
  underflowed when the mean position moved far compared to their width (e.g. `ik` with classical statistics),
  which made the force and the Hessian NaN; they are now computed relative to the largest one, which leaves
  their ratios, the only thing used, unchanged. The weights of the batches of samples, exp(-var(log w)),
  underflowed in the same way (e.g. starting 10% farther from the minimum in the tests), giving a NaN Hessian:
  the batches are averaged with weights relative to the largest one (the absolute ones are kept for the
  threshold), and those whose relative weight underflows, whose variance may be undefined, are left out of the
  error of the force, which was NaN. The duplicated displacement matrices (`tDm`, `itDm`, `iDm`,
  `sqtDm`), equal to the others and never used, are removed (the `sqtdm` "typo" had no effect). New tests,
  including the convergence to the optimal Gaussian of a quartic double well, known in closed form for
  classical statistics and from a one-dimensional equation for quantum statistics.
- Climbing-image NEB (and the string method): after each climbing step, the forces of the beads are stored so
  that they are not recomputed; for the climbing bead this stored the gradient of the climbing force (opposite
  sign, and inverted component along the path) instead of its physical force, which reached the force outputs,
  everything reading the forces afterwards, and the forces saved in restarts. New unit test.
- Fourth-order Suzuki-Chin path integrals (production runs always have alpha = 0, which is not an input):
  - The SC virial missed the force-product term of the strain derivative of |f|^2 / m, 2 (f/m) x f (Kapil et
    al., JCTC 15, 3237 (2019), Eq. 9; checked by strain finite differences): SC NPT and the `pressure_tdsc` and
    `vir_tdsc` outputs change (the missing term is O(1/P^2), so the pressure converged only at second order).
    SC NVT trajectories are unchanged.
  - With contracted components (RPC) and centered differences, the + and - displaced odd beads were put on
    one auxiliary ring, which the contraction mixed: the fourth-order term of contracted components was wrong.
    They displace the whole ring, as for alpha != 0 (whose limit alpha -> 0 this matches).
  - The scaled-coordinate estimators (`sc_scaledcoords`) kept the coefficient of |f|^2, proportional to beta^2,
    fixed when scaling beta (Yamamoto, JCP 123, 104101 (2005), Eqs. 33-34): on a harmonic oscillator they gave
    1.421 instead of 1.500 for the energy, and a heat capacity 3.6 times too large.
  - With alpha = 0, mixing components with forward and centered differences failed, as they shared one
    auxiliary system: there is one per number of beads.
  - The help of `sc_op_scaledcoords` gives the operator heat capacity, beta^2 kB (Cov(eps_v, eps_v^TD) +
    <eps_v'>) (Kapil et al., Eq. 13): note the covariance with the thermodynamic estimator, and the sign of
    eps_v', opposite to `sc_scaledcoords`.
  Strict xfails of the unit tests (the forward-difference virial cases use a smaller displacement), a test of
  mixed differences, and an exact test of all the SC energy estimators and heat capacities on a harmonic
  oscillator (Gaussian ring polymer: the thermodynamic ones reproduce the derivatives of ln Z_SC, the operator
  ones the quantum result to fourth order). No regression test uses Suzuki-Chin.
- The normal-mode transformation is done by default with a matrix product (`transform="matrix"`) instead of
  an FFT: it is faster up to ~128 beads on CPUs (numpy with or without PyFFTW, torch, jax; 1 or 12 threads), and
  on GPUs with fast double precision (H100: 3-6x up to 256 beads x 10000 atoms); the FFT wins from 256-512
  beads on CPUs, and from 32-64 beads on GPUs with slow double precision (Quadro P400). Results change at the
  level of rounding errors (the two agree to ~4e-13).
- When the momenta are held in normal modes (after the PILE thermostats), the motion of the centre of mass
  is removed from the centroid normal mode, which is the only one it changes (also with open paths), instead
  of transforming the momenta to beads and back: a PIMD NVT step with PILE_L does 4 normal-mode
  transformations instead of 6 (numpy PIMD-4/8, 1000 atoms: 14-17% faster). Results change at the level of
  rounding errors. The free ring-polymer step is `@compiled` (fused when the steps are replayed as CUDA graphs).
- With one bead, the positions are propagated in the bead representation (`qcstep`), which avoids copying
  positions and momenta between beads and normal modes at each step (identical results; numpy NVE 6-20%
  faster, NVT Langevin 2-11%, torch eager NVE 21% on the P400).
- The PILE_L and PILE_G thermostats use one Langevin thermostat for all the normal modes (all but the
  centroid for PILE_G), with one time scale and temperature per mode, instead of one thermostat per mode: the
  noise of all the modes is drawn at once, and a PIMD step launches the same number of kernels whatever the
  number of beads (Quadro P400, 100 atoms, no forces: 78 instead of 156 / 260 for 4 / 8 beads; `--jit`
  0.48 -> 0.28 / 0.88 -> 0.33 ms/step; numpy 17-35% faster at 100 atoms). `T_on_sm` has the shape of the
  momenta, `T2` and `S` one value per mode. With numpy the trajectories are the same as before (only the
  thermostat energy changes, at the level of rounding errors); with torch and JAX the random numbers change.
  The thermostat energy is no longer kept per mode: a stored value goes to the single thermostat on restart.
- The SVR thermostat (also on the centroid of PILE_G) runs on the device: the kinetic energy, the sign of the
  scaling and the random numbers stay on the device (the gamma variate from `Random.gammavec`, with
  `torch._standard_gamma` and the torch generator, or `jax.random.gamma`), so that its steps can be compiled
  with `--jit`. With numpy the random numbers are the same as before (also in order); with torch and JAX they
  come from the device generator, so trajectories change for a given seed. The thermostat no longer skips the
  draws when the kinetic energy is zero.
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
  kinetic stress estimators, the fourth-order Suzuki-Chin finite differences and SC kinetic stress, the
  centroid-virial estimators `kinetic_cv(atom)`, `kinetic_ij`, `kinetic_tens` and the `kinetic_cv`,
  `kinetic_od` and `r_gyration` trajectories. Results change at the level of rounding errors (the SC
  trajectories and these outputs are unchanged). `Cell.array_pbc` returns the wrapped positions instead of
  modifying its argument. The unused `stress` and `kstress` of the barostats (`stress` added the virials of
  each bead rather than their sum) and `coeffsc` of the integrators are removed.
- Depend objects that become computed after they are created (thermostats of the normal modes, constraints,
  contracted beads, bias scaling) are wired with `dsetfunc` or `dpipe` instead of setting `_func` by hand.
  `dsetfunc` raises an `AssertionError` if the object is already computed, synchronized, or a view. `dpipe`
  still replaces an earlier pipe (the time step of the normal modes follows the motion, then each dynamics
  bound to them), but now raises for synchronized objects and views, whose function it used to overwrite.
  Dependencies added after an object is created go through `ddeps`, which raises if the object is not
  computed; the `add_dependency`/`add_dependant` methods are removed.
- Depend objects are stored in attributes with a trailing underscore, `obj.x_` instead of `obj._x`, and
  `dproperties` exposes their value as `obj.x`. Entries of `dproperties` with no depend object behind them
  (they raised a `KeyError` when read) are removed.
- Depend objects follow one set of rules, stated in the docstring of `ipi/utils/depend.py`. The function of a
  computed object is held by the status shared with its views (`obj.x_._status.func`, it was `obj.x_._func`),
  and the conversions of a synchronized pair by the pair. Views are objects of their own classes
  (`depend_array_view`, `depend_immutable_view`). The stored array of a depend array is `_data` (it was
  `_wvalue`). All the writes (`set`, `obj.x = v`, item assignments, in-place operators, also through views) take
  the same steps, in the same order.
- Writing a computed depend object raises a `NameError` before changing anything (item assignments changed
  the value, and left it up to date). This also holds for its views, also those created before the object was
  made computed with `dsetfunc` or `dpipe` (e.g. the Atoms of contracted beads), which could be written, and
  the change was lost at the next recompute. `dsynchro` raises for computed objects and views.
- A write that keeps part of a depend array that is out of date (`obj.x_[i] = v`, `obj.x_ += v`, and any
  write to a view, e.g. `beads[0].p = v` after writing `nm.pnm`) brings the array up to date first: the new
  elements were merged into the old values, and writing a whole view lost the other elements. With `--jit`,
  the objects that a step updates in place or assigns in part (also through a view) are inputs of the recorded
  segments, so that values written between replayed steps (e.g. by `InteractiveSimulation.thermalize_momenta`,
  through `nm.pnm`) are used (H100: +1-3 us per replayed PIMD step). From Python 3.12, `obj.x_[a:b] op= v`
  assigned the updated view back to the array.
- Deep copies and pickles of depend objects are isolated (no dependants, no synchronization: a deep copy
  shared the dependants of the original, and an object with dependants could not be pickled), and views cannot
  be copied. Those of a synchronized object hold its current value, and can be written (they were computed from
  the peer of the original, and could not be written). Arrays are copied and pickled once (twice with numpy).
- Depend arrays behave as numpy arrays with all the backends. Indices with Ellipsis or None (`obj.x_[..., i]`,
  `obj.x_[None]`) select views connected to the array, as slices do (they were taken for scalar indices, and
  gave read-only copies); 0-d integer arrays select as integers do (they raised). With torch, `obj.x_[i, j] op= v`
  on a computed array modifies the stored tensor through the element, which shares memory with it, before the
  write is refused: the refused write now marks the array out of date, so that it is recomputed (the wrong
  value was kept). `reshape`/`flatten` raise
  if the view cannot be reshaped without a copy (writes through the copy were lost), which needs numpy 2.1.
  In-place operators raise a `TypeError` when they cannot update the stored array in place, e.g. a numpy array
  added to a torch tensor (which left the array unchanged), and with JAX when the result would change the shape
  of the array or its type beyond a `same_kind` cast (it replaced the array).
- New tool `i-pi-depend-graph` (`ipi.utils.tools.DependGraph`): writes the graph of the depend objects of a
  simulation in the dot format of Graphviz, optionally only around some objects (`--around`, `--depth`).
- Dead code removed: `Atoms.kin`/`kstress` (each replica of the beads had them as dependants of the momenta;
  `Beads.kins`, `kin` and `kstress` are computed from the momenta of all the beads), `Random.gfill`,
  `mintools.bracket_neb`/`min_brent_neb`/`L_BFGS_nls` (only used by their tests), the `k` argument of
  `L_BFGS`, `mathtools.sinch`/`mat_taylor`, `Forces.run`/`stop`, the `GradientMapper` and the unimplemented
  single-force-call step of the string method, and commented-out code. The SC stresses of the barostats
  (`kstress_sc`, `stress_sc`) are members of `BaroSCBZP` only. NEB and the string method printed the
  distances between the beads at every evaluation, and debugging information at each step: these are now
  messages at debug verbosity. The NEB, string and instanton modules no longer change the global print
  options of numpy when they are imported. The Al6xxx KMC motion prints its messages through i-PI, with the
  details of each event at debug verbosity and a summary of each step at high verbosity. Suzuki-Chin with
  an odd number of beads, or with a forward difference and ring-polymer contraction, raises an error
  instead of exiting.

The unit tests marked as strict xfail on main for the `b2tob1`, NMFD/eNMFD and `Fix.get_full_vector` (t=3, 4, 5)
bugs pass, and are no longer marked.

## Status and work left

### Pre-existing bugs found and not fixed (need a decision)
- Instanton: the spring term for non-uniform discretizations is not implemented (`SpringMapper` ignores the
  discretization; a disabled draft of it was removed), while the physical potential is weighted by the
  discretization coefficients.
- VSCF, not checked or changed: with `nprim > 1` the one-mode and the pair terms are scaled differently (the
  pair mean field is divided by `nprim` once more); `mptwo`, `threebody`, `print_mftpot`, `print_2b_map`,
  `print_vib_density`, `nkbt`, `nexc` and `ethresh` are accepted but not used (the MP2 correction and the
  three-body terms of Kapil et al. are not implemented); the solver uses a fixed basis for the one-mode
  problems (1.8e-4 Ha below the exact IMF free energy at kT ~ hbar w in the test model, which the IMF calculation
  converges), and has no limit on the number of SCF iterations.
- Constrained dynamics: see [Constrained dynamics](#constrained-dynamics) (angle and linear-group conditioning,
  tolerance units, batching).
- Socket interface (`ipi/interfaces/sockets.py`, `FFSocket`), also on main, left for the planned rewrite of the
  sockets with greenlets. They are fixed on `bugfix/clean-next-claude` (0a9a8f46, 5c1992ca, 4387eb32), whose
  regression tests of the shared-memory interface (72db7086: `SHM/plain`, `SHM/batched`, `SHM/noconsolidate`)
  can check the rewrite:
  - in worker mode, a free client that cannot take a request (e.g. one that answered HAVEDATA) makes the
    polling thread loop forever (the select mode stops when no client takes one), and in select mode the
    collection waits for replies also after the forcefield is stopped, so the interface cannot be closed
    while a client does not reply;
  - the drivers work on a duplicate of the accepted socket, whose descriptor names their shared-memory
    segments: it is freed at once and reused, so several shm clients get the same segments;
  - a dead client is closed only if EXIT can be sent, and the accepted socket, or the driver after a rejected
    handshake, is left open (in shm mode the driver finalizer does not close the socket);
  - a worker that finds its client in an unexpected state, or gets forces of the wrong size, leaves the
    request running and the client up, so the run hangs, and the exception is raised later in the polling
    thread; a client in an unexpected status is kept, with a warning at every poll;
  - removing a client waits for the workers of all the clients, and requeues requests that are already done;
  - a socket forcefield shared by systems with different numbers of atoms sends a larger system only the
    atoms of the first one (the others get zero forces), and fails in the polling thread for a smaller one:
    the indices of the active atoms, the socket batches and the shm segments assume a single size.
- Upstream, also on main: `test_metatomic.py` and `test_metatomic_reference.py` error at setup with Python 3.14,
  as TorchScript cannot compile metatomic's `AtomisticModel` with the deferred annotations of Python 3.14 (CI
  runs Python 3.10 and 3.13; the tests skip without metatomic).
- The metatomic LLPR test is a strict xfail for an upstream bug (metatrain).

### GPU performance and `torch.compile`
- The cell and barostat state stay on the device, so that whole steps can be compiled. The remaining host
  syncs in step code are graph breaks: the scalar acceptance logic of stochastic cell rescaling, the convergence
  loop of the constraint solver.
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
  `safe_stride` > 1. The state of the PRNG is kept as a copy and encoded as JSON only when the checkpoint is
  written or read back (`Random.snapshot`/`encode_state`): storing it takes ~50 us with all the backends, it
  took 120 us with numpy and 320 us with torch on the CPU (copying the state of the numpy MT19937 generator
  takes most of the 50 us). `InputArray.store` copies the arrays once.
- Reading the spring energy (`nm.vspring`, e.g. for the conserved quantity at every step) no longer computes
  the spring forces, which only the propagator of bosons and the `forces_spring` output read: 37 -> 20 us
  with 100 atoms and 8 beads, 890 -> 390 us with 1000 atoms and 32 beads, after a change of the positions.
  The virials of the beads and of the components are summed with `sum` instead of a loop; `FFDirect` applies
  the offset only if it is not zero (in place, it modified the energy returned by the driver);
  `ForceField.queue` copies the positions once.
- Host-side code that could move to the device if it becomes a bottleneck: properties estimators that loop over
  atoms and beads (on the arrays of the backend, with a host synchronization per element on GPUs), the bosonic
  exchange recursions (cycle energies and prefix/suffix potentials), the constraint solver (batch groups of the
  same shape).
- Benchmarks to repeat on a compute GPU (the numbers above come from a 2 GB Quadro P400 with slow float64).
- H100 with metatomic models through `ffdirect` (in-process steps without outputs, ms/step, i-PI overhead
  = step - model): with the torch backend the overhead is ~0.44 ms (NVE) and ~2.2 ms (4-bead PIMD), independent
  of size from 32 to 11000 atoms, i.e. kernel launches. Keeping the forcefield on the device lowered it by
  0.1-0.3 ms at small sizes and by 0.8 ms at 11000 atoms (LJ PIMD 20.3 -> 18.1 ms/step); PET-MAD-xs 4-bead
  PIMD of 96 water atoms went from 20.1 to 17.8 ms/step (numpy backend: 18.3). With the numpy backend and the
  model on the GPU, the overhead of 4-bead PIMD grows to 15 ms (1 BLAS thread) - 30 ms (16 threads) at 11000
  atoms, and was several times larger on a node shared with other jobs.

### CUDA graphs (`--jit` with torch on CUDA)
- `ipi/utils/cudagraphs.py` has the generic machinery: `GraphReplay` runs a function eagerly for the first
  calls (two, and at least one in the calling thread, as cuBLAS creates handles per thread), records the next
  one in segments separated by eager functions, and replays them. Code in the recorded function marks the
  points where the recording may be interrupted with `jit.split(obj)`, and calls functions eagerly
  with `jit.eager(func, *args)` (plain calls when not recording); these markers, the recording of the
  operations on the depend objects, `JitError` and the calls shared with jax.jit (`jit.Replay`: check of
  the fixed values, replay, eager calls, recording, final failure) are in `ipi/utils/jit.py`.
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
  thermostats are tensor scalars, so the temperature can change). The host values of all the segments are
  checked before the first one is replayed, so that a change does not leave half a step. After a replay,
  the taint status of all the objects downstream of those modified, the `synchro.manual` of synchronized
  pairs, and the objects held by tensor scalars are restored as at the end of the recorded segment.
- Each segment is captured in its own memory pool: with a pool shared by the segments, a tensor scalar written
  by a segment (e.g. the cell volume) could take the memory of the temporaries of an earlier segment, which
  overwrote it before a later segment of the next step read it. This happened when the positions were written
  before the recorded NPT step (the momentum of the barostat was off by 1e-4 after one replay); the graphs
  take ~2 MiB more per segment (P400: 2 -> 10 MiB for PIMD with the dummy driver, 5 segments).
- The recording is also interrupted at a force evaluation when the results of the forcefield are known, but
  not yet converted to device arrays (e.g. after the potential is read between two steps), as the conversion
  copies host data.
- A failed recording is final, as the kernels of the segment being captured have not run while the depend
  objects have moved on (with JAX, the objects written by the traced call hold tracers): the following calls
  raise `JitError`. Recording again crashed the CUDA caching allocator (`capture_begin` after an aborted
  capture into the same pool).
- Replayed: `Dynamics` with NVE/NVT/NPT/NST integrators (BAOAB and OBABO, MTS, PIMD) with no, Langevin, SVR,
  PILE_L or PILE_G thermostat, and the isotropic (BZP), flexible (MTK) or anisotropic (RGB) barostat with no
  or a Langevin thermostat; not replayed (eager, with a warning): stochastic cell rescaling, GLE,
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
- Replayed trajectories and outputs match eager ones (`test_cuda_graphs.py`), bit-identical without the
  functions compiled by torch.compile, to rounding errors with them (see kernel fusion below). The regression
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
  `jit.compiled` compiles a function with torch.compile when the steps are compiled with torch (`--jit`, also
  in the eager steps before the recording; plain call otherwise, or if torch.compile is unavailable, e.g.
  torch 2.9 on Python 3.14, or for GPUs older than compute capability 7.0); its arguments should be
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
- Possible extensions: stochastic cell rescaling (host logic and `float` syncs), evaluating the forces of all the components of a level at once
  (now each component is computed when its first bead is queued), cheaper staging of the inputs (host
  arrays such as `nmts` and `pdt` are compared with `np.array_equal` at each replay).

### Other backends
- cupy and jax: `xp` comes from `array_api_compat`, the PRNG has one implementation per backend, and all
  writes to the state go through the depend objects (`set`, `__setitem__`, in-place operators), which is where
  jax would need functional updates.
- JAX, eager: `set_backend("jax")` enables 64-bit types; depend arrays holding
  JAX arrays become `depend_immutable`, whose writes replace the stored array and whose views
  (`depend_immutable_view`: per-bead Atoms, thermostat `bind`s, `obj.x_[i] op= v`) read the parent at each
  access and write through it; masks and
  arrays of indices select copies, as with the other backends.
  `depend_array` writes go through two hooks (`_assign`, `_apply`), the only change for numpy/torch
  (~0.05 us per write, within noise in the benchmarks). Indexing views uses jitted take/put with the index as
  static argument (9-15 us), as eager JAX indexing costs 60-400 us. The PRNG splits a stored key per draw,
  in a single jitted call with the draw (40 us instead of 150 us). Only these array operations are jitted,
  not the algorithms. The MD path builds no arrays in place (see the behaviour changes above).
  Equivalence with numpy (same Gaussian draws, 1e-10): NVE, NVT Langevin/SVR/FFL/CL/multi, MTS with a bias,
  NPT isotropic/flexible (also with SVR on the barostat), NST, PIMD NVT/NPT with PILE_L/PILE_G, NVT-CC, fixed
  atoms, contracted components, Suzuki-Chin NVT/NPT (centred and forward differences), and the
  centroid-virial and gyration estimators; GLE xfails, also on the barostat (`stab_cholesky` and the
  auxiliary momenta `s` are modified in place). Motions that write into plain arrays (geometry optimization,
  NEB, string, phonons, IMF/VSCF, SCP, instantons, planetary, Al6xxx KMC, constrained dynamics) set
  `jax_support = False`, and `System.bind` stops with a `ValueError` (also within `multi`); dynamics,
  replay, ramps, alchemy and atom swaps run. Regression tests with `--ipi-backend jax`: 46 of 72 pass, with
  the MD restart tests (PLUMED is skipped); 21 stop at bind (geometry optimizers, NEB, phonons, instantons,
  constrained dynamics, and their restart tests), and the others at plain-array writes (bosonic exchange,
  open paths, `stab_cholesky` of the GLE thermostats).
  Eager JAX is slow: ms/step without forces, CPU, numpy / torch / jax, 100 atoms: NVE 0.049 / 0.12 / 0.78,
  NVT 0.11 / 0.29 / 1.8, PIMD-4 0.56 / 1.4 / 11.8, NPT 0.29 / 0.98 / 6.3, NPT flexible 0.93 / 2.6 / 13.4; at
  10000 atoms PIMD-4 28.9 / 26.5 / 34.8. Each eager JAX operation costs 11-16 us (arithmetic), 26 us (reshape), 40-85 us (indexing), 70 us
  (`zeros_like`), against <1 us with numpy; the depend helpers are <10% of the step. Performance needs
  jit of whole steps, which would need the depend state passed in and out of the traced function (as the
  CUDA-graph replay does).
- JAX, jit of the MD steps (`--jit` with JAX):
  `ipi.utils.jaxjit.JaxReplay` uses the recording of the CUDA graphs (`_Recorder`, the depend hooks, the
  restoration of the taint status, the `split`/`eager` markers, `jit.Replay`, `unsupported` and the warmup
  steps). The eager steps from the second are observed to find the depend objects read before being
  written, until two steps read the same ones (the second and the third, unless the state is written
  between them, e.g. positions written before a step make it recompute the forces; after 8 eager steps
  without agreement the steps stay eager, with a warning); the next one is traced, with these brought up
  to date (and the PRNG key) as arguments of `jax.jit`, returning the values of the objects it writes,
  compiled ahead of time (`lower().compile()`, so it never retraces), and replayed: inputs brought up to
  date with `get()`, host inputs (`dt`, `nmts`, `splitting`, ...) compared with their traced values
  (`JitError` if they change), outputs stored, statuses restored. depend changes: the `_value` of a
  `depend_immutable` notifies all the reads of its stored array (also those made by writes and views), and
  full assignments replace the stored array without reading it. Unlike CUDA graphs, a trace cannot be
  split: the driver calls marked with `eager` are traced, so only direct forcefields whose driver can be
  traced (`jax_traceable`, only `dummy`) are compiled; other forcefields (also metatomic, which converts JAX
  arrays to torch tensors) leave the steps eager, with a warning (they would
  need `jax.pure_callback`, i.e. a host round trip). Same trajectories as eager JAX to ~1e-15, same keys:
  NVE, NVT Langevin/SVR, PIMD with PILE_L/PILE_G, NPT isotropic and flexible, PIMD-NPT, MTS, NST. CPU,
  100 atoms, dummy forcefield, ms/step
  numpy / jax eager / jax jit: NVE 0.14 / 1.8 / 0.078, NVT 0.22 / 3.0 / 0.19, PIMD-4 0.93 / 16 / 0.78; at
  1000 atoms jit is slower than numpy with thermostats (NVT 0.43 vs 0.86, PIMD-4 2.6 vs 4.9), partly as
  XLA draws float64 Gaussian numbers 4x slower than numpy on the CPU (160 vs 40 us for 3000).
  H100 (kuma jobs 4679488, 4679541; Python 3.11, torch 2.12.1, JAX 0.10.2), ms/step without forces / with the
  dummy forcefield, 100 -> 10000 atoms: NVE graphs 0.031-0.036 / 0.078-0.084, jit 0.085-0.082 / 0.102-0.104;
  PIMD-8 (PILE_L) graphs 0.30-0.48 / 0.51-0.85, jit 0.30-0.35 / 0.43-0.48 (numpy 0.47-19.5 / 0.76-22.6, eager
  JAX 19-27). jit is 13-78x faster than eager JAX; CUDA graphs have the lower fixed cost (classical MD), jit
  scales better with beads and atoms (one XLA program, forces included). Jit PIMD is slower at 3162 atoms than
  at 10000 (not investigated). Jit vs eager JAX on the GPU: 2e-16 to 3e-15 after 10 steps, same keys.
  torch.compile and jax.jit do not overlap: with torch the steps are CUDA graphs and only the functions
  marked `@compiled` are compiled (whole-step torch.compile breaks on the depend locks and the lazy forces);
  with JAX the whole step is one XLA program and `@compiled` functions are traced inline.
- Normal modes on the H100 (job 4713726), ms/step without forces, 92a241d8 (FFT) -> matrix by default ->
  also the centre of mass removed from the centroid mode and the fused free ring-polymer step; NVT Langevin
  unchanged in all modes. PIMD-4 / PIMD-8 with PILE_L at 100 / 1000 / 10000 atoms: torch `--jit` 0.233 ->
  0.162 -> 0.129 / 0.258 -> 0.180 -> 0.151 / 0.308 -> 0.223 -> 0.194 and 0.289 -> 0.215 -> 0.187 / 0.337 ->
  0.245 -> 0.220 / 0.483 -> 0.346 -> 0.307; JAX `--jit` 0.235 -> 0.185 -> 0.181 / 0.255 -> 0.172 -> 0.176 /
  0.258 -> 0.202 -> 0.204 and 0.300 -> 0.253 -> 0.251 / 0.375 -> 0.264 -> 0.258 / 0.336 -> 0.290 -> 0.285;
  torch eager 1.46 -> 1.08 -> 1.04 and 2.06 -> 1.63 -> 1.59 (all sizes); numpy (1 thread) 0.312 -> 0.192 ->
  0.175 / 1.36 -> 0.81 -> 0.68 / 10.8 -> 6.8 -> 5.6 and 0.455 -> 0.322 -> 0.286 / 2.27 -> 1.59 -> 1.29 / 19.6 ->
  14.3 -> 11.4. 32-bead NPT at 100 / 1000 atoms: torch `--jit` 0.89 -> 0.78 -> 0.75 / 1.04 -> 0.91 -> 0.89,
  numpy 1.55 -> 1.30 -> 1.14 / 8.54 -> 7.75 -> 5.88.
- One option, `--jit` (`set_backend(..., jit=True)`, `InteractiveSimulation(jit=True)`), compiles the steps
  with either backend: CUDA graphs with torch on a CUDA device, jax.jit with JAX; with numpy, or torch on the
  CPU, it is an error (a usage error of `i-pi`). `InteractiveSimulation` raises if `device` or `jit` is given
  without `backend`, which it ignored.

### Multiple systems and replica exchange
- Systems are stepped in threads by default (`threading='true'`). With in-process forcefields this only adds
  overhead: classical NVT with the dummy forcefield, 100 atoms, numpy, 4 / 8 / 16 replicas: 0.76 / 1.54 / 3.0
  ms/step stepping the systems in turn, 4.0 / 8.6 / 17.7 ms/step in threads (workstation with the `powersave`
  governor; most of the cost is waking threads on other cores, the rest GIL contention). Threads cannot overlap
  in-process model evaluations either: `FFDirect.poll` evaluates under the forcefield lock (LJ metatomic model,
  8 replicas: no gain on CPU or CUDA). Batching the replicas (`threaded` `ffdirect` with `batch_size`) is the only
  in-process use of threads, and gained 0-11% with the LJ model.
- With socket drivers threads pay off when the driver cost is several ms: LJ model, 128 atoms, 4 / 8 replicas:
  19 / 39 ms/step in turn, 18 / 27 ms/step in threads; 2048 atoms: 140 / 285 vs 84 / 121 ms/step. With
  `matching='auto'` each replica keeps its driver, and stepping in turn leaves each driver idle most of the time
  (48 instead of 19 ms/step at 128 atoms, 4 replicas); `matching='any'` avoids it. The collect loop of the
  consolidated socket dispatch waits for all the jobs in flight before dispatching new requests.
- The exchange step tries all the pairs of replicas (with probability 1/stride each), each attempt swapping the
  ensembles, rescaling the momenta and recomputing `econs`/`lpens` through the depend graph: ~0.1-0.2 ms per pair,
  +3.4 ms/step for 8 replicas and +12 ms/step for 16 with stride 1.
- Options: step the systems in turn unless a forcefield runs outside the process or batches across systems;
  a single-threaded cooperative scheduler (e.g. greenlets) that interleaves the systems at the force evaluations,
  which would give batching and overlapping socket calls without threads; a cheaper acceptance test.

### Testing
- CI has no GPU: the torch job runs on the CPU, with implicit tensor->numpy conversions made to fail (as on
  GPUs). The sync-free CUDA test and CUDA runs must be done by hand on a GPU machine
  (`pytest --ipi-backend torch --ipi-device cuda`).
- The array-backend CI jobs run the unit tests and the regression tests with the torch backend, the regression
  tests with JAX (`tox -e backends`, with the `jax` extra), and the examples with torch (`tox -e torch-examples`);
  the unit tests use the backend of `--ipi-backend` for the simulations of the `simulation` and
  `reference_system` fixtures, the other tests select their backend (`select_backend`). With JAX, the
  regression tests of the features it does not support (`JAX_UNSUPPORTED` in `ipi_tests/test_tools.py`: they
  stop at the start, or write arrays in place) are strict expected failures.
- Opt-in tests are pytest options: `--ipi-mpi` (the ffmpi end-to-end tests, `tox -e mpitests`),
  `--ipi-metatomic-model PATH` (an exported PET-MAD-xs, otherwise exported with upet), `--ipi-fetch-driver`
  (fetching the Fortran driver from GitHub; the other `test_setup` tests compile it from the source tree).
- On CUDA, `test_restart.py::test_restart_matches_regression_reference[instanton]` fails about one time in
  three, with a 3e-10 difference in one element of the final Hessian (non-deterministic GPU arithmetic).
- PLUMED: the regression tests `PLUMED/METAD` and `PLUMED/BIAS` (they were `.NOAUTO`, with an outdated input
  tag and references that did not match main) are collected with the others, and skipped when PLUMED cannot be
  loaded; their references were regenerated (identical outputs on main and next). The `plumed.yml` workflow
  installs PLUMED from conda-forge and runs them with numpy and torch, together with the PLUMED unit tests.
- `test_softexit.py::test_signal_during_outputs_completes_the_step` failed now and then with "release unlocked
  lock": the test reset the soft exit while the monitoring thread was running its cleanup, and the reset
  replaced the lock it held. `softexit.reset` waits for a cleanup in progress. New unit test.

## Constrained dynamics

Scope: `ipi/utils/constrtools.py`, `ipi/engine/motion/constrained_dynamics.py`,
`ipi/inputs/motion/constrained_dynamics.py`.

### Fixed (no change of method)

- **Constrained dynamics could not run from an input file (also on main).**
  - `InputMotion` passes `fixatoms_dof=`, but `ConstrainedDynamics.__init__` accepted `fixatoms` and forwarded it
    to `Motion.__init__`, whose argument is now `fixatoms_dof`. It now takes `fixatoms_dof` (and still refuses
    fixed atoms).
  - `NVEConstrainedIntegrator.step_B` called the removed `Forces.forces_mts(level)`; it now uses
    `forces.mts_forces[level].f`, like the MD integrators.
  - There was no test: `ipi_tests/unit_tests/engine/test_constrained_dynamics.py` now runs bond, angle and
    Eckart constraints from XML, checks that the constraints hold, and checks numpy vs torch.
- **Small errors with obvious intent.**
  - An invalid ensemble now raises (the `ValueError` was created but not raised).
  - The "no convergence" warning of `ConstraintSolver.proj_manifold` was never printed
    (`if i == self.maxit` after `range(self.maxit)`); it now uses a `for ... else`.
  - `array.shape = ...` on the constructor arguments, which mutated the caller's arrays and is deprecated since
    numpy 2.5, is replaced by `reshape`.
  - `constraint_list=[]` (mutable default) is replaced by `None`.
  - `norm_order` -1 (the infinity norm) was stored as `inf`, so the restart could not be read; it is stored
    as -1.
- **Dead code and docs.** Removed `ConstraintList.get_iai` (it called a method that no constraint defines) and
  the `dq`/`dqprev` dproperties of `ConstraintList` (no such depend objects). Fixed the docstrings of `Gfunc`
  and `GCfunc`, and the comment on `_m3`. Also removed a debugging block in `RigidBondConstraint.gfunc` that
  printed "autsch" and called `exit()` if `q[0]` was infinite. Removed the all-zero `constraint_values` of
  `EckartConstraint`, on which its constraint and Jacobian depended without using them; scipy is imported once,
  in constrtools.
- **Index check, defaults, dependencies, coverage.**
  - `EckartConstraint` rejected any atom indices that were not sorted, as "repeated atom indices" (the reference
    positions are matched to the sorted atoms): the atoms are sorted, with their reference positions, and only
    repeated ones are rejected.
  - The solver defaults are those of the input (tolerance 1e-4, 1000 iterations) also in `ConstraintSolver` and
    in the fallback of `ConstrainedDynamics` (they were 1e-3/1000 and 1e-4/10000).
  - The Jacobians of bonds and angles no longer depend on the target values, so changing a target does not
    recompute them, the Gram matrix and its factorization.
  - Regression test with reference outputs: `CONSTRAINTS/rigid_water`, 64 rigid q-TIP4P/F molecules (two bonds
    and an angle each) in NVE. New unit tests for the Eckart indices and the bond Jacobian.
- **Array-backend port.** The per-constraint loops of `RigidBondConstraint`/`AngleConstraint` are vectorized
  (same values); scipy's Cholesky solve is kept for numpy, and torch solves the Gram systems with
  `xp.linalg.solve`. Trajectories match the pre-port code to ~5e-16 and numpy/torch CPU/torch CUDA to ~1e-15.

### To change (method or behaviour)

1. **Angle constraints are ill-conditioned near 0 and 180 degrees.** `AngleConstraint` constrains
   `cos(angle) - cos(target)`, whose gradient vanishes at 0 and 180 degrees: the Gram matrix becomes singular
   for (nearly) linear angles.
2. **Eckart constraints assume a non-linear group.** Six constraints are always imposed; for linear groups (or two
   atoms) the rotational rows of the Jacobian are linearly dependent and the Cholesky factorization fails.
3. **Inconsistent convergence measure.** `proj_manifold` tests `norm(g)` over a group, where `g` mixes units:
   squared lengths minus squared targets (bonds), cosine differences (angles), mass-weighted displacements
   (Eckart). The tolerance means different things for different constraints, and for bonds it scales with the
   bond length; a per-constraint relative tolerance would be more meaningful.
4. **Performance for many molecules, and on GPU.** The solver loops over constraint groups in Python, with a
   small dense solve and a host-side convergence check (`float(norm(g))`) per group and iteration. For many
   small molecules (rigid water: one group per molecule) this is slow with numpy and dominated by kernel
   launches and device syncs with torch. Groups of the same type and size could be batched
   (`xp.linalg.solve`/`cholesky` accept stacks of matrices), with one convergence check per iteration. On numpy
   the scipy `cho_factor` is cached and reused; on torch the Gram system is re-solved every time (a cached
   `torch.linalg.cholesky` + `torch.cholesky_solve` would match the numpy behaviour).
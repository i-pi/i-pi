# Constrained dynamics: review notes

Scope: `ipi/utils/constrtools.py`, `ipi/engine/motion/constrained_dynamics.py`,
`ipi/inputs/motion/constrained_dynamics.py`.

## Fixed in the array-backend PR (trivial, no change of method)

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
- **Dead code and docs.** Removed `ConstraintList.get_iai` (it called a method that no constraint defines) and
  the `dq`/`dqprev` dproperties of `ConstraintList` (no such depend objects). Fixed the docstrings of `Gfunc`
  and `GCfunc`, and the comment on `_m3`. Also removed a debugging block in `RigidBondConstraint.gfunc` that
  printed "autsch" and called `exit()` if `q[0]` was infinite.
- **Array-backend port.** The per-constraint loops of `RigidBondConstraint`/`AngleConstraint` are vectorized
  (same values); scipy's Cholesky solve is kept for numpy, and torch solves the Gram systems with
  `xp.linalg.solve`. Trajectories match the pre-port code to ~5e-16 and numpy/torch CPU/torch CUDA to ~1e-15.

## To change in a follow-up PR (method or behaviour)

1. **Misleading Eckart check.** `EckartConstraint` raises "repeated atom indices" whenever
   `constrained_indices != np.unique(...)`, i.e. also for unsorted but unique indices. Either sort the indices
   (keeping the reference coordinates consistent) or test for duplicates only.
2. **Angle constraints are ill-conditioned near 0 and 180 degrees.** `AngleConstraint` constrains
   `cos(angle) - cos(target)`, whose gradient vanishes at 0 and 180 degrees: the Gram matrix becomes singular
   for (nearly) linear angles.
3. **Eckart constraints assume a non-linear group.** Six constraints are always imposed; for linear groups (or two
   atoms) the rotational rows of the Jacobian are linearly dependent and the Cholesky factorization fails.
4. **Inconsistent convergence measure.** `proj_manifold` tests `norm(g)` over a group, where `g` mixes units:
   squared lengths minus squared targets (bonds), cosine differences (angles), mass-weighted displacements
   (Eckart). The tolerance means different things for different constraints, and for bonds it scales with the
   bond length; a per-constraint relative tolerance would be more meaningful.
5. **Inconsistent defaults.** `ConstrainedDynamics` builds a default `ConstraintSolver(tolerance=1e-4,
   maxit=10000)`, while `ConstraintSolver` defaults to `1e-3` and `1000`, and the input class has its own
   defaults.
6. **Unneeded recomputation.** `ValueConstraintBase.bind` makes `Dg` depend on `constraint_values`, but the
   Jacobians of bonds and angles do not depend on the target values, so `Dg`, `Gram` and its factorization are
   recomputed when the values change.
7. **Performance for many molecules, and on GPU.** The solver loops over constraint groups in Python, with a
   small dense solve and a host-side convergence check (`float(norm(g))`) per group and iteration. For many
   small molecules (rigid water: one group per molecule) this is slow with numpy and dominated by kernel
   launches and device syncs with torch. Groups of the same type and size could be batched
   (`xp.linalg.solve`/`cholesky` accept stacks of matrices), with one convergence check per iteration. On numpy
   the scipy `cho_factor` is cached and reused; on torch the Gram system is re-solved every time (a cached
   `torch.linalg.cholesky` + `torch.cholesky_solve` would match the numpy behaviour).
8. **Coverage.** The only example (`examples/temp/broken/lammps/h2o-constraints`) is in the "broken" folder,
   and there is no regression test with reference outputs.

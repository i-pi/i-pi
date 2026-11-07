Combined friction driver
========================

A friction driver can return both the bath coupling factors and the bath-induced
mean-field energy and forces through the ordinary i-PI socket interface. Model
evaluation, coupling construction and differentiation remain in the external
driver. One cached request supplies all these quantities.

Input
-----

Declare a dedicated socket forcefield::

   <ffsocket name="ace" mode="unix" pbc="false">
     <address>hcu_ace</address>
   </ffsocket>

Select it inside the existing ``friction`` element, alongside the bath parameters::

   <forcefield>ace</forcefield>

Do not also include this forcefield in the system's potential-energy ``forces``
list: the friction machinery adds its correction energy and forces once.
The physical potential (for example MACE) retains its usual force component.
Without this option, the existing internal mean-field calculation is retained.
The combined interface requires variable, non-Markovian friction, a dedicated
forcefield for each system, all coordinates active in the socket request,
zero energy offset, and no forcefield batching.

Ordinary socket requests
------------------------

The driver receives the whole ring in one POSDATA request: ``P*N`` positions,
ordered by bead, then atom, then Cartesian component. The cell follows the
usual i-PI wire convention (row-major matrix, lattice vectors in columns).
Periodic wrapping must be disabled to preserve a common image for the beads.

The reserved INIT parameter ``ring_polymer`` contains a base64-encoded UTF-8
JSON object. Encoding protects its commas from the existing parameter/checkpoint
parser. The object contains:

* ``schema``: ``ipi_friction_ring_v1``; ``units``: ``atomic``.
* ``nbeads`` and ``natoms``: the actual ring dimensions.
* ``alpha``: mean-field coefficients in i-PI normal-mode order.
* ``transform``: the orthonormal bead-to-normal-mode matrix ``C``.

The normal GETFORCE reply contains the total ring correction energy and its
negative Cartesian derivative, in atomic units, without division by P.
For example, for

.. math::

   U_{\mathrm{mf}} = \frac12\sum_{n,a}\alpha_n [(C F)_{na}]^2,

return :math:`-\partial U_{\mathrm{mf}}/\partial q_{bi}`, including moving
references, environment atoms and cross-bead dependencies.

The JSON extras contain ``friction_ring_schema="ipi_friction_ring_v1"`` and
the existing ``friction_coupling_jacobian``, ``friction_gamma`` and
``friction_meta`` entries. Factors and tensors have a leading bead dimension;
metadata identifies active atoms, channel labels and units. The driver may
also return ``friction_coupling``. The factor used by the stochastic bath and
the derivative used for the correction are separately specified quantities;
their relationship must be validated for the chosen model.

The H/Cu ACEfriction driver's ``--mode ring`` evaluates native Sigma, builds
the centroid-relative endpoint coupling and differentiates the resulting
correction completely. It reuses the ACE evaluations for both outputs.
Its existing bead and centroid modes remain available. The factor-two tensor
correction remains in the ACE driver (a square-root-of-two factor on Sigma).

Restart and integration
-----------------------

The external calculation must be stateless with respect to the trajectory.
Model caches are permitted; hidden trajectory or random state is not.
Checkpoints retain the forcefield selection, bath auxiliary variables and
random-number state. Restart initializes a fresh driver with the mode data.
Keep the driver, model and package versions available for reproducibility.
The current correction virial is not propagated; use fixed-cell dynamics.

The ``symmetric-gle`` splitting requires ``nve-f``, a non-Markovian bath,
closed-ring RPMD propagation, one force timestep and ``fixcom=False``.
Fixed substrate coordinates must have collapsed ring polymers, zero initial
momenta and no bath coupling. Their reaction forces are excluded from momentum
updates.


Array backends
--------------

Electronic friction uses the runtime-selected ``xp`` namespace for propagation
and supports NumPy, Torch and JAX arrays. Select a backend using
``i-pi --backend numpy input.xml`` (or ``torch`` / ``jax``). Socket payloads,
input validation and SciPy spectral integration and OU matrix preparation use
host NumPy arrays. OU matrices are cached for each timestep/temperature pair
and transferred to the selected backend. Checkpoints explicitly convert the
auxiliary state to host arrays and preserve the selected backend's RNG state.
Use the same backend when continuing a stochastic trajectory exactly.

Dependency objects follow the ``next`` conventions: ``obj.q`` reads a value,
while ``obj.q_`` supports tracked partial writes. Auxiliary-state updates use
array replacement, which also works for immutable JAX arrays. Model-specific
coupling differentiation remains in the external driver.

The tests in ``test_friction_backends.py`` compare matched-noise trajectories,
native-RNG restarts and conservative timestep convergence on all three CPU
backends. They include scalar and variable friction, the combined ring driver,
fixed substrate coordinates and short harmonic/double-well trajectories.
These are numerical regression checks, not converged correlation-function
benchmarks. CUDA execution has not been validated by these CPU tests.
Friction integrators currently use eager execution; the existing ``next`` JIT
selection excludes them and falls back to eager steps.

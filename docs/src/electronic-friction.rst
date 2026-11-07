Electronic friction
===================

The electronic-friction integrators combine physical forces with a Markovian
or auxiliary-variable bath. Configure the bath in the ``friction`` element
of ``motion/dynamics``. Physical potential evaluation remains in the usual
forcefield driver; model-specific coupling construction belongs in an external
driver.

Bath inputs
-----------

``bath_mode=markovian`` supports one bead. Its momentum update solves the
Ornstein--Uhlenbeck process at fixed positions. Spatially varying friction
requires a positive-semidefinite friction tensor and the associated metadata.
A scalar, constant coupling is selected with ``variable_friction=False`` and
``sigma_static``.

``bath_mode=non-markovian`` requires externally fitted ``Ap`` matrices for
every normal mode. For P beads, their shape is ``(P, 1+naux, 1+naux)``;
the centroid comes first, followed by the native i-PI normal-mode ordering.
A two-dimensional matrix is accepted only for one bead. Each block has the
form ``[[0, theta.T], [-theta, A]]``. Stability and the auxiliary noise
covariance are checked when the bath is prepared. The fitting procedure must
use the intended temperature, bead number, units and spectral normalization.
i-PI does not generate or fit the mode spectra.

``Lambda`` supplies rows ``[omega, Lambda(omega)]`` in atomic units for the
mean-field coefficients. These coefficients are integrated from the spectrum,
independently of the fitted matrices. The grid must start at zero; shape-preserving
interpolation is used without extrapolation. Check convergence with frequency spacing
and cutoff, as well as the accuracy of each fitted bath. A normalized spectrum
and a coupling factor must together reproduce the intended physical friction;
rescaling either changes that convention. ``debug_alpha_input`` overrides the
integration and is intended for explicit validation cases.

Propagation and restarts
------------------------

For non-Markovian, closed-ring RPMD, ``splitting=symmetric-gle`` uses symmetric
half updates of the momentum-to-auxiliary coupling around free-ring propagation.
It requires ``mode=nve-f``, one force timestep and ``fixcom=False``. Assess
timestep error using mode-resolved observables as well as centroid observables.

Fixed coordinates must have collapsed beads, zero momenta and zero bath
coupling. The combined external correction currently supports fixed-cell
simulations; it does not propagate a correction virial.

XML checkpoints save the auxiliary variables and random-number state, as well
as the nuclear state and bath inputs. Continue with the same backend, driver
and models to reproduce a stochastic trajectory. Tests cover restart both
before the first step and after propagation, including restored bath memory.

Coupling models
---------------

An exact scalar coupling F and its Jacobian dF/dQ define both the mean-field
force and the bath coupling. A learned factor of a friction tensor need not be
the gradient of any scalar coupling. When a driver constructs an approximate
F from that factor, validate the resulting derivative and the bath factor
together. Differentiating the correction completely makes its force
conservative; it does not by itself make a separately supplied bath factor
equal to dF/dQ.

The :doc:`external-mean-field` interface lets a driver return the correction
energy, its complete derivative and the bath factors in one ring request.
This keeps atomistic coupling approximations outside the integration engine.

Validation and array backends
-----------------------------

NumPy, Torch and JAX CPU propagation is exercised by the friction unit tests.
Array algebra uses the selected backend; socket serialization, spectral
quadrature and SciPy matrix preparation use host arrays. Friction propagation
currently runs eagerly, and these CPU tests do not establish CUDA support.

The tests include independent splitting updates, matched-noise trajectories,
checkpoint continuation, fixed coordinates and harmonic/double-well models.
Long correlation functions require separate statistical benchmarks. An exact
linear one-step map and its discrete stationary covariance also provide a
sampling-independent check for a harmonic potential and constant coupling.

Migration from development inputs
---------------------------------

The older ``trenins-rossi`` splitting and ``centroid_endpoint_trapezoid``
coupling option have been removed. Old inputs/checkpoints require explicit
migration; changing a splitting does not preserve an identical trajectory.
Use ``symmetric-gle`` for the supported symmetric RPMD-GLE propagation.
Supply exact coupling values and their derivatives, or construct and
differentiate an approximate correction in the external driver. The old
centroid-specific input keys and contracted-extras exception are also removed.
Frozen production runtimes can still read their original checkpoints.

The spectrum utility now writes Lambda by default, includes zero frequency,
and extends to twenty exponential decay frequencies. Select ``--quantity J``
explicitly for J and ``--omega-max`` for another grid cutoff. Spectrum grid
and cutoff convergence remain the user's responsibility.

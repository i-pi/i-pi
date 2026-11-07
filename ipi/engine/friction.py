"""Electronic-friction operators for friction-enabled dynamics in i-PI.

This module separates electronic friction into two layers.

`Friction`
    High-level operator owned by the `nve-f` / `nvt-f` integrators. It reads a
    canonical, channel-resolved coupling Jacobian and friction matrix from
    force-driver extras and evaluates the optional friction mean-field (MF)
    contribution (`friction_coupling_nm`, `force_mf`, `energy_mf`).

`FrictionGLE`
    Bath propagator used by `Friction.step(pdt)`. In the current implementation
    it covers the Markovian electronic-friction bath for both:
      - static friction (`sigma_static`)
      - variable friction (canonical driver Jacobian and Gamma)
    and Non-Markovian for both static and variable friction.

Current Markovian/Non-Markovian behavior:
    1. `Friction.step(pdt)` applies the MF momentum kick, if enabled.
    2. `FrictionGLE.step(pdt)` applies the dissipative/stochastic bath update.

For variable friction, drivers provide arrays in a model-independent wire
format. When only friction-active coordinates are returned, `active_dofs` or
`active_atoms` metadata is used to embed reduced arrays into the full Cartesian
system. Representation-specific packing belongs in the driver.
"""

# debug_mf_mode is a physical CL mean-field/counterterm toggle: on/off.


import base64
import json
import numpy as np
from scipy.linalg import expm, cholesky

from ipi.engine.motion import Motion
from ipi.engine.normalmodes import NormalModes
from ipi.engine.beads import Beads
from ipi.engine.atoms import Atoms
from ipi.engine.forces import ForceBead

from ipi.utils import nmtransform
from ipi.utils.depend import depend_value, dproperties, ddeps, dcompute
from ipi.utils.array_backend import xp, xp_size, to_numpy
from ipi.utils.messages import info, verbosity, warning


# Markovian static
def _apply_mass_scaled_ou(p, sm, drift, noise_scale, noise):
    """Applies an exact OU step in mass-scaled coordinates."""

    p_ms = p / sm
    deltah = xp.sum(p_ms * p_ms)
    p_ms = drift * p_ms + noise_scale * noise
    deltah -= xp.sum(p_ms * p_ms)
    return p_ms * sm, 0.5 * deltah


# Markovian variable
def _apply_mass_scaled_matrix_ou(s, A, dt, kbt, noise):
    """Applies an exact matrix OU step in mass-scaled coordinates.
    We diagonalise to independent modes similar to R. J. Maurer et al, PRL, 2017"""

    et = 0.5 * float(xp.sum(s * s))

    A = 0.5 * (
        A + A.T
    )  # Confusing notation. The A matrix here is not the same A matrix defined in GLE case.
    evals, evecs = xp.linalg.eigh(
        A
    )  # This is mode diagonalization to treat each mode seperately as independent OU process.

    evals = xp.clip(evals, 0.0, None)
    c = xp.exp(-evals * dt)
    s2 = xp.sqrt(1.0 - c * c)
    y = evecs.T @ s
    y = c * y + kbt**0.5 * s2 * noise
    s = evecs @ y

    et -= 0.5 * float(xp.sum(s * s))
    return s, et


# Non-markovian variable
def _compute_aux_ou_matrices(A, dt, kbt):
    """Build exact OU drift/noise matrices for auxiliary variables. i.e builds T and S.

    The auxiliary covariance is canonical with variance kBT in each auxiliary
    coordinate, so S S^T = kBT * (I - T T^T).
    """

    A = np.asarray(A, dtype=float)
    T = expm(-A * dt)
    C = float(kbt) * (np.eye(A.shape[0]) - T @ T.T)
    C = 0.5 * (C + C.T)
    try:
        S = cholesky(C, lower=True, check_finite=False)
    except Exception:
        evals, evecs = np.linalg.eigh(C)
        evals = np.clip(evals, 0.0, None)
        S = evecs @ (np.sqrt(evals)[:, None] * evecs.T)
    return T, S


class FrictionBath:
    """Base class for electronic-friction bath propagators."""

    def __init__(self):
        self.friction = None
        self.motion = None

    def bind(self, friction, motion: Motion | None) -> None:
        self.friction = friction
        self.motion = motion

    @property
    def ediss(self) -> float:
        return float(self.friction.ediss)

    def state_shape(self):
        return None

    def initialize_state(self) -> None:
        pass

    def step(self, pdt: float) -> None:
        raise NotImplementedError


class FrictionGLE(FrictionBath):
    """Unified electronic-friction bath propagator.

    Markovian friction uses an exact frozen-position momentum OU update.
    Non-Markovian auxiliary variables are kept separate from the physical
    system momentum: unlike `ThermoGLE`, `self.s` is reserved for auxiliary
    bath variables only.
    """

    def __init__(self):
        super().__init__()
        self.s = None
        self.ns = 0
        self.theta = None
        self.A_aux = None
        self.T_aux = None
        self.S_aux = None
        self._ou_key = None
        # Valid only while one non-Markovian GLE operator is being applied.
        # Coordinates and driver extras are unchanged within that operation.
        self._step_nm_jacobian = None
        self._step_nm_active_dofs = None
        self._cache_nm_jacobian_within_step = False

    def bind(self, friction, motion: Motion | None) -> None:
        friction._validate_bath_bead_count()
        super().bind(friction, motion)
        self.initialize_state()

    def initialize_state(self) -> None:
        if self.is_markovian():
            self.s = xp.zeros((0,), dtype=xp.float64)
            self.ns = 0
            self.theta = None
            self.A_aux = None
            self.T_aux = None
            self.S_aux = None
            return

        # Non-Markovian state:
        # `s` stores auxiliary bath variables only, not the physical momentum.
        # The user-facing input follows George's notation and supplies Ap,
        # whose first row gives theta and whose lower-right block is A.
        # Trenins--Rossi RPMD-GLE uses a different auxiliary embedding for
        # each free-ring-polymer normal mode.  ``Ap_modes[m]`` has the form
        # [[0, theta_m^T], [-theta_m, A_m]].  A two-dimensional Ap is accepted
        # only for the classical/single-bead limit by ``_validate_Ap``.
        Ap_modes = self.friction._mode_resolved_Ap()
        self.theta = xp.asarray(Ap_modes[:, 0, 1:], dtype=xp.float64)
        self.A_aux = np.asarray(Ap_modes[:, 1:, 1:], dtype=float).copy()
        self.ns = int(self.theta.shape[1])
        self.T_aux = None
        self.S_aux = None
        self._ou_key = None
        self._nm_transform_matrix = None
        self._step_nm_jacobian = None
        self._step_nm_active_dofs = None
        self._cache_nm_jacobian_within_step = False
        if self.friction.variable_friction:
            self.s = None
            return

        self._initialize_auxiliary_state(self.state_shape())

    def _initialize_auxiliary_state(self, shape) -> None:
        """Restore bath memory without consuming random numbers on restart."""
        saved = self.friction._initial_s
        if saved.size:
            if saved.shape != shape:
                raise ValueError(
                    f"Friction restart s has shape {saved.shape}; expected {shape}."
                )
            self.s = xp.asarray(saved, copy=True)
            self.friction._initial_s = np.zeros(0)
        else:
            self.s = xp.zeros(shape, dtype=xp.float64)
            if self.friction.prng is not None and self.ns > 0:
                self.s = self.friction._kbt_rp() ** 0.5 * self.friction.prng.gvec(shape)

    def is_markovian(self) -> bool:
        return str(self.friction.bath_mode) == "markovian"

    def state_shape(self):
        if self.is_markovian():
            return (0,)
        if self.friction.variable_friction:
            return None if self.s is None else self.s.shape
        return (
            int(self.friction.beads.nbeads),
            int(self.ns),
            3 * int(self.friction.beads.natoms),
        )

    def _get_non_markovian_jacobian(self) -> np.ndarray:
        coupling_jacobian = xp.asarray(
            self.friction.coupling_jacobian, dtype=xp.float64
        )
        active_dofs = self.friction._active_coupling_dofs(coupling_jacobian.shape[-1])
        coupling_jacobian = coupling_jacobian[..., active_dofs]
        self._step_nm_active_dofs = active_dofs
        if coupling_jacobian.ndim != 3:
            raise ValueError(
                "Variable non-markovian friction requires coupling_jacobian "
                "with shape (nbeads, nchannels, ndof)."
            )
        expected_shape = (
            int(self.friction.beads.nbeads),
            int(self.ns),
            int(coupling_jacobian.shape[1]),
        )
        if self.s is None:
            self._initialize_auxiliary_state(expected_shape)
        if self.s.shape != expected_shape:
            raise ValueError(
                "Non-markovian auxiliary state shape is inconsistent with the "
                f"current Jacobian payload. Expected {expected_shape}, got {self.s.shape}."
            )
        return coupling_jacobian

    def _get_nm_transform_matrix(self) -> np.ndarray:
        nbeads = int(self.friction.beads.nbeads)
        if self._nm_transform_matrix is None or self._nm_transform_matrix.shape != (
            nbeads,
            nbeads,
        ):
            self._nm_transform_matrix = xp.asarray(
                nmtransform.mk_nm_matrix(nbeads), dtype=xp.float64
            )
        return self._nm_transform_matrix

    def _get_non_markovian_nm_jacobian(self) -> np.ndarray:
        """Return the active-DOF dF_nm[n']/dQ_nm[n] transform.

        The cached value is reused only by substeps of one ``step`` call.
        """
        if self._cache_nm_jacobian_within_step and self._step_nm_jacobian is not None:
            return self._step_nm_jacobian
        coupling_jacobian = self._get_non_markovian_jacobian()
        cmat = self._get_nm_transform_matrix()
        transformed = xp.einsum("rb,nb,bci->rnci", cmat, cmat, coupling_jacobian)
        if self._cache_nm_jacobian_within_step:
            self._step_nm_jacobian = transformed
        return transformed

    def step(self, pdt: float) -> None:
        if self.is_markovian():
            self._step_markovian(pdt)
            return
        self._step_nm_jacobian = None
        self._step_nm_active_dofs = None
        self._cache_nm_jacobian_within_step = True
        try:
            self._step_non_markovian(pdt)
        finally:
            # Do not retain a coordinate-dependent object across outer
            # splitting operations, where positions may have changed.
            self._step_nm_jacobian = None
            self._step_nm_active_dofs = None
            self._cache_nm_jacobian_within_step = False

    def _step_markovian(self, pdt: float) -> None:
        if pdt <= 0.0:
            return
        if self.friction.variable_friction:
            self._step_markovian_variable(pdt)
        else:
            self._step_markovian_static(pdt)

    def _step_markovian_static(self, pdt: float) -> None:
        friction = self.friction
        sigma = float(friction.sigma_static)
        gamma = sigma * sigma
        if gamma < 0.0:
            raise ValueError(
                "gamma must be non-negative for Markovian linear coupling."
            )
        if gamma == 0.0:
            return

        m = friction.nm.dynm3
        p = friction.nm.pnm
        sm = xp.sqrt(m)
        gamma_nm = xp.full_like(m, gamma)
        drift = xp.exp(-(gamma_nm / m) * pdt)
        noise_scale = xp.sqrt(friction._kbt_rp() * (1.0 - drift * drift))
        p_new, bath_energy = _apply_mass_scaled_ou(
            p, sm, drift, noise_scale, friction.prng.gvec(p.shape)
        )
        friction.nm.pnm = p_new
        friction.beads.p = friction.nm.transform.nm2b(p_new)
        friction.ediss += bath_energy

    def _step_markovian_variable(self, pdt: float) -> None:
        friction = self.friction
        coupling_jacobian = friction._get_coupling_jacobian()
        nbeads = coupling_jacobian.shape[0]
        gamma = xp.asarray(friction.gamma, dtype=xp.float64)
        p = friction.beads.p
        updated = []
        m = friction.beads.m3
        sm = xp.sqrt(m)
        et = 0.0

        for b in range(nbeads):
            gam = gamma[b, :, :]
            inv_sm = 1.0 / sm[b, :]
            A = (inv_sm[:, None] * gam) * inv_sm[None, :]
            s = p[b, :] * inv_sm
            s, et_b = _apply_mass_scaled_matrix_ou(
                s,
                A,
                pdt,
                friction._kbt_rp(),
                friction.prng.gvec(s.shape),
            )
            et += et_b
            updated.append(s * sm[b, :])

        friction.beads.p = xp.stack(updated)
        friction.ediss += et

    def _step_non_markovian(self, pdt: float) -> None:
        # Symmetric non-Markovian bath update at fixed coordinates:
        #   O_s(dt/2)   : Eq. (S38)
        #   B_P,F(dt/2) : Eq. (S40)
        #   B_s(dt)     : Eq. (S42)
        #   B_P,F(dt/2) : Eq. (S40)
        #   O_s(dt/2)   : Eq. (S38)
        if pdt <= 0.0:
            return

        self.os_step(0.5 * pdt)
        self.bp_f_step(0.5 * pdt)
        self.bs_step(pdt)
        self.bp_f_step(0.5 * pdt)
        self.os_step(0.5 * pdt)

    def os_step(self, pdt: float) -> None:
        """Auxiliary OU step, Eq. (S38).
        Updates auxiliary momenta.

        Implements:
            s <- T_{dt} s + S_{dt} xi

        This acts on auxiliary bath variables only. Physical system momenta
        remain in `beads.p` / `nm.pnm`.
        """

        if pdt <= 0.0 or self.ns == 0:
            return
        if self.s is None:
            self._get_non_markovian_jacobian()
        key = (pdt, self.friction._kbt_rp())
        if getattr(self, "_ou_key", None) != key:
            matrices = [_compute_aux_ou_matrices(a, *key) for a in self.A_aux]
            self.T_aux = xp.asarray(np.stack([t for t, _ in matrices]))
            self.S_aux = xp.asarray(np.stack([s for _, s in matrices]))
            self._ou_key = key
        et = 0.5 * float(xp.sum(self.s * self.s))
        # s has shape (nmodes, naux, nchannel).  Each normal mode has its own
        # OU drift/noise matrices, fitted to K^(n)(t).
        self.s = xp.einsum("mab,mbd->mad", self.T_aux, self.s)
        noise = self.friction.prng.gvec(self.s.shape)
        self.s = self.s + xp.einsum("mab,mbd->mad", self.S_aux, noise)
        et -= 0.5 * float(xp.sum(self.s * self.s))
        self.friction.ediss += et

    def bp_f_step(self, pdt: float) -> None:
        """Bath (aux) -to- physical momentum coupling, Eq. (S40).

        Implements the momentum kick generated by the current auxiliary bath
        state through:
            P <- P - dt * dF/dQ * theta^T s

        """

        if pdt <= 0.0 or self.ns == 0:
            return
        p = self.friction.nm.pnm
        if self.friction.variable_friction:
            nm_jacobian = self._get_non_markovian_nm_jacobian()
            theta_s = xp.einsum("ra,rac->rc", self.theta, self.s)
            active_dofs = self._step_nm_active_dofs
            if active_dofs is None:
                raise RuntimeError(
                    "Missing active DOFs for variable non-Markovian friction."
                )
            self.friction.nm.pnm_[:, active_dofs] -= pdt * xp.einsum(
                "rc,rnci->ni", theta_s, nm_jacobian
            )
            self.friction.beads.p = self.friction.nm.transform.nm2b(
                self.friction.nm.pnm
            )
            return
        else:
            m = self.friction.nm.dynm3
            sm = xp.sqrt(m)
            sigma = float(self.friction.sigma_static)
            theta_s = xp.einsum("ma,mad->md", self.theta, self.s)
            p_ms = p / sm
            p_ms -= pdt * (sigma / sm) * theta_s
            p_new = p_ms * sm
        self.friction.nm.pnm = p_new
        self.friction.beads.p = self.friction.nm.transform.nm2b(p_new)

    def bs_step(self, pdt: float) -> None:
        """Physical momentum-to-bath (aux) coupling, Eq. (S42).

        The coupling Jacobian contracts with velocity (P / M), then theta
        maps the resulting channel velocity onto the auxiliary variables.
        """

        if pdt <= 0.0 or self.ns == 0:
            return
        p = self.friction.nm.pnm
        m = self.friction.nm.dynm3
        if self.friction.variable_friction:
            nm_jacobian = self._get_non_markovian_nm_jacobian()
            active_dofs = self._step_nm_active_dofs
            if active_dofs is None:
                raise RuntimeError(
                    "Missing active DOFs for variable non-Markovian friction."
                )
            drive = xp.einsum("rnci,ni->rc", nm_jacobian, (p / m)[:, active_dofs])
            self.s = self.s + pdt * self.theta[:, :, None] * drive[:, None, :]
        else:
            sigma = float(self.friction.sigma_static)
            drive = sigma * p / m
            self.s = self.s + pdt * self.theta[:, :, None] * drive[:, None, :]


class Friction:
    """
    Friction operator for friction-enabled dynamics.
    """

    # -------------------------
    # Input-configured
    # -------------------------
    variable_friction: bool
    bath_mode: str  # "none" | "markovian" | "non-markovian"
    debug_mf_mode: str  # "on" | "off"

    Lambda: np.ndarray  # [omega, J(omega)] for non-markovian OU fit
    Ap: np.ndarray  # mode-resolved George-style momentum + auxiliary drift matrices
    debug_alpha_input: np.ndarray  # optional [omega_k, alpha]

    sigma_static: float

    # Canonical force-driver extras
    coupling_jacobian_key: str
    gamma_key: str
    friction_meta_key: str

    # -------------------------
    # Runtime bound
    # -------------------------
    beads: Beads
    """Reference to the beads"""
    nm: NormalModes
    """Reference to the normal modes"""

    def __init__(
        self,
        variable_friction: bool = True,
        bath_mode: str = "non-markovian",  # can be 1. none (no dissipative, no random force),
        debug_mf_mode: str = "on",
        # on  - apply Caldeira-Leggett mean-field/counterterm.
        # off - no friction contribution to conservative force.
        Lambda=np.zeros((0, 2), float),
        Ap=np.zeros((0, 0), float),
        # Trenins--Rossi mode-resolved auxiliary drift matrices.  For P beads,
        # the faithful RPMD-GLE representation has shape
        #   (P, 1+naux, 1+naux)
        # with Ap[n] = [[0, theta_n^T], [-theta_n, A_n]].  The 2-D
        # shape is retained only for P=1 (classical GLE).
        debug_alpha_input=np.zeros((0, 2), float),
        sigma_static: float = 1.0,
        # Constant factor for independent Cartesian channels.
        coupling_jacobian_key: str = "friction_coupling_jacobian",
        gamma_key: str = "friction_gamma",
        friction_meta_key: str = "friction_meta",
        coupling_key: str = "friction_coupling",
        s=None,
        ediss: float = 0.0,
        bath_contract: str = "",
        forcefield: str = "",
    ):
        """Initialises the friction object.
        Args:
            Lambda: Friction spectrum J(omega)/omega in atomic units, supplied
                as columns [omega, Lambda(omega)] on a grid starting at zero.
                Defaults to np.zeros(0, float).
            debug_alpha_input: Normal-mode coefficients in expression for the frictional mean-field
                potential [Eq. (8b) in https://doi.org/10.1103/PhysRevLett.134.226201].
                Defaults to np.zeros(0, float).
            variable_friction (bool, optional): True if the gradient of the friction coupling F(q)
                [introduced in Eq. (5) of https://doi.org/10.1103/PhysRevLett.134.226201]
                depends on position.
                Defaults to False.
        """

        # Choices
        self.variable_friction = bool(variable_friction)
        self.bath_mode = str(bath_mode)
        self.debug_mf_mode = self._normalize_debug_mf_mode(debug_mf_mode)

        # Optional whole-ring forcefield: ordinary i-PI requests, model work external.
        self.forcefield = str(forcefield)
        self._ring_force = None
        self.driver_data_ = depend_value(
            name="friction_driver_data", func=self._get_driver_data
        )

        # Kernel shape
        self.Lambda = np.asanyarray(Lambda, dtype=float).copy()
        self.Ap = np.asanyarray(Ap, dtype=float).copy()
        self.debug_alpha_input = np.asanyarray(debug_alpha_input, dtype=float).copy()

        self.coupling_jacobian_ = depend_value(
            name="coupling_jacobian", func=self._get_coupling_jacobian
        )
        self.gamma_ = depend_value(
            name="gamma", func=self._get_gamma, dependencies=[self.coupling_jacobian_]
        )

        #     # Friction coupling: F(q), such that Σ{i,α} = ∂F(q) / ∂q{i,α}
        self.friction_coupling_nm_ = depend_value(
            name="friction_coupling_nm",
            func=self.get_friction_coupling_nm,
            dependencies=[self.coupling_jacobian_],
        )
        # Frictional mean-field force
        self.force_mf_nm_ = depend_value(
            name="force_mf_nm",
            func=self.get_force_mf_nm,
            dependencies=[
                self.driver_data_ if self.forcefield else self.friction_coupling_nm_
            ],
        )

        # force_meanfield
        self.force_mf_ = depend_value(
            name="force_mf", func=self.get_force_mf, dependencies=[self.force_mf_nm_]
        )

        # Conserved mean-field potential
        self.energy_mf_ = depend_value(name="energy_mf", value=0.0)

        self.sigma_static = float(sigma_static)

        self.coupling_jacobian_key = str(coupling_jacobian_key)
        self.gamma_key = str(gamma_key)
        self.friction_meta_key = str(friction_meta_key)
        self.coupling_key = str(coupling_key)
        self._friction_meta = {}
        self._friction_atoms_idx: np.ndarray | None = None
        self._friction_dof_idx: np.ndarray | None = None
        self._channel_labels: list[str] = []
        self._friction_contract_signature = None
        if bath_contract:
            dofs, labels, jacobian_units, gamma_units = json.loads(bath_contract)
            self._friction_contract_signature = (
                tuple(dofs),
                tuple(labels),
                jacobian_units,
                gamma_units,
            )
        self._initial_s = np.array(s if s is not None else [], dtype=float, copy=True)
        if self._initial_s.size:
            if self.bath_mode != "non-markovian":
                raise ValueError("Friction restart s requires a non-markovian bath.")
            if self._initial_s.ndim != 3 or not np.all(np.isfinite(self._initial_s)):
                raise ValueError(
                    "Friction restart s must be a finite three-dimensional array."
                )
        self._gamma_reconstruction_error = 0.0
        self._nm_transform_matrix = None
        self.bath: FrictionBath | None = None

        # runtime handles
        self.alpha: np.ndarray | None = None
        self.forces = None
        self.ensemble = None
        self.prng = None

        # bookkeeping: cumulative energy exchange with bath (Markovian)
        self.ediss_ = depend_value(
            name="ediss", value=float(ediss)
        )  # positive = system -> bath

    # ==========================================================================
    # bind
    # ==========================================================================

    def bind(self, motion: Motion) -> None:
        self.beads = motion.beads
        self.nm = motion.nm
        self.ensemble = motion.ensemble
        self.forces = motion.ensemble.forces
        self.prng = motion.prng
        self._validate_bath_bead_count()
        if self.bath_mode not in ("none", "markovian", "non-markovian"):
            raise ValueError(
                "bath_mode must be one of: 'none', 'markovian', 'non-markovian'."
            )

        self.debug_mf_mode = self._normalize_debug_mf_mode(self.debug_mf_mode)
        if self.debug_mf_mode not in ("on", "off"):
            raise ValueError("debug_mf_mode must be one of: 'on', 'off'.")

        # Non-Markovian OU fit requires spectral density
        if self.bath_mode == "non-markovian":
            if self.Ap.size == 0:
                raise ValueError(
                    "non-markovian requires externally fitted Ap matrices. "
                    "i-PI does not fit Ap from Lambda."
                )
            self._validate_Ap()

        # Setup alpha for MF (may be zeros if MF disabled)
        self.alpha = xp.asarray(self._setup_alpha())
        if self.forcefield:
            self._bind_ring_forcefield()

        info(
            "Friction.bind:\n"
            f"  variable_friction = {self.variable_friction}\n"
            f"  bath_mode           = {self.bath_mode}\n"
            f"  debug_mf_mode       = {self.debug_mf_mode}\n"
            f"  sigma_static     = {self.sigma_static}\n"
            f"  Ap shape           = {self.Ap.shape}\n"
            f"  coupling_jacobian_key = '{self.coupling_jacobian_key}'\n"
            f"  gamma_key             = '{self.gamma_key}'\n",
            verbosity.low,
        )

        self.bath = self._build_bath()
        if self.bath is not None:
            self.bath.bind(self, motion)

        # Dependencies
        ddeps(
            self.coupling_jacobian_,
            [self.driver_data_ if self.forcefield else self.forces.extras_],
        )
        ddeps(self.friction_coupling_nm_, [self.beads.q_])
        dcompute(
            self.energy_mf_,
            self.get_energy_mf,
            [self.driver_data_ if self.forcefield else self.friction_coupling_nm_],
        )

    def _validate_bath_bead_count(self) -> None:
        if self.bath_mode == "markovian" and self.beads.nbeads > 1:
            raise ValueError(
                "Markovian electronic friction requires nbeads=1. "
                "For nbeads>1, use bath_mode='non-markovian' with externally "
                "fitted Ap matrices for each ring-polymer normal mode; "
                "a Markovian physical bath can produce internal-mode memory."
            )

    def _build_bath(self) -> FrictionBath | None:
        if self.bath_mode == "none":
            return None
        if self.bath_mode in ("markovian", "non-markovian"):
            return FrictionGLE()
        raise RuntimeError("bath_mode must be one of none, markovian or non-markovian")

    def _mode_resolved_Ap(self) -> np.ndarray:
        """Return ``Ap`` with an explicit leading normal-mode dimension.

        The leading index follows i-PI's native normal-mode ordering, with the
        centroid at index zero. Fits for degenerate modes must be repeated in
        the corresponding positions of this array by the external fitter.
        """

        Ap = np.asarray(self.Ap, dtype=float)
        if Ap.ndim == 2:
            return Ap[np.newaxis, ...]
        if Ap.ndim == 3:
            return Ap
        raise ValueError(
            "Ap must have shape (1+naux, 1+naux) for P=1 or "
            "(P, 1+naux, 1+naux) for mode-resolved RPMD-GLE."
        )

    def _validate_Ap(self) -> None:
        Ap_input = np.asarray(self.Ap, dtype=float)
        nbeads = int(self.beads.nbeads)
        if Ap_input.ndim == 2 and nbeads != 1:
            raise ValueError(
                "A single 2-D Ap describes only the centroid/classical kernel and "
                "cannot be broadcast over a multi-bead ring polymer. Supply "
                "mode-resolved Ap with shape (P, 1+naux, 1+naux), fitted to "
                "K^(n)(t) for the current temperature and bead count."
            )

        Ap_modes = self._mode_resolved_Ap()
        if Ap_modes.shape[0] != nbeads:
            raise ValueError(
                f"Mode-resolved Ap has {Ap_modes.shape[0]} modes, but the ring "
                f"polymer has {nbeads} beads."
            )
        if Ap_modes.shape[1] != Ap_modes.shape[2] or Ap_modes.shape[1] < 2:
            raise ValueError(
                "Each mode-resolved Ap must be square with shape " "(1+naux, 1+naux)."
            )
        if not np.all(np.isfinite(Ap_modes)):
            raise ValueError("Ap contains non-finite values.")

        for mode, Ap in enumerate(Ap_modes):
            A = Ap[1:, 1:]
            if A.size == 0:
                raise ValueError("Ap must contain at least one auxiliary variable.")
            if not np.isclose(Ap[0, 0], 0.0, rtol=0.0, atol=1e-14):
                raise ValueError(
                    f"Ap[{mode},0,0] must be zero for the CL/GLE convention."
                )
            if not np.allclose(Ap[1:, 0], -Ap[0, 1:], rtol=1e-10, atol=1e-14):
                raise ValueError(
                    "Each mode-resolved Ap must follow the CL/GLE convention "
                    "Ap[n,1:,0] == -Ap[n,0,1:]."
                )
            if np.any(np.real(np.linalg.eigvals(A)) <= 0.0):
                raise ValueError(
                    f"Auxiliary drift matrix A for normal mode {mode} is not stable."
                )
            symmetric_A = 0.5 * (A + A.T)
            if np.min(np.linalg.eigvalsh(symmetric_A)) < -1e-12:
                raise ValueError(
                    "The symmetric part of auxiliary drift matrix A for normal "
                    f"mode {mode} is not positive semidefinite. It cannot define "
                    "the requested canonical OU noise covariance."
                )

    @staticmethod
    def _normalize_debug_mf_mode(mode: str) -> str:
        mode = str(mode).strip().lower()
        if mode == "none":
            warning(
                "debug_mf_mode='none' is deprecated; use 'off'. Treating it as 'off'.",
                verbosity.low,
            )
            return "off"
        if mode not in ("on", "off"):
            raise ValueError("debug_mf_mode must be one of: 'on', 'off'.")
        return mode

    def _ensure_bath_bound(self) -> None:
        if self.bath is not None or self.bath_mode == "none":
            return
        self.bath = self._build_bath()
        if self.bath is not None:
            self.bath.bind(self, None)

    # ==========================================================================
    # temperature helper
    # ==========================================================================

    def _kbt_rp(self) -> float:
        """kB * (P*T) for ring-polymer effective classical temperature."""
        from ipi.utils.units import Constants

        kb = Constants.kb
        return kb * float(self.ensemble.temp) * float(self.beads.nbeads)

    # ==========================================================================
    # Alpha setup (MF)
    # ==========================================================================

    def _setup_alpha(self) -> np.ndarray:
        wk = np.asarray(to_numpy(self.nm.omegak), dtype=float)
        nmodes = wk.size

        # If MF disabled, keep alpha zeros.
        if self.debug_mf_mode == "off":
            return np.zeros(nmodes, dtype=float)

        # Explicit debug coefficients override the physical bath spectrum.
        if self.debug_alpha_input.size > 0:
            if self.debug_alpha_input.ndim != 2 or self.debug_alpha_input.shape[1] != 2:
                raise ValueError(
                    "debug_alpha_input must have shape (nmodes,2) [omega_k, alpha]."
                )
            if self.debug_alpha_input.shape[0] != nmodes:
                raise ValueError(
                    f"debug_alpha_input rows ({self.debug_alpha_input.shape[0]}) != nmodes ({nmodes})."
                )
            if not np.allclose(self.debug_alpha_input[:, 0], wk):
                raise ValueError(
                    "debug_alpha_input omega_k does not match current nm.omegak."
                )
            info("Friction: using alpha from debug_alpha_input table.", verbosity.low)
            alpha = np.asarray(self.debug_alpha_input[:, 1], dtype=float)
            if not np.all(np.isfinite(alpha)) or np.any(alpha < 0.0):
                raise ValueError(
                    "debug_alpha_input coefficients must be finite and nonnegative."
                )
            if np.any(alpha[wk == 0.0] != 0.0):
                raise ValueError("debug_alpha_input centroid coefficient must be zero.")
            return alpha

        # MF comes from the physical spectrum (paper Eq. 8b), not the GLE fit.
        if self.Lambda.size > 0:
            if self.Lambda.ndim != 2 or self.Lambda.shape[1] != 2:
                raise ValueError(
                    "Lambda must have shape (nfrequency,2) [omega, Lambda(omega)]."
                )
            if not np.all(np.isfinite(self.Lambda)) or np.any(self.Lambda[:, 1] < 0.0):
                raise ValueError(
                    "Lambda must contain finite frequencies and nonnegative spectral values."
                )
            omega = np.asarray(self.Lambda[:, 0], dtype=float)
            if omega.size < 2:
                raise ValueError("Lambda must contain at least two points.")
            if np.any(omega < 0.0) or np.any(np.diff(omega) <= 0.0):
                raise ValueError(
                    "Lambda omega must be nonnegative and strictly increasing."
                )
            alpha = get_alpha_numeric(Lambda=self.Lambda[:, 1], omega=omega, omegak=wk)
            info("Friction: computed alpha^(n) from Lambda.", verbosity.low)
            return alpha

        raise ValueError(
            "Friction: MF requires Lambda [omega, Lambda(omega)] or debug_alpha_input. Ap matrices specify dynamics only."
        )

    # ==========================================================================
    # Canonical driver payload
    # ==========================================================================

    @staticmethod
    def _decode_json(value):
        if isinstance(value, str):
            try:
                return json.loads(value)
            except json.JSONDecodeError:
                return value
        return value

    def _get_consistent_meta(self, key: str, bead_resolved: bool) -> dict:
        """Read metadata and require identical dictionaries across beads."""
        meta = self._driver_extras().get(key)
        if meta is None:
            raise KeyError(f"Missing required canonical friction metadata '{key}'.")
        meta = self._decode_json(meta)
        if isinstance(meta, (list, tuple)):
            if len(meta) == 0:
                raise ValueError(f"Canonical friction metadata '{key}' is empty.")
            decoded = [self._decode_json(item) for item in meta]
            if bead_resolved and len(decoded) != int(self.beads.nbeads):
                raise ValueError(
                    f"'{key}' must contain one entry per bead; got {len(decoded)} "
                    f"for {self.beads.nbeads} beads."
                )
            first = decoded[0]
            if any(item != first for item in decoded[1:]):
                raise ValueError(
                    f"Canonical friction metadata '{key}' differs across beads."
                )
            meta = first
        if not isinstance(meta, dict):
            raise ValueError(
                f"Canonical friction metadata '{key}' must be a dictionary."
            )
        if meta.get("schema") != "ipi_friction_v1":
            raise ValueError(
                f"'{key}.schema' must be 'ipi_friction_v1', got {meta.get('schema')!r}."
            )
        return meta

    def _configure_active_dofs(self, meta: dict) -> None:
        natoms = int(self.beads.natoms)
        ndof = 3 * natoms
        if "active_dofs" in meta:
            dofs = np.asarray(meta["active_dofs"], dtype=int).reshape(-1)
            atoms = np.unique(dofs // 3)
        elif "active_atoms" in meta:
            atoms = np.asarray(meta["active_atoms"], dtype=int).reshape(-1)
            dofs = (
                np.concatenate(
                    [np.arange(3 * atom, 3 * atom + 3, dtype=int) for atom in atoms]
                )
                if atoms.size
                else np.zeros(0, dtype=int)
            )
        else:
            raise ValueError(
                f"'{self.friction_meta_key}' must define zero-based active_dofs or active_atoms."
            )
        if dofs.size == 0 or np.any(dofs < 0) or np.any(dofs >= ndof):
            raise ValueError(
                f"Active Cartesian DOFs must be unique indices in [0, {ndof - 1}]."
            )
        if np.unique(dofs).size != dofs.size:
            raise ValueError(
                "Canonical friction metadata contains duplicate active DOFs."
            )
        if np.any(atoms < 0) or np.any(atoms >= natoms):
            raise ValueError(
                f"Active atoms must be zero-based indices in [0, {natoms - 1}]."
            )
        self._friction_dof_idx = dofs
        self._friction_atoms_idx = atoms

    def _embed_jacobian(self, reduced: np.ndarray) -> np.ndarray:
        ndof = 3 * int(self.beads.natoms)
        if reduced.shape[-1] == ndof:
            return reduced
        if self._friction_dof_idx is None or reduced.shape[-1] != len(
            self._friction_dof_idx
        ):
            raise ValueError(
                f"Reduced coupling Jacobian shape {reduced.shape} is inconsistent with active DOFs."
            )
        full = np.zeros(reduced.shape[:-1] + (ndof,), dtype=float)
        full[..., self._friction_dof_idx] = reduced
        return full

    def _active_coupling_dofs(self, full_ndof: int) -> np.ndarray:
        """Return Cartesian columns that can carry variable friction.

        The canonical driver Jacobian is embedded in the full system with
        exactly zero columns for inactive atoms.  Expensive normal-mode
        contractions can therefore operate on this reduced set and scatter
        their result back into full Cartesian space afterwards.
        """
        if self._friction_dof_idx is None:
            return np.arange(full_ndof, dtype=int)
        active = np.asarray(self._friction_dof_idx, dtype=int)
        if active.ndim != 1 or np.any(active < 0) or np.any(active >= full_ndof):
            raise ValueError("Invalid friction-active DOF indices.")
        return active

    def _embed_gamma(self, reduced: np.ndarray) -> np.ndarray:
        ndof = 3 * int(self.beads.natoms)
        if reduced.shape[-2:] == (ndof, ndof):
            return reduced
        nactive = 0 if self._friction_dof_idx is None else len(self._friction_dof_idx)
        if reduced.shape[-2:] != (nactive, nactive):
            raise ValueError(
                f"Reduced Gamma shape {reduced.shape} is inconsistent with {nactive} active DOFs."
            )
        full = np.zeros(reduced.shape[:-2] + (ndof, ndof), dtype=float)
        for bead in range(reduced.shape[0]):
            full[bead][np.ix_(self._friction_dof_idx, self._friction_dof_idx)] = (
                reduced[bead]
            )
        return full

    def _get_coupling_jacobian(self) -> np.ndarray:
        """Return dF_c(Q_b)/dQ_i as (nbeads, nchannels, full_ndof)."""
        if not self.variable_friction:
            return float(self.sigma_static)
        meta = self._get_consistent_meta(self.friction_meta_key, bead_resolved=True)
        self._configure_active_dofs(meta)
        labels = meta.get("channel_labels")
        if not isinstance(labels, list) or not all(
            isinstance(label, str) for label in labels
        ):
            raise ValueError(
                f"'{self.friction_meta_key}.channel_labels' must be a list of strings."
            )
        signature = (
            tuple(int(dof) for dof in self._friction_dof_idx),
            tuple(labels),
            meta.get("jacobian_units"),
            meta.get("gamma_units"),
        )
        if self._friction_contract_signature is None:
            self._friction_contract_signature = signature
        elif signature != self._friction_contract_signature:
            raise ValueError(
                "Canonical friction active DOFs, channel labels, or units changed "
                "during the simulation."
            )
        payload = self._decode_json(
            self._driver_extras().get(self.coupling_jacobian_key)
        )
        if payload is None:
            raise KeyError(
                f"Missing canonical driver extra '{self.coupling_jacobian_key}'."
            )
        jacobian = np.asarray(payload, dtype=float)
        expected_prefix = (int(self.beads.nbeads), len(labels))
        if jacobian.ndim != 3 or jacobian.shape[:2] != expected_prefix:
            raise ValueError(
                f"'{self.coupling_jacobian_key}' must have shape "
                f"(nbeads, nchannels, nactive_dof); got {jacobian.shape}."
            )
        if not np.all(np.isfinite(jacobian)):
            raise ValueError(
                f"'{self.coupling_jacobian_key}' contains non-finite values."
            )
        self._friction_meta = meta
        self._channel_labels = list(labels)
        return xp.asarray(self._embed_jacobian(jacobian))

    def _get_gamma(self):
        """Return driver-supplied Gamma after validating its canonical factor
        The validation is done by constructing Gamma from J^T*J."""
        if not self.variable_friction:
            return float(self.sigma_static) ** 2
        jacobian = to_numpy(self.coupling_jacobian)
        payload = self._decode_json(self._driver_extras().get(self.gamma_key))
        if payload is None:
            raise KeyError(f"Missing canonical driver extra '{self.gamma_key}'.")
        gamma = np.asarray(payload, dtype=float)
        if gamma.ndim != 3 or gamma.shape[0] != int(self.beads.nbeads):
            raise ValueError(
                f"'{self.gamma_key}' must have shape (nbeads, nactive_dof, nactive_dof); "
                f"got {gamma.shape}."
            )
        if not np.all(np.isfinite(gamma)):
            raise ValueError(f"'{self.gamma_key}' contains non-finite values.")
        gamma = self._embed_gamma(gamma)
        reconstructed = np.einsum("bci,bcj->bij", jacobian, jacobian)
        diff = np.linalg.norm(gamma - reconstructed, axis=(1, 2))
        scale = np.maximum(np.linalg.norm(gamma, axis=(1, 2)), 1.0e-30)
        self._gamma_reconstruction_error = float(np.max(diff / scale))
        if not np.allclose(gamma, reconstructed, rtol=1.0e-8, atol=1.0e-12):
            raise ValueError(
                f"Driver Gamma is inconsistent with coupling Jacobian J^T J; "
                f"maximum relative Frobenius error={self._gamma_reconstruction_error:.3e}."
            )
        return xp.asarray(gamma)

    # ==========================================================================
    # Friction forces and coupling
    # ==========================================================================

    def _get_nm_transform_matrix(self) -> np.ndarray:
        nbeads = int(self.beads.nbeads)
        if self._nm_transform_matrix is None or self._nm_transform_matrix.shape != (
            nbeads,
            nbeads,
        ):
            self._nm_transform_matrix = xp.asarray(
                nmtransform.mk_nm_matrix(nbeads), dtype=xp.float64
            )
        return self._nm_transform_matrix

    def _get_friction_coupling(self) -> np.ndarray:
        coupling = self._driver_extras().get(self.coupling_key)
        if coupling is None:
            raise KeyError(
                f"Did not find '{self.coupling_key}' among the force extras = {self.forces.extras}"
            )
        if isinstance(coupling, str):
            try:
                coupling = json.loads(coupling)
            except json.JSONDecodeError:
                pass
        coupling = np.asarray(coupling, dtype=float)
        if coupling.ndim == 1:
            coupling = coupling[:, np.newaxis]
        if coupling.ndim != 2:
            raise ValueError(
                f"{self.coupling_key} must have shape (nbeads, nbath). Got {coupling.shape}."
            )
        if coupling.shape[0] != int(self.beads.nbeads):
            raise ValueError(
                f"{self.coupling_key} shape {coupling.shape} incompatible with nbeads={self.beads.nbeads}."
            )
        if not np.all(np.isfinite(coupling)):
            raise ValueError(f"{self.coupling_key} contains non-finite values.")
        return xp.asarray(coupling)

    def _get_nm_coupling_jacobian(self) -> np.ndarray:
        coupling_jacobian = xp.asarray(self.coupling_jacobian, dtype=xp.float64)
        active_dofs = self._active_coupling_dofs(coupling_jacobian.shape[-1])
        coupling_jacobian = coupling_jacobian[..., active_dofs]
        cmat = self._get_nm_transform_matrix()
        return xp.einsum("rb,nb,bci->rnci", cmat, cmat, coupling_jacobian)

    def get_friction_coupling_nm(self):
        """Compute the friction coupling for each normal-mode index"""
        if self.variable_friction:
            # Refresh the upstream dependency before evaluating F. Otherwise a
            # later Jacobian read can invalidate F during force evaluation,
            # leaving a cached MF force clean while its dependency is stale.
            _ = self.coupling_jacobian
            # In the variable friction case, we must get the coupling F(Q) from the driver. A future example of this could be
            # an ML model that diretly provides the coupling for atomistic systems.
            return self._get_nm_transform_matrix() @ self._get_friction_coupling()
        else:
            # Match the scalar bath's independent Cartesian channels:
            # F_d(q) = sigma_static * q_d. Summing coordinates here would
            # introduce a collective bath and spurious cross-coordinate forces.
            return float(self.sigma_static) * self.nm.qnm

    def _bind_ring_forcefield(self):
        """Use the existing force request machinery with bead-major coordinates."""
        if not self.variable_friction or self.bath_mode != "non-markovian":
            raise ValueError(
                "A ring friction forcefield requires variable non-Markovian friction."
            )
        if self.forcefield not in self.forces.ff:
            raise ValueError(f"Unknown friction forcefield {self.forcefield!r}.")
        ff = self.forces.ff[self.forcefield]
        if any(component.ff is ff for component in self.forces.mforces):
            raise ValueError(
                "The ring friction forcefield must not also be a potential-energy force component."
            )
        if getattr(ff, "_ring_friction_owner", None) is not None:
            raise ValueError(
                "A ring friction forcefield cannot be shared between systems."
            )
        if ff.dopbc or list(ff.active) != [-1] or ff.offset != 0:
            raise ValueError(
                "Ring friction requires pbc=False, all coordinates active, and zero energy offset."
            )
        if getattr(ff, "batch_size", 1) != 1:
            raise ValueError(
                "A whole-ring request cannot be combined with forcefield batching."
            )
        ff._ring_friction_owner = self
        config = dict(
            schema="ipi_friction_ring_v1",
            units="atomic",
            nbeads=int(self.beads.nbeads),
            natoms=int(self.beads.natoms),
            alpha=to_numpy(self.alpha).tolist(),
            transform=to_numpy(self._get_nm_transform_matrix()).tolist(),
        )
        # Ordinary socket INIT data; no new socket message types.
        ff.pars = dict(ff.pars or {})
        ff.pars["ring_polymer"] = base64.b64encode(
            json.dumps(config, allow_nan=False).encode("utf-8")
        ).decode("ascii")
        self._ring_atoms = Atoms(self.beads.natoms * self.beads.nbeads)
        dcompute(self._ring_atoms.q_, lambda: self.beads.q.reshape(-1), [self.beads.q_])
        self._ring_force = ForceBead()
        self._ring_force.bind(
            self._ring_atoms, self.forces.cell, ff, self.forces.output_maker
        )
        ddeps(self.driver_data_, [self._ring_force.ufvx_])

    def _get_driver_data(self):
        energy, force, _, extras = self._ring_force.ufvx
        energy = float(energy)
        force = np.asarray(force, dtype=float)
        expected = xp_size(self.beads.q)
        if (
            force.size != expected
            or not np.isfinite(energy)
            or not np.all(np.isfinite(force))
        ):
            raise ValueError("Ring friction driver returned invalid energy or forces.")
        if (
            not isinstance(extras, dict)
            or extras.get("friction_ring_schema") != "ipi_friction_ring_v1"
        ):
            raise ValueError(
                "Ring friction driver requires friction_ring_schema='ipi_friction_ring_v1'."
            )
        return energy, force.reshape(self.beads.q.shape), extras

    def _driver_extras(self):
        return self.driver_data_.get()[2] if self.forcefield else self.forces.extras

    def get_energy_mf(self):
        """Compute the frictional potential of mean field, Eq. (S19) of https://doi.org/10.1103/PhysRevLett.134.226201"""

        if self.debug_mf_mode == "off":
            return 0.0

        if self.forcefield:
            return self.driver_data_.get()[0]

        coupling_nm = xp.asarray(self.friction_coupling_nm, dtype=xp.float64)
        if coupling_nm.ndim == 1:
            weighted_coupling2 = self.alpha * coupling_nm**2
        else:
            weighted_coupling2 = self.alpha[:, xp.newaxis] * coupling_nm**2
        return xp.sum(weighted_coupling2) / 2

    def get_force_mf_nm(self):
        """Negative derivative of the frictional potential of mean field with respect to normal modes"""
        if self.forcefield:
            return self._get_nm_transform_matrix() @ xp.asarray(
                self.driver_data_.get()[1]
            )
        if self.variable_friction:
            reduced_force = -xp.einsum(
                "r,rc,rnci->ni",
                self.alpha,
                xp.asarray(self.friction_coupling_nm, dtype=xp.float64),
                self._get_nm_coupling_jacobian(),
            )
            ndof = self.nm.pnm.shape[-1]
            active_dofs = self._active_coupling_dofs(ndof)
            # Compact scatter map: O(nactive * ndof), without an ndof-square identity.
            scatter = xp.asarray(
                xp.asarray(active_dofs)[:, None] == xp.arange(ndof)[None, :],
                dtype=xp.float64,
            )
            return reduced_force @ scatter

        return (
            -self.alpha[:, xp.newaxis]
            * self.friction_coupling_nm
            * float(self.sigma_static)
        )

    def get_force_mf(self):
        """Negative derivative of the frictional potential of mean field with respect to bead positions"""
        return self.nm.transform.nm2b(self.force_mf_nm)

    # ==========================================================================
    # main step
    # ==========================================================================

    def mf_step(self, pdt: float) -> None:
        """Apply only the conservative mean-force momentum kick."""
        if self.debug_mf_mode == "on":
            self.beads.p = self.beads.p + self.force_mf * pdt

    def step(self, pdt: float) -> None:
        """
        Apply friction operator over time interval pdt:
          - MF kick (if enabled)
          - bath kick (markovian/non-markovian) - dissipative and random forces
        """

        # MF
        self._validate_bath_bead_count()
        self.mf_step(pdt)

        self._ensure_bath_bound()  # Just makes sure a bath is attached, if not , tries to initialize a bath
        if self.bath is None:
            return
        self.bath.step(pdt)


dproperties(
    Friction,
    [
        "coupling_jacobian",
        "gamma",
        "friction_coupling_nm",
        "energy_mf",
        "ediss",
        "force_mf_nm",
        "force_mf",
    ],
)


def get_alpha_numeric(
    Lambda: np.ndarray, omega: np.ndarray, omegak: np.ndarray
) -> np.ndarray:
    """Integrate a nonnegative spectrum on an explicit zero-to-cutoff grid.

    Shape-preserving interpolation of the integrand avoids negative lobes.
    No extrapolation is performed; callers must converge the grid and cutoff.
    """
    from scipy.interpolate import PchipInterpolator

    omega = np.asarray(omega, dtype=float)
    spectrum = np.asarray(Lambda, dtype=float)
    wk = np.asarray(omegak, dtype=float)
    if (
        omega.ndim != 1
        or spectrum.shape != omega.shape
        or omega.size < 2
        or not np.all(np.isfinite(omega))
        or not np.all(np.isfinite(spectrum))
        or np.any(spectrum < 0.0)
        or np.any(np.diff(omega) <= 0.0)
        or omega[0] != 0.0
    ):
        raise ValueError(
            "Lambda requires finite nonnegative values on a strictly increasing grid starting at zero."
        )
    if not np.all(np.isfinite(wk)) or np.any(wk < 0.0):
        raise ValueError("Normal-mode frequencies must be finite and nonnegative.")
    alpha = np.zeros(wk.shape, dtype=float)
    for index in np.ndindex(wk.shape):
        frequency = wk[index]
        if frequency == 0.0:
            continue
        integrand = spectrum / (1.0 + (omega / frequency) ** 2)
        interpolant = PchipInterpolator(omega, integrand, extrapolate=False)
        alpha[index] = (2.0 / np.pi) * interpolant.integrate(0.0, omega[-1])
    return alpha

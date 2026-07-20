
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

# TODO: Decide on representation (normal mode etc) for all branches
# debug_mf_mode is a physical CL mean-field/counterterm toggle: on/off.




import json
import numpy as np
from scipy.linalg import expm, cholesky

from ipi.engine.motion import Motion
from ipi.engine.normalmodes import NormalModes
from ipi.engine.beads import Beads

from ipi.utils import nmtransform
from ipi.utils.depend import depend_value,  dproperties
from ipi.utils.messages import info, verbosity, warning

# Markovian static
def _apply_mass_scaled_ou(p, sm, drift, noise_scale, noise):
    """Applies an exact OU step in mass-scaled coordinates."""

    p_ms = p * (drift / sm)
    deltah = np.sum((p_ms * p_ms) / (drift * drift))
    p_ms += noise_scale * noise
    deltah -= np.sum(p_ms * p_ms)
    return p_ms * sm, 0.5 * deltah

# Markovian variable
def _apply_mass_scaled_matrix_ou(s, A, dt, kbt, noise):
    """Applies an exact matrix OU step in mass-scaled coordinates.
    We diagonalise to independent modes similar to R. J. Maurer et al, PRL, 2017 """

    et = 0.5 * float(np.dot(s, s))

    A = 0.5 * (A + A.T) # Confusing notation. The A matrix here is not the same A matrix defined in GLE case. 
    evals, evecs = np.linalg.eigh(A) # This is mode diagonalization to treat each mode seperately as independent OU process.

    evals = np.clip(evals, 0.0, None)
    c = np.exp(-evals * dt)
    s2 = np.sqrt(1.0 - c * c)
    y = evecs.T @ s
    y = c * y + np.sqrt(kbt) * s2 * noise
    s = evecs @ y

    et -= 0.5 * float(np.dot(s, s))
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

    Markovian friction is treated as the zero-auxiliary limit of the GLE bath.
    Non-Markovian auxiliary variables are kept separate from the physical
    system momentum: unlike `ThermoGLE`, `self.s` is reserved for auxiliary
    bath variables only.

    TODO:  Markovian should run through identical treatment as Non-markovian
    """

    def __init__(self):
        super().__init__()
        self.s = None
        self.ns = 0
        self.theta = None
        self.A_aux = None
        self.T_aux = None
        self.S_aux = None

    def bind(self, friction, motion: Motion | None) -> None:
        super().bind(friction, motion)
        self.initialize_state()

    def initialize_state(self) -> None:
        if self.is_markovian():
            self.s = np.zeros((0,), dtype=float)
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
        self.theta = np.asarray(self.friction.Ap[0, 1:], dtype=float).copy()
        self.A_aux = np.asarray(self.friction.Ap[1:, 1:], dtype=float).copy()
        self.ns = int(self.theta.size)
        self.T_aux = None
        self.S_aux = None
        self._nm_transform_matrix = None
        if self.friction.variable_friction:
            self.s = None
            return

        self.s = np.zeros(self.state_shape(), dtype=float)
        if self.friction.prng is not None and self.ns > 0:
            self.s[:] = np.sqrt(self.friction._kbt_rp()) * self.friction.prng.gvec(
                self.s.shape
            )

    def is_markovian(self) -> bool:
        return str(self.friction.bath_mode) == "markovian"

    def state_shape(self):
        if self.is_markovian():
            return (0,)
        if self.friction.variable_friction:
            return None if self.s is None else self.s.shape
        return (
            int(self.friction.beads.nbeads),
            int(self.friction.Ap.shape[0] - 1),
            3 * int(self.friction.beads.natoms),
        )

    def _get_non_markovian_jacobian(self) -> np.ndarray:
        coupling_jacobian = np.asarray(self.friction._get_coupling_jacobian(), dtype=float)
        if coupling_jacobian.ndim != 3:
            raise ValueError(
                "Variable non-markovian friction requires coupling_jacobian "
                "with shape (nbeads, nchannels, ndof)."
            )
        expected_shape = (
            int(self.friction.beads.nbeads),
            int(self.friction.Ap.shape[0] - 1),
            int(coupling_jacobian.shape[1]),
        )
        if self.s is None:
            self.s = np.zeros(expected_shape, dtype=float)
            if self.friction.prng is not None and self.ns > 0:
                self.s[:] = np.sqrt(self.friction._kbt_rp()) * self.friction.prng.gvec(
                    self.s.shape
                )
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
            self._nm_transform_matrix = np.asarray(nmtransform.mk_nm_matrix(nbeads), dtype=float)
        return self._nm_transform_matrix

    def _get_non_markovian_nm_jacobian(self) -> np.ndarray:
        """Return dF_nm[n']/dQ_nm[n] from the bead-space Jacobian."""
        coupling_jacobian = self._get_non_markovian_jacobian()
        cmat = self._get_nm_transform_matrix()
        return np.einsum("rb,nb,bci->rnci", cmat, cmat, coupling_jacobian)

    def step(self, pdt: float) -> None:
        if self.is_markovian():
            self._step_markovian(pdt)
            return
        self._step_non_markovian(pdt)

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
            raise ValueError("gamma must be non-negative for Markovian linear coupling.")
        if gamma == 0.0:
            return

        m = friction.nm.dynm3.copy()
        p = friction.nm.pnm.copy()
        sm = np.sqrt(m)
        gamma_nm = np.full_like(m, gamma)
        drift = np.exp(-(gamma_nm / m) * pdt)
        noise_scale = np.sqrt(friction._kbt_rp() * (1.0 - drift * drift))
        p_new, bath_energy = _apply_mass_scaled_ou(
            p, sm, drift, noise_scale, friction.prng.gvec(p.shape)
        )
        friction.nm.pnm[:] = p_new
        friction.beads.p = friction.nm.transform.nm2b(p_new)
        friction.ediss += bath_energy

    def _step_markovian_variable(self, pdt: float) -> None:
        friction = self.friction
        coupling_jacobian = friction._get_coupling_jacobian()
        nbeads = coupling_jacobian.shape[0]
        gamma = np.asarray(friction.gamma, dtype=float)
        p = friction.beads.p
        m = friction.beads.m3
        sm = np.sqrt(m)
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
            p[b, :] = s * sm[b, :]

        friction.beads.p[:] = p
        friction.ediss += et

    def _step_non_markovian(self, pdt: float) -> None:
        # Symmetric non-Markovian bath splitting scaffold:
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
        self.T_aux, self.S_aux = _compute_aux_ou_matrices(
            self.A_aux, pdt, self.friction._kbt_rp()
        )
        et = 0.5 * float(np.sum(self.s * self.s))
        # s has shape (nmodes, naux, nchannel). Apply the same auxiliary OU
        # embedding independently to each normal mode and trailing channel.
        self.s[:] = np.einsum("ab,mbd->mad", self.T_aux, self.s)
        noise = self.friction.prng.gvec(self.s.shape)
        self.s[:] += np.einsum("ab,mbd->mad", self.S_aux, noise)
        et -= 0.5 * float(np.sum(self.s * self.s))
        self.friction.ediss += et

    def bp_f_step(self, pdt: float) -> None:
        """Bath (aux) -to- physical momentum coupling, Eq. (S40).

        Implements the momentum kick generated by the current auxiliary bath
        state through:
            P <- P - dt * dF/dQ * theta^T s

        """

        if pdt <= 0.0 or self.ns == 0:
            return
        p = self.friction.nm.pnm.copy()
        if self.friction.variable_friction:
            nm_jacobian = self._get_non_markovian_nm_jacobian()
            theta_s = np.einsum("a,rac->rc", self.theta, self.s)
            p_new = p - pdt * np.einsum("rc,rnci->ni", theta_s, nm_jacobian)
        else:
            m = self.friction.nm.dynm3.copy()
            sm = np.sqrt(m)
            sigma = float(self.friction.sigma_static)
            theta_s = np.einsum("a,mad->md", self.theta, self.s)
            p_ms = p / sm
            p_ms -= pdt * (sigma / sm) * theta_s
            p_new = p_ms * sm
        self.friction.nm.pnm[:] = p_new
        self.friction.beads.p = self.friction.nm.transform.nm2b(p_new)

    def bs_step(self, pdt: float) -> None:
        """Physical momentum-to-bath (aux) coupling, Eq. (S42).

        Implements the auxiliary update driven by physical momentum through:
            s <- s + dt * dF/dQ * theta * P


        """

        if pdt <= 0.0 or self.ns == 0:
            return
        p = self.friction.nm.pnm.copy()
        m = self.friction.nm.dynm3.copy()
        if self.friction.variable_friction:
            nm_jacobian = self._get_non_markovian_nm_jacobian()
            drive = np.einsum("rnci,ni->rc", nm_jacobian, p / m)
            self.s[:] += pdt * self.theta[None, :, None] * drive[:, None, :]
        else:
            sigma = float(self.friction.sigma_static)
            drive = sigma * p / m
            self.s[:] += pdt * self.theta[None, :, None] * drive[:, None, :]


class Friction:
    """
    Friction operator for friction-enabled dynamics.
    """

    # -------------------------
    # Input-configured
    # -------------------------
    variable_friction: bool
    bath_mode: str       # "none" | "markovian" | "non-markovian"
    debug_mf_mode: str         # "on" | "off"

    Lambda: np.ndarray  # [omega, J(omega)] for non-markovian OU fit
    Ap: np.ndarray  # George-style momentum + auxiliary drift matrix
    debug_alpha_input: np.ndarray      # optional [omega_k, alpha]

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
        bath_mode: str = "non-markovian", # can be 1. none (no dissipative, no random force),  
        # todo:Make a boolean, only_conservative - true or false. 

        #todo: block that is identity then it si automatically markovian.  Ap

        debug_mf_mode: str = "on",
        # on  - apply Caldeira-Leggett mean-field/counterterm.
        # off - no friction contribution to conservative force.

        Lambda=np.zeros((0, 2), float),
        #todo: switch back to spectral density.  with some extrapolation to zero. use cubic spline and linear extrapolation to zero. and check
        # rename spectral_density
        Ap=np.zeros((0, 0), float),
        # George-style auxiliary drift matrix:
        #   Ap = [[0, theta^T], [-theta, A]]
        # The current non-Markovian implementation uses this direct embedding
        # and does not yet fit Ap from Lambda.

        debug_alpha_input=np.zeros((0, 2), float),

        sigma_static: float = 1.0,
        # if vartiable_friction is false.. then gamma = s * s   (s is a float)

        coupling_jacobian_key: str = "friction_coupling_jacobian",
        gamma_key: str = "friction_gamma",
        friction_meta_key: str = "friction_meta",
        coupling_key: str = "friction_coupling",
        coupling_mode: str = "driver",
        centroid_coupling_jacobian_key: str = "centroid_friction_coupling_jacobian",
        centroid_friction_meta_key: str = "centroid_friction_meta",
        coupling_friction_atom: int = -1,
    ):
        """Initialises the friction object.
        Args:
            Lambda: Cosine transform of the time-dependent factor in the friction kernel,
                divided by frequncy. Supplied as a 2d array of two columns containing frequency and
                spectral density, respectively.
                Defaults to np.zeros(0, float).
            debug_alpha_input: Normal-mode coefficients in expression for the frictional mean-field
                potential [Eq. (8b) in https://doi.org/10.1103/PhysRevLett.134.226201].
                Defaults to np.zeros(0, float).
            variable_friction (bool, optional): True if the gradient of the friction coupling F(q)
                [introduced in Eq. (5) of https://doi.org/10.1103/PhysRevLett.134.226201]
                depends on position.
                Defaults to False.
        """
        #todo: more descriptive.

        # Choices
        self.variable_friction = bool(variable_friction)
        self.bath_mode = str(bath_mode)
        self.debug_mf_mode = self._normalize_debug_mf_mode(debug_mf_mode)

        # Kernel shape
        self.Lambda = np.asanyarray(Lambda, dtype=float).copy()
        self.Ap = np.asanyarray(Ap, dtype=float).copy()
        self.debug_alpha_input = np.asanyarray(debug_alpha_input, dtype=float).copy()

        self._coupling_jacobian = depend_value(
            name="coupling_jacobian", func=self._get_coupling_jacobian
        )
        self._gamma = depend_value(
            name="gamma", func=self._get_gamma, dependencies=[self._coupling_jacobian]
        )
    
    #     # Friction coupling: F(q), such that Σ{i,α} = ∂F(q) / ∂q{i,α}
        self._friction_coupling_nm = depend_value(
            name="friction_coupling_nm",
            func=self.get_friction_coupling_nm,
            dependencies=[self._coupling_jacobian],
        )
        # Frictional mean-field force
        self._force_mf_nm = depend_value(
            name="force_mf_nm",
            func=self.get_force_mf_nm,
            dependencies=[self._friction_coupling_nm],
        )

        #force_meanfield
        self._force_mf = depend_value(
            name="force_mf", func=self.get_force_mf, dependencies=[self._force_mf_nm]
        )

        #Conserved mean-field potential 
        self._energy_mf = depend_value(name="energy_mf", value=0.0)

        self.sigma_static = float(sigma_static)

        self.coupling_jacobian_key = str(coupling_jacobian_key)
        self.gamma_key = str(gamma_key)
        self.friction_meta_key = str(friction_meta_key)
        self.coupling_key = str(coupling_key)
        self.coupling_mode = str(coupling_mode)
        self.centroid_coupling_jacobian_key = str(centroid_coupling_jacobian_key)
        self.centroid_friction_meta_key = str(centroid_friction_meta_key)
        self.coupling_friction_atom = int(coupling_friction_atom)
        self._friction_meta = {}
        self._centroid_friction_meta = {}
        self._friction_atoms_idx: np.ndarray | None = None
        self._friction_dof_idx: np.ndarray | None = None
        self._channel_labels: list[str] = []
        self._friction_contract_signature = None
        self._gamma_reconstruction_error = 0.0
        self._nm_transform_matrix = None
        self.bath: FrictionBath | None = None

        # runtime handles
        self.alpha: np.ndarray | None = None
        self.forces = None
        self.ensemble = None
        self.prng = None

        # bookkeeping: cumulative energy exchange with bath (Markovian)
        self._ediss = depend_value(name="ediss", value=0.0)     # positive = system -> bath via friction



    # ==========================================================================
    # bind
    # ==========================================================================

    def bind(self, motion: Motion) -> None:
        self.beads = motion.beads
        self.nm = motion.nm
        self.ensemble = motion.ensemble
        self.forces = motion.ensemble.forces
        self.prng = motion.prng
        if self.bath_mode not in ("none", "markovian", "non-markovian"):
            raise ValueError("bath_mode must be one of: 'none', 'markovian', 'non-markovian'.")

        self.debug_mf_mode = self._normalize_debug_mf_mode(self.debug_mf_mode)
        if self.debug_mf_mode not in ("on", "off"):
            raise ValueError("debug_mf_mode must be one of: 'on', 'off'.")

        if self.coupling_mode not in ("driver", "centroid_endpoint_trapezoid"):
            raise ValueError(
                "coupling_mode must be one of: 'driver', 'centroid_endpoint_trapezoid'."
            )


        # Non-Markovian OU fit requires spectral density
        if self.bath_mode == "non-markovian":
            if self.Ap.size == 0:
                raise ValueError(
                    "non-markovian requires an explicit Ap matrix. "
                    "Fitting Ap from Lambda is not implemented yet."
                )
            self._validate_Ap()

        # Setup alpha for MF (may be zeros if MF disabled)
        self.alpha = self._setup_alpha()

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
        self._coupling_jacobian.add_dependency(self.forces._extras)
        self._friction_coupling_nm.add_dependency(self.beads._q)
        self._energy_mf.add_dependency(self._friction_coupling_nm)
        self._energy_mf._func = self.get_energy_mf

    def _build_bath(self) -> FrictionBath | None:
        if self.bath_mode == "none":
            return None
        if self.bath_mode in ("markovian", "non-markovian"):
            return FrictionGLE()
        raise RuntimeError("bath_mode must be one of none, markovian or non-markovian")

    def _validate_Ap(self) -> None:
        Ap = np.asarray(self.Ap, dtype=float)
        if Ap.ndim != 2 or Ap.shape[0] != Ap.shape[1] or Ap.shape[0] < 2:
            raise ValueError(
                "Ap must be a square 2D matrix with shape (1+naux, 1+naux)."
            )
        A = Ap[1:, 1:]
        if A.size == 0:
            raise ValueError("Ap must contain at least one auxiliary variable.")
        if not np.all(np.isfinite(Ap)):
            raise ValueError("Ap contains non-finite values.")
        if not np.isclose(Ap[0, 0], 0.0, rtol=0.0, atol=1e-14):
            raise ValueError("Ap[0,0] must be zero for the CL/GLE convention.")
        if not np.allclose(Ap[1:, 0], -Ap[0, 1:], rtol=1e-10, atol=1e-14):
            raise ValueError(
                "Ap must follow the CL/GLE convention Ap[1:,0] == -Ap[0,1:]."
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
    
    def _kbt(self) -> float:
        from ipi.utils.units import Constants
        return Constants.kb * float(self.ensemble.temp)
    
    @staticmethod
    def _kinetic(p: np.ndarray, m3: np.ndarray) -> float:
        # p and m3 same shape
        return 0.5 * float(np.sum((p * p) / m3))



    # ==========================================================================
    # Alpha setup (MF)
    # ==========================================================================

    def _setup_alpha(self) -> np.ndarray:
        wk = np.asarray(self.nm.omegak, dtype=float)
        nmodes = wk.size

        # If MF disabled, keep alpha zeros.
        if self.debug_mf_mode == "off":
            return np.zeros(nmodes, dtype=float)

        # Ap is the physical Caldeira-Leggett/GLE input. Development-only
        # debug_alpha_input/Lambda fallbacks are used only when Ap is absent.
        if self.Ap.size > 0 and self.bath_mode == "non-markovian":
            alpha = self._alpha_from_Ap(wk)
            info("Friction: computed alpha^(n) from Ap.", verbosity.low)
            return alpha

        # debug_alpha_input is a development fallback.
        if self.debug_alpha_input.size > 0:
            if self.debug_alpha_input.ndim != 2 or self.debug_alpha_input.shape[1] != 2:
                raise ValueError("debug_alpha_input must have shape (nmodes,2) [omega_k, alpha].")
            if self.debug_alpha_input.shape[0] != nmodes:
                raise ValueError(f"debug_alpha_input rows ({self.debug_alpha_input.shape[0]}) != nmodes ({nmodes}).")
            if not np.allclose(self.debug_alpha_input[:, 0], wk):
                raise ValueError("debug_alpha_input omega_k does not match current nm.omegak.")
            info("Friction: using alpha from debug_alpha_input table.", verbosity.low)
            return np.asarray(self.debug_alpha_input[:, 1], dtype=float)

        # Lambda is a development fallback for old workflows.
        if self.Lambda.size > 0:
            omega = np.asarray(self.Lambda[:, 0], dtype=float)
            if omega.size < 2:
                raise ValueError("Lambda must contain at least two points.")
            if np.any(omega <= 0.0) or np.any(np.diff(omega) <= 0.0):
                raise ValueError("Lambda omega must be strictly positive and increasing.")
            alpha = get_alpha_numeric(Lambda=self.Lambda[:,1], omega=omega, omegak=wk)
            info("Friction: computed alpha^(n) from Lambda.", verbosity.low)
            return alpha

        raise ValueError(
            "Friction: MF requested but no Ap, debug_alpha_input, or Lambda was supplied."
        )

    def _alpha_from_Ap(self, omegak: np.ndarray) -> np.ndarray:
        """Compute CL counterterm coefficients from the Ap GLE kernel."""
        Ap = np.asarray(self.Ap, dtype=float)
        theta = np.asarray(Ap[0, 1:], dtype=float)
        A = np.asarray(Ap[1:, 1:], dtype=float)
        eye = np.eye(A.shape[0], dtype=float)
        alpha = np.zeros(np.asarray(omegak, dtype=float).shape, dtype=float)
        for i, wk in enumerate(np.asarray(omegak, dtype=float)):
            if np.isclose(wk, 0.0, rtol=0.0, atol=1e-14):
                alpha[i] = 0.0
                continue
            alpha[i] = float(wk * theta @ np.linalg.solve(A + wk * eye, theta))
        if not np.all(np.isfinite(alpha)):
            raise ValueError("Alpha derived from Ap contains non-finite values.")
        if np.any(alpha < -1e-12):
            raise ValueError(f"Alpha derived from Ap contains negative values: {alpha}")
        alpha[alpha < 0.0] = 0.0
        return alpha

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
        meta = self.forces.extras.get(key)
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
                raise ValueError(f"Canonical friction metadata '{key}' differs across beads.")
            meta = first
        if not isinstance(meta, dict):
            raise ValueError(f"Canonical friction metadata '{key}' must be a dictionary.")
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
            dofs = np.concatenate(
                [np.arange(3 * atom, 3 * atom + 3, dtype=int) for atom in atoms]
            ) if atoms.size else np.zeros(0, dtype=int)
        else:
            raise ValueError(
                f"'{self.friction_meta_key}' must define zero-based active_dofs or active_atoms."
            )
        if dofs.size == 0 or np.any(dofs < 0) or np.any(dofs >= ndof):
            raise ValueError(f"Active Cartesian DOFs must be unique indices in [0, {ndof - 1}].")
        if np.unique(dofs).size != dofs.size:
            raise ValueError("Canonical friction metadata contains duplicate active DOFs.")
        if np.any(atoms < 0) or np.any(atoms >= natoms):
            raise ValueError(f"Active atoms must be zero-based indices in [0, {natoms - 1}].")
        self._friction_dof_idx = dofs
        self._friction_atoms_idx = atoms

    def _embed_jacobian(self, reduced: np.ndarray) -> np.ndarray:
        ndof = 3 * int(self.beads.natoms)
        if reduced.shape[-1] == ndof:
            return reduced
        if self._friction_dof_idx is None or reduced.shape[-1] != len(self._friction_dof_idx):
            raise ValueError(
                f"Reduced coupling Jacobian shape {reduced.shape} is inconsistent with active DOFs."
            )
        full = np.zeros(reduced.shape[:-1] + (ndof,), dtype=float)
        full[..., self._friction_dof_idx] = reduced
        return full

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
            full[bead][np.ix_(self._friction_dof_idx, self._friction_dof_idx)] = reduced[bead]
        return full

    def _get_coupling_jacobian(self) -> np.ndarray:
        """Return dF_c(Q_b)/dQ_i as (nbeads, nchannels, full_ndof)."""
        if not self.variable_friction:
            return float(self.sigma_static)
        meta = self._get_consistent_meta(self.friction_meta_key, bead_resolved=True)
        self._configure_active_dofs(meta)
        labels = meta.get("channel_labels")
        if not isinstance(labels, list) or not all(isinstance(label, str) for label in labels):
            raise ValueError(f"'{self.friction_meta_key}.channel_labels' must be a list of strings.")
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
        payload = self._decode_json(self.forces.extras.get(self.coupling_jacobian_key))
        if payload is None:
            raise KeyError(f"Missing canonical driver extra '{self.coupling_jacobian_key}'.")
        jacobian = np.asarray(payload, dtype=float)
        expected_prefix = (int(self.beads.nbeads), len(labels))
        if jacobian.ndim != 3 or jacobian.shape[:2] != expected_prefix:
            raise ValueError(
                f"'{self.coupling_jacobian_key}' must have shape "
                f"(nbeads, nchannels, nactive_dof); got {jacobian.shape}."
            )
        if not np.all(np.isfinite(jacobian)):
            raise ValueError(f"'{self.coupling_jacobian_key}' contains non-finite values.")
        self._friction_meta = meta
        self._channel_labels = list(labels)
        return self._embed_jacobian(jacobian)

    def _get_centroid_coupling_jacobian(self) -> np.ndarray:
        payload = self._decode_json(
            self.forces.extras.get(self.centroid_coupling_jacobian_key)
        )
        if payload is None:
            raise KeyError(
                f"coupling_mode='centroid_endpoint_trapezoid' requires "
                f"'{self.centroid_coupling_jacobian_key}'."
            )
        meta = self._get_consistent_meta(
            self.centroid_friction_meta_key, bead_resolved=False
        )
        for key in ("active_atoms", "active_dofs", "channel_labels"):
            if meta.get(key) != self._friction_meta.get(key):
                raise ValueError(f"Centroid friction metadata mismatch for '{key}'.")
        jacobian = np.asarray(payload, dtype=float)
        if jacobian.ndim == 3 and jacobian.shape[0] == 1:
            jacobian = jacobian[0]
        if jacobian.ndim != 2 or jacobian.shape[0] != len(self._channel_labels):
            raise ValueError(
                f"'{self.centroid_coupling_jacobian_key}' must have shape "
                f"(nchannels, nactive_dof); got {jacobian.shape}."
            )
        if not np.all(np.isfinite(jacobian)):
            raise ValueError(
                f"'{self.centroid_coupling_jacobian_key}' contains non-finite values."
            )
        self._centroid_friction_meta = meta
        return self._embed_jacobian(jacobian[np.newaxis, ...])

    def _get_gamma(self):
        """Return driver-supplied Gamma after validating its canonical factor."""
        if not self.variable_friction:
            return float(self.sigma_static) ** 2
        jacobian = np.asarray(self.coupling_jacobian, dtype=float)
        payload = self._decode_json(self.forces.extras.get(self.gamma_key))
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
        return gamma
    

    # ==========================================================================
    # Friction forces and coupling
    # ==========================================================================

    def _get_nm_transform_matrix(self) -> np.ndarray:
        nbeads = int(self.beads.nbeads)
        if self._nm_transform_matrix is None or self._nm_transform_matrix.shape != (
            nbeads,
            nbeads,
        ):
            self._nm_transform_matrix = np.asarray(nmtransform.mk_nm_matrix(nbeads), dtype=float)
        return self._nm_transform_matrix

    def _get_friction_coupling(self) -> np.ndarray:
        if self.coupling_mode == "centroid_endpoint_trapezoid":
            return self._get_centroid_endpoint_trapezoid_coupling()

        coupling = self.forces.extras.get(self.coupling_key)
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
        return coupling

    def _infer_coupling_atom(self) -> int:
        """Return 0-based atom index for centroid-relative displacement."""
        if self.coupling_friction_atom >= 0:
            if self.coupling_friction_atom >= int(self.beads.natoms):
                raise ValueError(
                    f"coupling_friction_atom={self.coupling_friction_atom} outside [0, {self.beads.natoms - 1}]"
                )
            return int(self.coupling_friction_atom)

        if self._friction_atoms_idx is None or self._friction_atoms_idx.size != 1:
            raise ValueError(
                "coupling_mode='centroid_endpoint_trapezoid' requires exactly one friction atom "
                "in friction_meta.active_atoms, or an explicit coupling_friction_atom."
            )
        return int(self._friction_atoms_idx[0])

    def _get_centroid_endpoint_trapezoid_coupling(self) -> np.ndarray:
        """Compute the endpoint-trapezoid coupling for each bead and channel."""
        jacobian_b = np.asarray(self.coupling_jacobian, dtype=float)
        jacobian_c = np.asarray(self._get_centroid_coupling_jacobian(), dtype=float)
        if jacobian_b.ndim != 3:
            raise ValueError(
                "centroid_endpoint_trapezoid requires a channel-resolved coupling Jacobian."
            )
        if jacobian_c.shape[1:] != jacobian_b.shape[1:]:
            raise ValueError(
                f"Centroid Jacobian shape {jacobian_c.shape} incompatible with "
                f"bead Jacobian shape {jacobian_b.shape}."
            )

        atom = self._infer_coupling_atom()
        dof = slice(3 * atom, 3 * atom + 3)
        dq = np.asarray(self.beads.q[:, dof] - self.beads.qc[dof], dtype=float).copy()
        if hasattr(self.forces, "cell") and self.forces.cell is not None:
            dq_flat = dq.reshape(-1).copy()
            self.forces.cell.array_pbc(dq_flat)
            dq = dq_flat.reshape((-1, 3))

        jacobian_h = jacobian_b[:, :, dof]
        jacobian_ch = jacobian_c[0, :, dof]
        coupling = 0.5 * np.einsum(
            "bci,bi->bc", jacobian_h + jacobian_ch[np.newaxis, :, :], dq
        )
        if not np.all(np.isfinite(coupling)):
            raise ValueError("centroid_endpoint_trapezoid coupling contains non-finite values.")
        return coupling

    def _get_nm_coupling_jacobian(self) -> np.ndarray:
        coupling_jacobian = np.asarray(self._get_coupling_jacobian(), dtype=float)
        cmat = self._get_nm_transform_matrix()
        return np.einsum("rb,nb,bci->rnci", cmat, cmat, coupling_jacobian)

    def get_friction_coupling_nm(self):
        """Compute the friction coupling for each normal-mode index"""
        if self.variable_friction:
            # In the variable friction case, we must get the coupling F(Q) from the driver. A future example of this could be 
            # an ML model that diretly provides the coupling for atomistic systems.
            return self._get_nm_transform_matrix() @ self._get_friction_coupling()
        else:
            # Here we assume that the interaction potential, F(Q) in https://doi.org/10.1103/PhysRevLett.134.226201,
            # is of the form F(q) = SUM[ c{i,α} q{i,α}, {{i,0,n_atom-1}, {α,0,2}} ] where α indexes Cartesian components
            # The diffusion coefficients for bead index l returned by the driver are expected to be packed as
            # Σ{i,α} = ∂F(q) / ∂q{i,α} = diffusion_coeff[l, 3*i+α].
            return np.sum(float(self.sigma_static) * self.nm.qnm, axis=-1)

    def get_energy_mf(self):
        """Compute the frictional potential of mean field, Eq. (S19) of https://doi.org/10.1103/PhysRevLett.134.226201"""

        if self.debug_mf_mode == "off":
            return 0.0

        coupling_nm = np.asarray(self.friction_coupling_nm, dtype=float)
        if coupling_nm.ndim == 1:
            weighted_coupling2 = self.alpha * coupling_nm**2
        else:
            weighted_coupling2 = self.alpha[:, np.newaxis] * coupling_nm**2
        return np.sum(weighted_coupling2) / 2
        
    def get_force_mf_nm(self):
        """Negative derivative of the frictional potential of mean field with respect to normal modes"""
        if self.variable_friction:
            return -np.einsum(
                "r,rc,rnci->ni",
                self.alpha,
                np.asarray(self.friction_coupling_nm, dtype=float),
                self._get_nm_coupling_jacobian(),
            )

        return -(self.alpha * self.friction_coupling_nm)[:, np.newaxis] * float(
            self.sigma_static
        )

    def get_force_mf(self):
        """Negative derivative of the frictional potential of mean field with respect to bead positions"""
        return self.nm.transform.nm2b(self.force_mf_nm)

    # ==========================================================================
    # main step
    # ==========================================================================

    def step(self, pdt: float) -> None:
        """
        Apply friction operator over time interval pdt:
          - MF kick (if enabled)
          - bath kick (markovian/non-markovian) - dissipative and random forces
        """

        # MF
        if self.debug_mf_mode == "on":
            self.beads.p += self.force_mf * pdt

        self._ensure_bath_bound()
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



def get_alpha_numeric(Lambda: np.ndarray, omega: np.ndarray, omegak: np.ndarray) -> np.ndarray:
    """Numerically compute alpha^(n) from Lambda(omega)."""
    try:
        from scipy.interpolate import CubicSpline
    except ModuleNotFoundError as e:
        raise ModuleNotFoundError(
            "Friction: scipy is required to compute alpha from Lambda. "
            "Install scipy or provide debug_alpha_input explicitly."
        ) from e

    omega = np.asarray(omega, dtype=float)
    Lambda = np.asarray(Lambda, dtype=float)
    omegak = np.asarray(omegak, dtype=float)

    alpha = np.zeros(omegak.shape, dtype=float)
    for i, wk in enumerate(omegak):
        f = CubicSpline(omega, Lambda * (wk**2) / (omega**2 + wk**2))
        alpha[i] = (2.0 / np.pi) * f.integrate(0.0, omega[-1])
        info(f"Friction: wk={wk} alpha={alpha[i]}", verbosity.high)

    # alpha = np.zeros(omegak.shape)
    # for idx, omegak in enumerate(omegak):
    #     # TODO: what if omega[0] > 0?
    #     f = CubicSpline(omega, Lambda * omegak**2 / (omega**2 + omegak**2))
    #     alpha[idx] = 2 / np.pi * f.integrate(0, omega[-1])
    #     info(
    #         f"for normal mode {omegak} alpha is {alpha[idx]}",
    #         verbosity.high,
    #     )
    return alpha

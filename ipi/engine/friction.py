
"""Electronic-friction operators for friction-enabled dynamics in i-PI.

This module separates electronic friction into two layers.

`Friction`
    High-level operator owned by the `nve-f` / `nvt-f` integrators. It parses
    force-driver extras, constructs `sigma` and `gamma`, manages metadata such
    as `sigma_meta`, and evaluates the optional friction mean-field (MF)
    contribution (`friction_coupling_nm`, `force_mf`, `energy_mf`).

`FrictionGLE`
    Bath propagator used by `Friction.step(pdt)`. In the current implementation
    it covers the Markovian electronic-friction bath for both:
      - static friction (`sigma_static`)
      - variable friction (`sigma(q)` from driver extras)
    and Non-Markovian for both static and variable friction.

Current Markovian/Non-Markovian behavior:
    1. `Friction.step(pdt)` applies the MF momentum kick, if enabled.
    2. `FrictionGLE.step(pdt)` applies the dissipative/stochastic bath update.

For variable friction, drivers provide `sigma`; i-PI converts it to
`gamma = sigma^T sigma` and uses an exact frozen-geometry OU update over the
current substep. When the driver only returns friction-active atoms,
`sigma_meta["friction_atoms"]` is used to embed the reduced matrix into the
full Cartesian system (n_atoms).
"""

# TODO: Decide on representation (normal mode etc) for all branches
# TODO: debug_mf_mode has outdated options (e.g linear),  for now we assume the driver will provide friction coupling
# so that debug_mf_mode will just be on or off (default on)




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

    def _get_non_markovian_sigma(self) -> np.ndarray:
        sigma = np.asarray(self.friction._get_sigma(), dtype=float)
        if sigma.ndim != 3:
            raise ValueError(
                "Variable non-markovian friction requires sigma with shape "
                "(nbeads, nbath, ndof)."
            )
        sigma_mode = str(self.friction._sigma_meta.get("sigma_mode", "column")).lower()
        if sigma_mode == "row":
            if self.friction._sigma_blocks is None:
                raise ValueError(
                    "Variable non-markovian row-mode friction requires a dict "
                    "sigma payload so row-packed blocks are available."
                )
            sigma = np.asarray(
                [
                    np.concatenate(
                        [np.asarray(m, dtype=float).T for m in mats], axis=0
                    )
                    for mats in self.friction._sigma_blocks
                ],
                dtype=float,
            )
            sigma = self.friction._embed_if_needed(
                sigma, ndof=3 * int(self.friction.beads.natoms)
            )
        elif sigma_mode == "pairwise":
            raise ValueError(
                "Variable non-markovian pairwise sigma_mode is not implemented. "
                "Use markovian friction or provide a column/row-packed Sigma."
            )
        elif sigma_mode != "column":
            raise ValueError(
                f"Unsupported {self.friction.sigma_meta_key}.sigma_mode='{sigma_mode}'. "
                "Supported non-markovian values are 'column' and 'row'."
            )
        expected_shape = (
            int(self.friction.beads.nbeads),
            int(self.friction.Ap.shape[0] - 1),
            int(sigma.shape[1]),
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
                f"current sigma payload. Expected {expected_shape}, got {self.s.shape}."
            )
        return sigma

    def _get_nm_transform_matrix(self) -> np.ndarray:
        nbeads = int(self.friction.beads.nbeads)
        if self._nm_transform_matrix is None or self._nm_transform_matrix.shape != (
            nbeads,
            nbeads,
        ):
            self._nm_transform_matrix = np.asarray(nmtransform.mk_nm_matrix(nbeads), dtype=float)
        return self._nm_transform_matrix

    def _get_non_markovian_nmdsigma(self) -> np.ndarray:
        """Returns dF_nm[n'] / dQ_nm[n] from bead-space sigma=dF/dq."""
        sigma = self._get_non_markovian_sigma()
        cmat = self._get_nm_transform_matrix()
        return np.einsum("rb,nb,bci->rnci", cmat, cmat, sigma)

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
        sigma = friction._get_sigma()
        nbeads = sigma.shape[0]
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
            self._get_non_markovian_sigma()
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
            nmdsigma = self._get_non_markovian_nmdsigma()
            theta_s = np.einsum("a,rac->rc", self.theta, self.s)
            p_new = p - pdt * np.einsum("rc,rnci->ni", theta_s, nmdsigma)
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
            nmdsigma = self._get_non_markovian_nmdsigma()
            drive = np.einsum("rnci,ni->rc", nmdsigma, p / m)
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
    debug_mf_mode: str         # "none" | "linear"

    Lambda: np.ndarray  # [omega, J(omega)] for non-markovian OU fit
    Ap: np.ndarray  # George-style momentum + auxiliary drift matrix
    debug_alpha_input: np.ndarray      # optional [omega_k, alpha]

    sigma_static: float

    # Extras parsing
    sigma_key: str

    # -------------------------
    # Runtime bound
    # -------------------------
    beads: Beads
    """Reference to the beads"""
    nm: NormalModes
    """Reference to the normal modes"""

    def __init__(
        self,
        variable_friction: bool = True,   #Variable_friction true means sigma changes with position. Otherwise use static_sigma
        bath_mode: str = "non-markovian", # can be 1. none (no dissipative, no random force),  
        # todo:Make a boolean, only_conservative - true or false. 

        #todo: block that is identity then it si automatically markovian.  Ap

        debug_mf_mode: str = "none",
        #none - no friction contribution to consevative force (basically same as alpha=0).
        ##linear - 

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

        sigma_key: str = "sigma", # points to dictionary key where sigma AKA diffusion coefficient is stored.
        coupling_key: str = "friction_coupling",
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
        self.debug_mf_mode = str(debug_mf_mode)

        # Kernel shape
        self.Lambda = np.asanyarray(Lambda, dtype=float).copy()
        self.Ap = np.asanyarray(Ap, dtype=float).copy()
        self.debug_alpha_input = np.asanyarray(debug_alpha_input, dtype=float).copy()

        self._sigma = depend_value(name="sigma", func=self._get_sigma)
        self._gamma = depend_value(
            name="gamma", func=self._get_gamma, dependencies=[self._sigma]
        )
    
    #     # Friction coupling: F(q), such that Σ{i,α} = ∂F(q) / ∂q{i,α}
        self._friction_coupling_nm = depend_value(
            name="friction_coupling_nm",
            func=self.get_friction_coupling_nm,
            dependencies=[self._sigma],
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

        self.sigma_key = str(sigma_key)
        self.coupling_key = str(coupling_key)
        self.sigma_meta_key = "sigma_meta"
        self._sigma_meta = {}
        self._sigma_blocks = None
        self._friction_atoms_idx: np.ndarray | None = None
        self._friction_dof_idx: np.ndarray | None = None
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

        if self.debug_mf_mode not in ("none", "linear"):
            raise ValueError("debug_mf_mode must be one of: 'none', 'linear'.")


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
            f"  sigma_key           = '{self.sigma_key}'\n",
            verbosity.low,
        )

        self.bath = self._build_bath()
        if self.bath is not None:
            self.bath.bind(self, motion)

        # Dependencies
        self._sigma.add_dependency(self.forces._extras)
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

        # debug_alpha_input overrides everything
        if self.debug_alpha_input.size > 0:
            if self.debug_alpha_input.ndim != 2 or self.debug_alpha_input.shape[1] != 2:
                raise ValueError("debug_alpha_input must have shape (nmodes,2) [omega_k, alpha].")
            if self.debug_alpha_input.shape[0] != nmodes:
                raise ValueError(f"debug_alpha_input rows ({self.debug_alpha_input.shape[0]}) != nmodes ({nmodes}).")
            if not np.allclose(self.debug_alpha_input[:, 0], wk):
                raise ValueError("debug_alpha_input omega_k does not match current nm.omegak.")
            info("Friction: using alpha from debug_alpha_input table.", verbosity.low)
            return np.asarray(self.debug_alpha_input[:, 1], dtype=float)

        # If MF disabled, keep alpha zeros
        if self.debug_mf_mode == "none":
            return np.zeros(nmodes, dtype=float)

        # If MF requested but no debug_alpha_input:
        # - for non-markovian compute from Lambda
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
            "Friction: MF requested but no debug_alpha_input (and no Lambda)."
        )

    # ==========================================================================
    # Parse Sigma
    # ==========================================================================

    def _set_friction_atoms_from_meta(self, sigma_meta: dict, natoms: int) -> None:
        atoms = sigma_meta.get("friction_atoms")
        if atoms is None:
            self._friction_atoms_idx = None
            self._friction_dof_idx = None
            return

        atoms = np.asarray(atoms, dtype=int).flatten()
        if atoms.size == 0:
            self._friction_atoms_idx = None
            self._friction_dof_idx = None
            return
        if np.any(atoms < 0) or np.any(atoms >= natoms):
            raise ValueError(
                f"{self.sigma_meta_key}.friction_atoms must be 0-based indices in [0, {natoms - 1}], got {atoms}."
            )
        if np.unique(atoms).size != atoms.size:
            raise ValueError(
                f"{self.sigma_meta_key}.friction_atoms contains duplicate indices: {atoms}"
            )

        self._friction_atoms_idx = atoms
        self._friction_dof_idx = np.concatenate(
            [np.arange(3 * a, 3 * a + 3, dtype=int) for a in atoms]
        )

    def _embed_if_needed(self, sigma3d: np.ndarray, ndof: int) -> np.ndarray:
        ndof_reduced = int(sigma3d.shape[2])
        if ndof_reduced == ndof:
            return sigma3d
        if self._friction_dof_idx is None:
            raise ValueError(
                f"{self.sigma_key} shape {sigma3d.shape} incompatible with ndof={ndof}, "
                f"and no {self.sigma_meta_key}.friction_atoms metadata was provided."
            )
        if ndof_reduced != self._friction_dof_idx.size:
            raise ValueError(
                f"Reduced {self.sigma_key} ndof={ndof_reduced} does not match "
                f"3*len({self.sigma_meta_key}.friction_atoms)={self._friction_dof_idx.size}."
            )
        sigma_full = np.zeros((sigma3d.shape[0], sigma3d.shape[1], ndof), dtype=sigma3d.dtype)
        sigma_full[:, :, self._friction_dof_idx] = sigma3d
        return sigma_full

    # AKA get_diffusion_coefficient
    def _get_sigma(self) -> np.ndarray:
        """ 

        Reads i-PI combined extras:
            extras[key] is expected to be a 
            each entry being the per-bead payload (a dict).

        Returns:
            sigma : (shape)

        """


        if (not self.variable_friction):
            return float(self.sigma_static)

        nbeads = int(self.beads.nbeads)
        ndof = 3 * int(self.beads.natoms)

        sigma = self.forces.extras.get(self.sigma_key)
        sigma_meta = self._get_sigma_meta()
        self._sigma_meta = sigma_meta
        self._sigma_blocks = None
        self._set_friction_atoms_from_meta(sigma_meta, natoms=int(self.beads.natoms))

        if sigma is None:
            raise KeyError(
                f"Did not find '{self.sigma_key}' among the force extras = {self.forces.extras}"
            )

        # Accept JSON-string payloads and decode before shape handling.
        if isinstance(sigma, str):
            try:
                sigma = json.loads(sigma)
            except json.JSONDecodeError:
                pass

        # Handle nested ACE/Julia payload:
        #   {"inv": {"1": M, ...}, "equ"/"eqv": {"1": M, ...}}
        # or per-bead list/tuple of such dicts.
        # IMPORTANT: independent representation blocks are concatenated along the
        # bath dimension (rows), not summed. Summing would introduce cross terms
        # in Sigma^T Sigma and distort Gamma.
        def _rep_mats_from_dict(dsig: dict, meta: dict | None = None) -> list[np.ndarray]:
            mats = []
            # Keep historical keys first for deterministic channel order.
            # Driver can override ordering via extras[sigma_meta_key]["rep_order"].
            rep_order = None
            if isinstance(meta, dict):
                ro = meta.get("rep_order")
                if isinstance(ro, (list, tuple)) and all(isinstance(k, str) for k in ro):
                    rep_order = list(ro)
            if rep_order is None:
                pref_keys = ("equ", "eqv", "inv", "cov")
                extra_keys = sorted([k for k in dsig.keys() if k not in pref_keys], key=lambda x: str(x))
                rep_order = list(pref_keys) + list(extra_keys)

            for rep_key in rep_order:
                rep_data = dsig.get(rep_key)
                if rep_data is None:
                    continue
                if isinstance(rep_data, dict):
                    # Preserve deterministic ordering of channels.
                    for k in sorted(rep_data.keys(), key=lambda x: str(x)):
                        mats.append(np.asarray(rep_data[k], dtype=float))
                else:
                    mats.append(np.asarray(rep_data, dtype=float))

            if len(mats) == 0:
                raise ValueError(
                    f"{self.sigma_key} dict payload does not contain 'equ'/'eqv' or 'inv' entries."
                )

            ndof0 = None
            for m in mats:
                if m.ndim != 2:
                    raise ValueError(
                        f"Each matrix in {self.sigma_key} dict payload must be 2D. Got {m.shape}."
                    )
                if ndof0 is None:
                    ndof0 = int(m.shape[1])
                elif int(m.shape[1]) != ndof0:
                    raise ValueError(
                        f"Inconsistent ndof across {self.sigma_key} dict payload: "
                        f"{[mm.shape for mm in mats]}"
                    )

            return mats

        if isinstance(sigma, dict):
            if nbeads != 1:
                raise ValueError(
                    f"{self.sigma_key} received a single dict payload but nbeads={nbeads}. "
                    f"Provide one dict per bead (list length must equal nbeads)."
                )
            mats = _rep_mats_from_dict(sigma, sigma_meta)
            sigma_eff = np.concatenate(mats, axis=0)
            sigma = sigma_eff[np.newaxis, :, :].copy()
            self._sigma_blocks = [mats]
        elif isinstance(sigma, (list, tuple)) and len(sigma) > 0 and all(
            isinstance(s, dict) for s in sigma
        ):
            if len(sigma) != nbeads:
                raise ValueError(
                    f"{self.sigma_key} list-of-dicts length {len(sigma)} incompatible with nbeads={nbeads}."
                )
            sigma_blocks = [_rep_mats_from_dict(s, sigma_meta) for s in sigma]
            sigma = np.asarray([np.concatenate(mats, axis=0) for mats in sigma_blocks], dtype=float)
            self._sigma_blocks = sigma_blocks
        else: # plain array like double well driver
            info(str(sigma), verbosity.low)
            sigma = np.asarray(sigma, dtype=float)

        
        if sigma.ndim != 3:
            raise ValueError(f"{self.sigma_key} must have ndim=3 (nbeads, nbath, ndof). Got shape {sigma.shape}.")
        if sigma.shape[0] != nbeads:
            raise ValueError(
                f"{self.sigma_key} shape {sigma.shape} incompatible with nbeads={nbeads}."
            )
        sigma = self._embed_if_needed(sigma, ndof=ndof)

        return sigma

    def _get_sigma_meta(self) -> dict:
        """Fetches optional sigma metadata dictionary from force extras."""
        meta = self.forces.extras.get(self.sigma_meta_key)
        if meta is None:
            return {}
        # Extras can be bead-resolved lists/tuples.
        if isinstance(meta, (list, tuple)):
            if len(meta) == 0:
                return {}
            nbeads = int(self.beads.nbeads)
            if len(meta) == nbeads:
                # Use bead-0 metadata; require consistency if multiple beads.
                m0 = meta[0]
                for mb in meta[1:]:
                    if type(mb) != type(m0):
                        raise ValueError(
                            f"{self.sigma_meta_key} payload types differ across beads: "
                            f"{[type(x) for x in meta]}"
                        )
                meta = m0
            else:
                # Non-bead list: try first item as a best-effort fallback.
                meta = meta[0]
        if isinstance(meta, str):
            try:
                meta = json.loads(meta)
            except json.JSONDecodeError:
                return {}
        return meta if isinstance(meta, dict) else {}

    def _get_gamma(self):
        """Returns Gamma from Sigma.

        - static friction: gamma = sigma_static^2 (scalar)
        - variable friction: Gamma[b] = Sigma[b]^T Sigma[b] (ndof x ndof)
        """
        sigma = self.sigma
        if np.isscalar(sigma):
            s = float(sigma)
            return s * s

        sarr = np.asarray(sigma, dtype=float)
        if sarr.ndim != 3:
            raise ValueError(
                f"friction.sigma has unsupported ndim={sarr.ndim}, expected scalar or 3."
            )
        # Metadata-controlled sigma mode:
        # - "column" (default): Gamma = Sigma^T Sigma from packed rows/channels.
        # - "row": Gamma = sum_k (M_k M_k^T) using per-channel blocks M_k.
        # - "pairwise": ACE PWC-style block square on 3x3 atom blocks.
        sigma_mode = str(self._sigma_meta.get("sigma_mode", "column")).lower()

        if sigma_mode in ("row", "pairwise"):
            if self._sigma_blocks is None:
                raise ValueError(
                    f"{self.sigma_meta_key}.sigma_mode='{sigma_mode}' requires sigma payload with dict blocks."
                )
            nbeads, _, ndof = sarr.shape
            gamma = np.zeros((nbeads, ndof, ndof), dtype=float)
            for b, mats in enumerate(self._sigma_blocks):
                for m in mats:
                    mm = np.asarray(m, dtype=float)

                    def _square_pairwise_block(M: np.ndarray) -> np.ndarray:
                        """ACE PWC-style square on dense 3x3 atom blocks."""
                        if M.ndim != 2 or M.shape[0] != M.shape[1] or (M.shape[0] % 3) != 0:
                            raise ValueError(
                                f"{self.sigma_meta_key}.sigma_mode='pairwise' requires square 3N x 3N blocks. "
                                f"Got shape {M.shape}."
                            )
                        nat = M.shape[0] // 3
                        G = np.zeros_like(M)
                        for i in range(nat):
                            si = slice(3 * i, 3 * i + 3)
                            for j in range(i, nat):
                                sj = slice(3 * j, 3 * j + 3)
                                sij = M[si, sj]
                                sji = M[sj, si]
                                G[si, sj] += sij @ sji.T
                                G[sj, si] += sji @ sij.T
                                G[si, si] += sij @ sij.T
                                G[sj, sj] += sji @ sji.T
                        return G

                    # Full-dof square block: use directly.
                    if mm.shape == (ndof, ndof):
                        gamma[b] += (
                            _square_pairwise_block(mm) if sigma_mode == "pairwise" else mm @ mm.T
                        )
                        continue

                    if self._friction_dof_idx is not None:
                        nred = int(len(self._friction_dof_idx))
                        if mm.shape == (nred, nred):
                            gamma[b][np.ix_(self._friction_dof_idx, self._friction_dof_idx)] += (
                                _square_pairwise_block(mm) if sigma_mode == "pairwise" else mm @ mm.T
                            )
                            continue

                    raise ValueError(
                        f"{self.sigma_meta_key}.sigma_mode='{sigma_mode}' got unsupported block shape {mm.shape}. "
                        f"Expected ({ndof},{ndof})"
                        + (
                            ""
                            if self._friction_dof_idx is None
                            else f" or ({len(self._friction_dof_idx)},{len(self._friction_dof_idx)})"
                        )
                        + "."
                    )
            return gamma

        if sigma_mode != "column":
            raise ValueError(
                f"Unsupported {self.sigma_meta_key}.sigma_mode='{sigma_mode}'. "
                "Supported values are 'column', 'row', 'pairwise'."
            )

        # (nbeads, nbath, ndof) -> (nbeads, ndof, ndof)
        return np.einsum("bai,baj->bij", sarr, sarr)
    

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

    def _get_nmdsigma(self) -> np.ndarray:
        sigma = np.asarray(self.sigma, dtype=float)
        cmat = self._get_nm_transform_matrix()
        return np.einsum("rb,nb,bci->rnci", cmat, cmat, sigma)

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
            return np.sum(self.sigma * self.nm.qnm, axis=-1)

    def get_energy_mf(self):
        """Compute the frictional potential of mean field, Eq. (S19) of https://doi.org/10.1103/PhysRevLett.134.226201"""

        if self.debug_mf_mode == "none":
            return 0.0

        # debug
        coupling_nm = np.asarray(self.friction_coupling_nm, dtype=float)
        if coupling_nm.ndim == 1:
            weighted_coupling2 = self.alpha * coupling_nm**2
        else:
            weighted_coupling2 = self.alpha[:, np.newaxis] * coupling_nm**2
        info("alpha "+str(self.alpha), verbosity.low)
        info("coupling "+str(coupling_nm**2))
        info("EMF "+str(np.sum(weighted_coupling2)), verbosity.low) 


        return np.sum(weighted_coupling2) / 2
        
    def get_force_mf_nm(self):
        """Negative derivative of the frictional potential of mean field with respect to normal modes"""
        if self.variable_friction:
            return -np.einsum(
                "r,rc,rnci->ni",
                self.alpha,
                np.asarray(self.friction_coupling_nm, dtype=float),
                self._get_nmdsigma(),
            )

        return -(self.alpha * self.friction_coupling_nm)[:, np.newaxis] * self.sigma

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
        if self.debug_mf_mode != "none":
            self.beads.p += self.force_mf * pdt

        self._ensure_bath_bound()
        if self.bath is None:
            return
        self.bath.step(pdt)



dproperties(Friction, ["sigma", "gamma", "friction_coupling_nm", "energy_mf", "ediss", "force_mf_nm", "force_mf"])



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

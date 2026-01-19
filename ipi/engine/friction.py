from __future__ import annotations

import numpy as np

from ipi.engine.motion import Motion
from ipi.engine.normalmodes import NormalModes
from ipi.engine.beads import Beads
from ipi.utils.depend import depend_value
from ipi.utils.messages import info, verbosity, warning


class Friction:
    """
    Friction operator for friction-enabled dynamics.

    Objects
    -------
    r : Cartesian coordinates, shape (nbeads, ndof)
    F(r) : coupling coordinates, shape (nbeads, nbath)
    sigma = dF/dr : Jacobian, shape (nbeads, nbath, ndof)

    Markovian (memoryless)
    ----------------------
    General coupling (sigma provided):
        Gamma(r) = sigma^T sigma    (Cartesian friction tensor)
        dp = -(Gamma v) dt + sqrt(2 kBT) sigma^T dW
    Implemented efficiently in coupling space without building Gamma:
        u = sigma v
        dp_drift  = - sigma^T u dt
        dp_noise  = sqrt(2 kBT dt) sigma^T g,  g~N(0,I_nbath)

    Linear coupling (sigma not needed):
        dp_nm = -gamma v_nm dt + sqrt(2 kBT gamma dt) N
      with gamma = friction_static.

    Mean-field (MF)
    ---------------
    Requires alpha^(n) per RP normal mode.
    Policy:
      - If bath_mode == 'markovian': MF is OFF unless alpha_input is provided.
      - If bath_mode == 'non-markovian': alpha can come from alpha_input or spectral_density.

    friction_static
    ---------------
    - General coupling: scales coupling strength via sigma_eff = friction_static * sigma
      (so friction and noise scale consistently).
    - Linear coupling + Markovian: friction_static is interpreted as Markovian gamma.
    """

    # -------------------------
    # Input-configured
    # -------------------------
    use_linear_coupling: bool
    bath_mode: str       # "none" | "markovian" | "non-markovian"
    mf_mode: str         # "none" | "linear" | "reconstruct"

    spectral_density: np.ndarray  # [omega, J(omega)] (required for non-markovian OU fit)
    alpha_input: np.ndarray      # optional [omega_k, alpha]

    friction_static: float

    # OU fitting (non-markovian)
    ou_fit_kind: str
    ou_nterms: int
    ou_tmax: float
    ou_nt: int
    ou_print: bool
    ou_propagator: str           # "exact" | "euler"

    sigma_key: str
    friction_key: str

    # -------------------------
    # Runtime bound
    # -------------------------
    beads: Beads
    nm: NormalModes

    def __init__(
        self,
        use_linear_coupling: bool = False,
        bath_mode: str = "non-markovian",
        mf_mode: str = "reconstruct",
        spectral_density=np.zeros((0, 2), float),
        alpha_input=np.zeros((0, 2), float),
        efric: float = 0.0,
        friction_static: float = 1.0,
        ou_fit_kind: str = "damped_cosine",
        ou_nterms: int = 4,
        ou_tmax: float = 200.0,
        ou_nt: int = 2000,
        ou_print: bool = True,
        ou_propagator: str = "exact",
        sigma_key: str = "sigma",
        friction_key: str = "friction",
    ):
        self.use_linear_coupling = bool(use_linear_coupling)
        self.bath_mode = str(bath_mode)
        self.mf_mode = str(mf_mode)

        self.spectral_density = np.asanyarray(spectral_density, dtype=float)
        self.alpha_input = np.asanyarray(alpha_input, dtype=float)

        self.friction_static = float(friction_static)

        self.ou_fit_kind = str(ou_fit_kind)
        self.ou_nterms = int(ou_nterms)
        self.ou_tmax = float(ou_tmax)
        self.ou_nt = int(ou_nt)
        self.ou_print = bool(ou_print)
        self.ou_propagator = str(ou_propagator)

        self.sigma_key = str(sigma_key)
        self.friction_key = str(friction_key)

        # MF energy only (conservative term)
        self._efric = depend_value(name="efric", value=float(efric))

        # runtime handles
        self.alpha: np.ndarray | None = None
        self.forces = None
        self.ensemble = None
        self.prng = None

        # OU params: (nmodes, nterms, 3) = [c, gamma, omega]
        self._ou_params: np.ndarray | None = None
        # OU auxiliary states: (nmodes, nbath, nterms, 2)
        self._ou_y: np.ndarray | None = None

        # Reconstruction memory for F
        self._F_beads: np.ndarray | None = None  # (nbeads, nbath)
        self._q_prev: np.ndarray | None = None   # (nbeads, ndof)

        self._istep: int = 0

    @property
    def efric(self) -> float:
        return float(self._efric.value)

    @efric.setter
    def efric(self, v: float) -> None:
        self._efric.value = float(v)

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

        if self.mf_mode not in ("none", "linear", "reconstruct"):
            raise ValueError("mf_mode must be one of: 'none', 'linear', 'reconstruct'.")

        if self.ou_propagator not in ("exact", "euler"):
            raise ValueError("ou_propagator must be one of: 'exact', 'euler'.")

        # Markovian MF policy: off unless alpha_input exists
        if self.bath_mode == "markovian" and self.mf_mode != "none":
            if self.alpha_input.size == 0:
                warning(
                    "Friction: bath_mode='markovian' -> MF is disabled unless alpha_input is provided. "
                    "Setting mf_mode='none'.",
                    verbosity.low,
                )
                self.mf_mode = "none"

        # Setup alpha for MF (may be zeros if MF disabled)
        self.alpha = self._setup_alpha()

        # Non-Markovian OU fit requires spectral density
        if self.bath_mode == "non-markovian":
            if self.spectral_density.size == 0:
                raise ValueError("non-markovian requires spectral_density to fit OU embedding (or provide alpha_input + OU params).")
            self._fit_ou_params()

        # init reconstruction memory of positions
        self._q_prev = np.array(self.beads.q, copy=True)

        info(
            "Friction.bind:\n"
            f"  use_linear_coupling = {self.use_linear_coupling}\n"
            f"  bath_mode           = {self.bath_mode}\n"
            f"  mf_mode             = {self.mf_mode}\n"
            f"  friction_static     = {self.friction_static}\n"
            f"  ou_fit_kind         = {self.ou_fit_kind}\n"
            f"  ou_propagator       = {self.ou_propagator}\n"
            f"  sigma_key           = '{self.sigma_key}'\n",
            verbosity.low,
        )

        if (not self.use_linear_coupling) and self.mf_mode == "linear":
            warning(
                "mf_mode='linear' is only meaningful with use_linear_coupling=True. Treating as mf_mode='none'.",
                verbosity.low,
            )
            self.mf_mode = "none"

        if self.bath_mode == "non-markovian" and self.ou_propagator == "euler":
            warning(
                "Friction: ou_propagator='euler' is legacy/debug and can be unstable for stiff gamma*dt >> 1.",
                verbosity.low,
            )

    # ==========================================================================
    # temperature helper
    # ==========================================================================

    def _kbt_rp(self) -> float:
        """kB * (P*T) for ring-polymer effective classical temperature."""
        try:
            from ipi.utils.units import Constants
            kb = Constants.kb
            return kb * float(self.ensemble.temp) * float(self.beads.nbeads)
        except Exception:
            return float(self.ensemble.temp) * float(self.beads.nbeads)

    # ==========================================================================
    # Alpha setup (MF)
    # ==========================================================================

    def _setup_alpha(self) -> np.ndarray:
        wk = np.asarray(self.nm.omegak, dtype=float)
        nmodes = wk.size

        # alpha_input overrides everything
        if self.alpha_input.size > 0:
            if self.alpha_input.ndim != 2 or self.alpha_input.shape[1] != 2:
                raise ValueError("alpha_input must have shape (nmodes,2) [omega_k, alpha].")
            if self.alpha_input.shape[0] != nmodes:
                raise ValueError(f"alpha_input rows ({self.alpha_input.shape[0]}) != nmodes ({nmodes}).")
            if not np.allclose(self.alpha_input[:, 0], wk):
                raise ValueError("alpha_input omega_k does not match current nm.omegak.")
            info("Friction: using alpha from alpha_input table.", verbosity.low)
            return np.asarray(self.alpha_input[:, 1], dtype=float)

        # If MF disabled, keep alpha zeros
        if self.mf_mode == "none":
            return np.zeros(nmodes, dtype=float)

        # If MF requested but no alpha_input:
        # - for non-markovian we can compute from spectral_density
        # - for markovian we should already have disabled MF in bind()
        if self.spectral_density.size > 0:
            omega = np.asarray(self.spectral_density[:, 0], dtype=float)
            J = np.asarray(self.spectral_density[:, 1], dtype=float)
            if omega.size < 2:
                raise ValueError("spectral_density must contain at least two points.")
            if np.any(omega <= 0.0) or np.any(np.diff(omega) <= 0.0):
                raise ValueError("spectral_density omega must be strictly positive and increasing.")
            Lambda = J / omega
            alpha = get_alpha_numeric(Lambda=Lambda, omega=omega, omegak=wk)
            info("Friction: computed alpha^(n) from spectral_density.", verbosity.low)
            return alpha

        # fall back: disable MF (shouldn't happen under the policy, but be safe)
        warning(
            "Friction: MF requested but no alpha_input (and no spectral_density). Disabling MF.",
            verbosity.low,
        )
        self.mf_mode = "none"
        return np.zeros(nmodes, dtype=float)

    # ==========================================================================
    # Extras parsing: sigma = dF/dr
    # ==========================================================================

    def _extras_dict(self) -> dict:
        return getattr(self.forces, "extras", {}) or {}

    @staticmethod
    def _to_array(obj) -> np.ndarray:
        return np.array(obj, dtype=float)

    def _get_sigma(self) -> np.ndarray:
        """
        Return sigma with shape (nbeads, nbath, ndof), interpreted as sigma = dF/dr.

        Accepted inputs:
          - (nbeads, nbath, ndof)
          - (nbath, ndof) broadcast over beads
        """
        extras = self._extras_dict()
        if self.sigma_key not in extras:
            raise KeyError(f"extras['{self.sigma_key}'] not found (sigma=dF/dr required).")

        nbeads = self.beads.nbeads
        ndof = 3 * self.beads.natoms
        raw = self._to_array(extras[self.sigma_key])

        if raw.ndim == 3 and raw.shape[0] == nbeads:
            if raw.shape[2] != ndof:
                raise ValueError(f"sigma last dim must be ndof={ndof}. Got {raw.shape}.")
            return raw

        if raw.ndim == 2:
            if raw.shape[1] != ndof:
                raise ValueError(f"sigma second dim must be ndof={ndof}. Got {raw.shape}.")
            sigma = np.zeros((nbeads, raw.shape[0], raw.shape[1]), float)
            sigma[:] = raw
            return sigma

        raise ValueError(f"sigma must be 2D or 3D numeric array. Got shape {raw.shape}.")

    # ==========================================================================
    # MF reconstruction and forces (general coupling)
    # ==========================================================================

    def _ensure_F_beads(self, nbath: int) -> None:
        nbeads = self.beads.nbeads
        if self._F_beads is None or self._F_beads.shape != (nbeads, nbath):
            self._F_beads = np.zeros((nbeads, nbath), float)

    def _update_reconstructed_F(self, sigma: np.ndarray) -> None:
        """
        Reconstruct coupling coordinate in coupling space:
            F <- F + sigma_eff * dr
        where sigma_eff = friction_static * sigma.
        """
        if self._q_prev is None:
            self._q_prev = np.array(self.beads.q, copy=True)

        q_now = np.array(self.beads.q, copy=False)
        dr = q_now - self._q_prev  # (nbeads, ndof)

        nbath = sigma.shape[1]
        self._ensure_F_beads(nbath)

        sigma_eff = self.friction_static * sigma
        dF = np.einsum("baj,bj->ba", sigma_eff, dr)
        self._F_beads += dF

        self._q_prev = np.array(q_now, copy=True)

    def _mf_forces_reconstruct(self, sigma: np.ndarray) -> np.ndarray:
        """
        MF for general coupling using reconstructed F.

        Energy:
          efric = 1/2 sum_{n,a} alpha[n] * |F_nm[n,a]|^2
        Force:
          g_nm[n,a] = alpha[n] * F_nm[n,a]
          f_beads = - sigma_eff^T g_beads
        """
        if self.alpha is None:
            raise RuntimeError("MF called before bind().")

        nbeads = self.beads.nbeads
        ndof = 3 * self.beads.natoms
        nbath = sigma.shape[1]
        nmodes = int(np.asarray(self.nm.omegak).size)

        self._ensure_F_beads(nbath)

        # Transform F to NM per channel
        F_nm = np.zeros((nmodes, nbath), float)
        for a in range(nbath):
            bead_vec = self._F_beads[:, a].reshape(nbeads, 1)
            nm_vec = self.nm.transform.b2nm(bead_vec)
            F_nm[:, a] = nm_vec[:, 0]

        # g_nm = alpha * F_nm
        g_nm = self.alpha[:, None] * F_nm

        # back to beads
        g_beads = np.zeros((nbeads, nbath), float)
        for a in range(nbath):
            nm_vec = g_nm[:, a].reshape(nmodes, 1)
            bead_vec = self.nm.transform.nm2b(nm_vec)
            g_beads[:, a] = bead_vec[:, 0]

        sigma_eff = self.friction_static * sigma
        f_beads = -np.einsum("baj,ba->bj", sigma_eff, g_beads)

        self.efric = 0.5 * float(np.sum(self.alpha[:, None] * (F_nm * F_nm)))

        return f_beads.reshape(nbeads, ndof)

    def meanfield_forces(self) -> np.ndarray:
        """
        MF force in bead representation.

        - mf_mode='none' : zero
        - linear coupling + mf_mode='linear' : f_nm = -alpha q_nm, transformed to beads
        - general coupling + mf_mode='reconstruct' : reconstructed-F MF (requires sigma)
        """
        nbeads = self.beads.nbeads
        ndof = 3 * self.beads.natoms

        if self.mf_mode == "none":
            self.efric = 0.0
            return np.zeros((nbeads, ndof), float)

        # linear coupling MF
        if self.use_linear_coupling:
            if self.mf_mode != "linear":
                self.efric = 0.0
                return np.zeros((nbeads, ndof), float)

            # NOTE: do NOT scale qnm by friction_static here, because friction_static is used as gamma in Markovian linear.
            qnm = self.nm.qnm
            fnm = -self.alpha[:, None] * qnm
            f_beads = self.nm.transform.nm2b(fnm)

            self.efric = 0.5 * np.einsum("n,nm,nm->", self.alpha, qnm, qnm)
            return f_beads

        # general coupling MF
        if self.mf_mode == "reconstruct":
            sigma = self._get_sigma()
            self._update_reconstructed_F(sigma)
            return self._mf_forces_reconstruct(sigma)

        self.efric = 0.0
        return np.zeros((nbeads, ndof), float)

    # ==========================================================================
    # Markovian bath
    # ==========================================================================

    def _step_markovian(self, pdt: float) -> None:
        """
        Markovian (memoryless) dissipation + noise.

        Linear coupling:
          gamma = friction_static
          dp_nm = -gamma v_nm dt + sqrt(2 kBT gamma dt) N

        General coupling:
          sigma_eff = friction_static * sigma
          dp_drift = -dt sigma_eff^T (sigma_eff v)
          dp_noise = sqrt(2 kBT dt) sigma_eff^T g, g~N(0,I)
        """
        if pdt <= 0.0:
            return

        kbt = self._kbt_rp()

        # Linear coupling: use friction_static as gamma
        if self.use_linear_coupling:
            gamma = float(self.friction_static)
            if gamma < 0.0:
                raise ValueError("friction_static (gamma) must be non-negative for Markovian linear coupling.")

            vnm = self.nm.pnm / self.nm.dynm3  # (nmodes, ndof_nm)
            nmodes, ndof_nm = vnm.shape

            # Drift
            self.nm.pnm += -(gamma * vnm) * pdt

            # Noise
            amp = np.sqrt(2.0 * kbt * gamma * pdt)
            for n in range(nmodes):
                self.nm.pnm[n, :] += amp * self.prng.gvec(ndof_nm)
            return

        # General coupling: requires sigma
        sigma = self._get_sigma()  # (nbeads, nbath, ndof)
        nbeads, nbath, ndof = sigma.shape

        sigma_eff = self.friction_static * sigma

        v = self.beads.p / self.beads.m3  # (nbeads, ndof)

        # u = sigma_eff v (coupling velocity)
        u = np.einsum("baj,bj->ba", sigma_eff, v)  # (nbeads, nbath)

        # drift: dp = - sigma_eff^T u dt
        self.beads.p += -np.einsum("baj,ba->bj", sigma_eff, u) * pdt

        # noise: dp = sqrt(2 kBT dt) sigma_eff^T g
        amp = np.sqrt(2.0 * kbt * pdt)
        for b in range(nbeads):
            g = self.prng.gvec(nbath)  # N(0,1) in coupling space
            self.beads.p[b, :] += amp * (sigma_eff[b, :, :].T @ g)

    # ==========================================================================
    # Non-Markovian OU: exact block update
    # ==========================================================================

    def _ou_block_exact_update(self, y: np.ndarray, u: np.ndarray, g: float, w: float, dt: float, kbt: float) -> None:
        """
        Exact discrete-time update for one damped-cosine OU block:

            d/dt y0 = -g y0 + w y1 + u + noise
            d/dt y1 = -w y0 - g y1     + noise

        Uses exact exp(M dt) + exact forcing integral + exact discrete-time noise covariance:
          cov_increment = kBT(1 - exp(-2 g dt)) I
        """
        if dt <= 0.0:
            return

        ed = np.exp(-g * dt)
        c = np.cos(w * dt)
        s = np.sin(w * dt)

        y0 = y[:, 0].copy()
        y1 = y[:, 1].copy()

        # homogeneous
        y[:, 0] = ed * (c * y0 + s * y1)
        y[:, 1] = ed * (-s * y0 + c * y1)

        # forcing integral: M^{-1}(T - I) [u,0]
        db0 = (ed * c - 1.0) * u
        db1 = (ed * (-s)) * u

        denom = g * g + w * w
        if denom > 0.0:
            y[:, 0] += (-g * db0 - w * db1) / denom
            y[:, 1] += (w * db0 - g * db1) / denom
        else:
            y[:, 0] += u * dt

        # exact discrete-time noise
        if kbt > 0.0 and g > 0.0:
            sig = np.sqrt(kbt * (1.0 - np.exp(-2.0 * g * dt)))
            y[:, 0] += sig * self.prng.gvec(y.shape[0])
            y[:, 1] += sig * self.prng.gvec(y.shape[0])

    # ==========================================================================
    # Non-Markovian OU: fit + step
    # ==========================================================================

    def _fit_ou_params(self) -> None:
        """
        Fit mode-dependent K^(n)(t) derived from Lambda(omega)=J(omega)/omega to:
          exp:           sum_j c exp(-gamma t)
          damped_cosine: sum_j c exp(-gamma t) cos(omega t)

        Stores:
          self._ou_params[n,j] = (c, gamma, omega)
        """
        try:
            import scipy.optimize  # noqa: F401
        except ModuleNotFoundError as e:
            raise ModuleNotFoundError("non-markovian OU fitting requires scipy.") from e

        omega = np.asarray(self.spectral_density[:, 0], float)
        J = np.asarray(self.spectral_density[:, 1], float)
        if omega.size < 2:
            raise ValueError("spectral_density must contain at least two points.")
        if np.any(omega <= 0.0) or np.any(np.diff(omega) <= 0.0):
            raise ValueError("spectral_density omega must be strictly positive and increasing.")
        Lambda = J / omega

        wk = np.asarray(self.nm.omegak, float)
        nmodes = wk.size

        t = np.linspace(0.0, self.ou_tmax, int(self.ou_nt), dtype=float)

        dw = np.diff(omega)
        trap_w = np.zeros_like(omega)
        trap_w[0] = dw[0] / 2.0
        trap_w[-1] = dw[-1] / 2.0
        trap_w[1:-1] = (dw[:-1] + dw[1:]) / 2.0

        pref = 2.0 / np.pi
        Kt = np.zeros((nmodes, t.size), float)

        for n in range(nmodes):
            if n == 0:
                phase = np.outer(omega, t)
                integrand = Lambda[:, None] * np.cos(phase)
            else:
                fac = (omega * omega) / (omega * omega + wk[n] * wk[n])
                wtil = np.sqrt(omega * omega + wk[n] * wk[n])
                phase = np.outer(wtil, t)
                integrand = (Lambda * fac)[:, None] * np.cos(phase)
            Kt[n] = pref * np.sum(integrand * trap_w[:, None], axis=0)

        import scipy.optimize

        self._ou_params = np.zeros((nmodes, self.ou_nterms, 3), float)

        for n in range(nmodes):
            y = Kt[n]

            if self.ou_fit_kind == "exp":
                gam0 = np.logspace(-2, 2, self.ou_nterms)
                c0 = np.maximum(y[0] / self.ou_nterms, 1e-8) * np.ones(self.ou_nterms)
                x0 = np.concatenate([c0, gam0])

                def model(x):
                    c = x[: self.ou_nterms]
                    g = np.abs(x[self.ou_nterms :])
                    return np.sum(c[:, None] * np.exp(-g[:, None] * t[None, :]), axis=0)

                res = scipy.optimize.least_squares(lambda x: model(x) - y, x0, max_nfev=2000)
                c = res.x[: self.ou_nterms]
                g = np.abs(res.x[self.ou_nterms :])
                self._ou_params[n, :, 0] = c
                self._ou_params[n, :, 1] = g
                self._ou_params[n, :, 2] = 0.0

            elif self.ou_fit_kind == "damped_cosine":
                gam0 = np.logspace(-2, 2, self.ou_nterms)
                omg0 = np.linspace(omega[0], omega[-1], self.ou_nterms + 2)[1:-1]
                c0 = np.maximum(y[0] / self.ou_nterms, 1e-8) * np.ones(self.ou_nterms)
                x0 = np.concatenate([c0, gam0, omg0])

                def model(x):
                    c = x[: self.ou_nterms]
                    g = np.abs(x[self.ou_nterms : 2 * self.ou_nterms])
                    om = np.abs(x[2 * self.ou_nterms :])
                    return np.sum(
                        c[:, None] * np.exp(-g[:, None] * t[None, :]) * np.cos(om[:, None] * t[None, :]),
                        axis=0,
                    )

                res = scipy.optimize.least_squares(lambda x: model(x) - y, x0, max_nfev=4000)
                c = res.x[: self.ou_nterms]
                g = np.abs(res.x[self.ou_nterms : 2 * self.ou_nterms])
                om = np.abs(res.x[2 * self.ou_nterms :])
                self._ou_params[n, :, 0] = c
                self._ou_params[n, :, 1] = g
                self._ou_params[n, :, 2] = om
            else:
                raise ValueError("ou_fit_kind must be 'exp' or 'damped_cosine'.")

        if self.ou_print:
            info("Friction: OU fit parameters (per normal mode):", verbosity.low)
            for n in range(nmodes):
                c = self._ou_params[n, :, 0]
                g = self._ou_params[n, :, 1]
                om = self._ou_params[n, :, 2]
                info(
                    f"  mode n={n:3d}  wk={wk[n]:.6g}\n"
                    f"    c     = {np.array2string(c, precision=4)}\n"
                    f"    gamma = {np.array2string(g, precision=4)}\n"
                    f"    omega = {np.array2string(om, precision=4)}",
                    verbosity.low,
                )

        self._ou_y = None

    def _ensure_ou_state(self, nbath: int) -> None:
        nmodes = int(np.asarray(self.nm.omegak).size)
        if self._ou_y is None or self._ou_y.shape[1] != nbath:
            self._ou_y = np.zeros((nmodes, nbath, self.ou_nterms, 2), float)

    def _step_non_markovian(self, pdt: float) -> None:
        if self._ou_params is None:
            raise RuntimeError("non-markovian requested but OU params not initialized.")

        kbt = self._kbt_rp()
        c = self._ou_params[:, :, 0]
        gam = self._ou_params[:, :, 1]
        omg = self._ou_params[:, :, 2]

        nmodes = int(np.asarray(self.nm.omegak).size)

        # Linear coupling: act directly in NM space
        if self.use_linear_coupling:
            nbath = 3 * self.beads.natoms  # NM DOF count
            self._ensure_ou_state(nbath)

            u_nm = self.nm.pnm / self.nm.dynm3  # (nmodes, nbath)
            Phi_nm = np.zeros((nmodes, nbath), float)

            for n in range(nmodes):
                for j in range(self.ou_nterms):
                    g = float(gam[n, j])
                    w = float(omg[n, j])
                    y = self._ou_y[n, :, j, :]

                    if self.ou_propagator == "exact":
                        self._ou_block_exact_update(y=y, u=u_nm[n, :], g=g, w=w, dt=pdt, kbt=kbt)
                    else:
                        # legacy Euler (not recommended)
                        drift0 = -(g * y[:, 0] - w * y[:, 1])
                        drift1 = -(w * y[:, 0] + g * y[:, 1])
                        y[:, 0] += (drift0 + u_nm[n, :]) * pdt
                        y[:, 1] += drift1 * pdt
                        dW0 = self.prng.gvec(nbath) * np.sqrt(pdt)
                        dW1 = self.prng.gvec(nbath) * np.sqrt(pdt)
                        y[:, 0] += np.sqrt(2.0 * kbt * g) * dW0
                        y[:, 1] += np.sqrt(2.0 * kbt * g) * dW1

                    self._ou_y[n, :, j, :] = y
                    Phi_nm[n, :] += c[n, j] * y[:, 0]

            # bath kick in NM space
            self.nm.pnm += -(Phi_nm) * pdt
            return

        # General coupling: coupling space nbath comes from sigma
        sigma = self._get_sigma()
        nbeads = sigma.shape[0]
        nbath = sigma.shape[1]
        self._ensure_ou_state(nbath)

        v_beads = self.beads.p / self.beads.m3
        sigma_eff = self.friction_static * sigma
        u_beads = np.einsum("baj,bj->ba", sigma_eff, v_beads)  # (nbeads, nbath)

        # Transform u to NM per coupling channel
        u_nm = np.zeros((nmodes, nbath), float)
        for a in range(nbath):
            bead_vec = u_beads[:, a].reshape(nbeads, 1)
            nm_vec = self.nm.transform.b2nm(bead_vec)
            u_nm[:, a] = nm_vec[:, 0]

        Phi_nm = np.zeros((nmodes, nbath), float)
        for n in range(nmodes):
            for j in range(self.ou_nterms):
                g = float(gam[n, j])
                w = float(omg[n, j])
                y = self._ou_y[n, :, j, :]

                if self.ou_propagator == "exact":
                    self._ou_block_exact_update(y=y, u=u_nm[n, :], g=g, w=w, dt=pdt, kbt=kbt)
                else:
                    drift0 = -(g * y[:, 0] - w * y[:, 1])
                    drift1 = -(w * y[:, 0] + g * y[:, 1])
                    y[:, 0] += (drift0 + u_nm[n, :]) * pdt
                    y[:, 1] += drift1 * pdt
                    dW0 = self.prng.gvec(nbath) * np.sqrt(pdt)
                    dW1 = self.prng.gvec(nbath) * np.sqrt(pdt)
                    y[:, 0] += np.sqrt(2.0 * kbt * g) * dW0
                    y[:, 1] += np.sqrt(2.0 * kbt * g) * dW1

                self._ou_y[n, :, j, :] = y
                Phi_nm[n, :] += c[n, j] * y[:, 0]

        # Transform Phi back to beads per coupling channel
        Phi_beads = np.zeros((nbeads, nbath), float)
        for a in range(nbath):
            nm_vec = Phi_nm[:, a].reshape(nmodes, 1)
            bead_vec = self.nm.transform.nm2b(nm_vec)
            Phi_beads[:, a] = bead_vec[:, 0]

        # backproject to Cartesian
        self.beads.p += -np.einsum("baj,ba->bj", sigma_eff, Phi_beads) * pdt

    # ==========================================================================
    # main step
    # ==========================================================================

    def step(self, pdt: float) -> None:
        """
        Apply friction operator over time interval pdt:
          - MF kick (if enabled)
          - bath kick (markovian/non-markovian)
        """
        self._istep += 1

        # MF
        f_mf = self.meanfield_forces()
        if np.any(f_mf):
            self.beads.p += f_mf * pdt

        # Bath
        if self.bath_mode == "none":
            return
        if self.bath_mode == "markovian":
            self._step_markovian(pdt)
            return
        if self.bath_mode == "non-markovian":
            self._step_non_markovian(pdt)
            return

        raise RuntimeError("Unreachable bath_mode branch.")


def get_alpha_numeric(Lambda: np.ndarray, omega: np.ndarray, omegak: np.ndarray) -> np.ndarray:
    """Numerically compute alpha^(n) from Lambda(omega)."""
    try:
        from scipy.interpolate import CubicSpline
    except ModuleNotFoundError as e:
        raise ModuleNotFoundError(
            "Friction: scipy is required to compute alpha from spectral_density. "
            "Install scipy or provide alpha_input explicitly."
        ) from e

    omega = np.asarray(omega, dtype=float)
    Lambda = np.asarray(Lambda, dtype=float)
    omegak = np.asarray(omegak, dtype=float)

    alpha = np.zeros(omegak.shape, dtype=float)
    for i, wk in enumerate(omegak):
        f = CubicSpline(omega, Lambda * (wk ** 2) / (omega ** 2 + wk ** 2))
        alpha[i] = (2.0 / np.pi) * f.integrate(0.0, omega[-1])
        info(f"Friction: wk={wk} alpha={alpha[i]}", verbosity.high)

    return alpha

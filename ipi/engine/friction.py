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

    Contract with i-PI extras (IMPORTANT)
    ------------------------------------
    i-PI calls the force provider once per bead and collects each bead's `extra`
    dict. In ForceComponent.extra_gather() it builds a dict-of-lists:
        forces.extras[key] == [payload_bead0, payload_bead1, ...]  (length nbeads)

    This class therefore assumes sigma/friction is provided PER BEAD, and will
    raise if the extras do not have length == nbeads.

    Objects
    -------
    r : Cartesian coordinates, shape (nbeads, ndof)
    sigma = dF/dr : Jacobian, shape (nbeads, nbath, ndof)

    Markovian (memoryless)
    ----------------------
    General coupling:
        u = sigma_eff v
        dp_drift = -sigma_eff^T u dt
        dp_noise = sqrt(2 kBT dt) sigma_eff^T g, g~N(0,I_nbath)

    Linear coupling:
        dp_nm = -gamma v_nm dt + sqrt(2 kBT gamma dt) N
        with gamma = friction_static.

    Mean-field (MF)
    ---------------
    For general coupling with reconstruction:
        F <- F + sigma_eff dr
        E = 1/2 sum_{n,a} alpha[n] |F_nm[n,a]|^2
        f = - sigma_eff^T g_beads, g_nm = alpha F_nm

    friction_static
    ---------------
    - General coupling: scales sigma as sigma_eff = friction_static * sigma
    - Linear coupling + Markovian: friction_static is gamma

    """

    # -------------------------
    # Input-configured
    # -------------------------
    use_linear_coupling: bool
    bath_mode: str       # "none" | "markovian" | "non-markovian"
    mf_mode: str         # "none" | "linear" | "reconstruct"

    spectral_density: np.ndarray  # [omega, J(omega)] for non-markovian OU fit
    alpha_input: np.ndarray      # optional [omega_k, alpha]

    friction_static: float

    # OU fitting (non-markovian)
    ou_fit_kind: str
    ou_nterms: int
    ou_tmax: float
    ou_nt: int
    ou_print: bool
    ou_propagator: str           # "exact" | "euler"

    # Extras parsing
    sigma_key: str
    friction_key: str

    # Optional ACE sigma selector (if your payload is {"sigma": {"equ": {"1": ...}}})
    sigma_rep: str | None
    sigma_index: str | int | None

    # Optional: if driver returns reduced matrix for a subset of atoms
    friction_atoms: np.ndarray | list | None  # 1-based atom indices

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
        sigma_rep: str | None = None,
        sigma_index: str | int | None = None,
        friction_atoms = np.zeros(0, dtype=int),
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

        self.sigma_rep = sigma_rep
        self.sigma_index = sigma_index

        self.friction_atoms = None if friction_atoms is None else np.asarray(friction_atoms)

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
                    "Friction: bath_mode='markovian' -> MF disabled unless alpha_input is provided. "
                    "Setting mf_mode='none'.",
                    verbosity.low,
                )
                self.mf_mode = "none"

        # Setup alpha for MF (may be zeros if MF disabled)
        self.alpha = self._setup_alpha()

        # Non-Markovian OU fit requires spectral density
        if self.bath_mode == "non-markovian":
            if self.spectral_density.size == 0:
                raise ValueError(
                    "non-markovian requires spectral_density to fit OU embedding "
                    "(or provide alpha_input + pre-set OU params)."
                )
            # self._fit_ou_params()

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
            f"  sigma_key           = '{self.sigma_key}'\n"
            f"  friction_key        = '{self.friction_key}'\n",
            verbosity.low,
        )

        if (not self.use_linear_coupling) and self.mf_mode == "linear":
            warning(
                "mf_mode='linear' is only meaningful with use_linear_coupling=True. Setting mf_mode='none'.",
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
        # - for non-markovian compute from spectral_density
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

        warning(
            "Friction: MF requested but no alpha_input (and no spectral_density). Disabling MF.",
            verbosity.low,
        )
        self.mf_mode = "none"
        return np.zeros(nmodes, dtype=float)

    # ==========================================================================
    # Extras parsing: sigma = dF/dr (STRICT PER-BEAD)
    # ==========================================================================



    @staticmethod
    def _to_array(obj) -> np.ndarray:
        return np.asarray(obj, dtype=float)



    def _embed_if_needed(self, M2d: np.ndarray, natoms: int) -> np.ndarray:
        """
        If M2d is reduced (3N_fric × 3N_fric), embed into full (3Nat × 3Nat)
        using self.friction_atoms (1-based atom indices).
        """
        ndof = 3 * natoms
        if M2d.shape == (ndof, ndof):
            return M2d

        if self.friction_atoms is None or np.size(self.friction_atoms) == 0:
            raise ValueError(
                f"Received matrix shape {M2d.shape} but full ndof is {ndof}. "
                "Provide friction_atoms (1-based atom indices) to embed, or send full 3Nat×3Nat."
            )

        fa = [int(a) - 1 for a in np.asarray(self.friction_atoms, dtype=int).reshape(-1)]
        nfa = len(fa)
        expected = 3 * nfa

        if M2d.shape != (expected, expected):
            raise ValueError(
                f"Reduced matrix shape {M2d.shape} does not match expected {(expected, expected)} "
                f"from friction_atoms (len={nfa})."
            )

        if any(a < 0 or a >= natoms for a in fa):
            raise ValueError(f"friction_atoms out of range for natoms={natoms}: {self.friction_atoms}")

        out = np.zeros((ndof, ndof), dtype=float)
        for bi, ai in enumerate(fa):
            ri = slice(3 * ai, 3 * ai + 3)
            sri = slice(3 * bi, 3 * bi + 3)
            for bj, aj in enumerate(fa):
                rj = slice(3 * aj, 3 * aj + 3)
                srj = slice(3 * bj, 3 * bj + 3)
                out[ri, rj] = M2d[sri, srj]
        return out

    def _get_sigma(self) -> np.ndarray:
        """ 

        Reads i-PI combined extras:
            extras[key] is expected to be a list/tuple of length nbeads,
            each entry being the per-bead payload (a dict).

        Returns:
            sigma : (nbeads, nbath, ndof)

        Convention:
          - We interpret the received 2D matrix for each bead as sigma[b] with shape (nbath, ndof).
          - Common case nbath == ndof, so sigma[b] is square.
        """
        info("Inside get_sigma", verbosity.low) 

        nbeads = int(self.beads.nbeads)
        natoms = int(self.beads.natoms)
        ndof = 3 * natoms

        sigma = self.forces.extras.get(self.sigma_key)
        if sigma is None:
            raise KeyError(
                f"Did not find 'sigma' among the force extras = {self.forces.extras}"
            )
        else:
            info(str(sigma), verbosity.low)
            return np.asarray(sigma)

    # ==========================================================================
    # MF reconstruction and forces (general coupling)
    # ==========================================================================

    def _ensure_F_beads(self, nbath: int) -> None:
        nbeads = int(self.beads.nbeads)
        if self._F_beads is None or self._F_beads.shape != (nbeads, nbath):
            self._F_beads = np.zeros((nbeads, nbath), float)

    def _update_reconstructed_F(self, sigma: np.ndarray) -> None:
        """
        Reconstruct coupling coordinate:
            F <- F + sigma_eff dr
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

        E = 1/2 sum_{n,a} alpha[n] |F_nm[n,a]|^2
        g_nm[n,a] = alpha[n] * F_nm[n,a]
        f_beads = - sigma_eff^T g_beads
        """
        if self.alpha is None:
            raise RuntimeError("MF called before bind().")

        nbeads = int(self.beads.nbeads)
        ndof = 3 * int(self.beads.natoms)
        nbath = sigma.shape[1]
        nmodes = int(np.asarray(self.nm.omegak).size)

        self._ensure_F_beads(nbath)

        # Transform F to NM per channel
        F_nm = np.zeros((nmodes, nbath), float)
        for a in range(nbath):
            bead_vec = self._F_beads[:, a].reshape(nbeads, 1)
            nm_vec = self.nm.transform.b2nm(bead_vec)
            F_nm[:, a] = nm_vec[:, 0]

        g_nm = self.alpha[:, None] * F_nm

        # Back to beads
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
        nbeads = int(self.beads.nbeads)
        ndof = 3 * int(self.beads.natoms)

        if self.mf_mode == "none":
            self.efric = 0.0
            return np.zeros((nbeads, ndof), float)

        # linear coupling MF
        if self.use_linear_coupling:
            if self.mf_mode != "linear":
                self.efric = 0.0
                return np.zeros((nbeads, ndof), float)

            qnm = self.nm.qnm
            fnm = -self.alpha[:, None] * qnm
            f_beads = self.nm.transform.nm2b(fnm)
            self.efric = 0.5 * np.einsum("n,nm,nm->", self.alpha, qnm, qnm)
            return f_beads

        # general coupling MF
        if self.mf_mode == "reconstruct":
            sigma = self._get_sigma()  # (nbeads, nbath, ndof)
            self._update_reconstructed_F(sigma)
            return self._mf_forces_reconstruct(sigma)

        self.efric = 0.0
        return np.zeros((nbeads, ndof), float)

    # ==========================================================================
    # Markovian bath
    # ==========================================================================

    def _step_markovian(self, pdt: float) -> None:
        """
        Markovian dissipation + noise.

        Linear coupling:
          gamma = friction_static
          dp_nm = -gamma v_nm dt + sqrt(2 kBT gamma dt) N

        General coupling:
          sigma_eff = friction_static * sigma
          dp_drift = -dt sigma_eff^T (sigma_eff v)
          dp_noise = sqrt(2 kBT dt) sigma_eff^T g, g~N(0,I)
        """

        info("Inside Markovian step:", verbosity.low)
        
        if pdt <= 0.0:
            return

        kbt = self._kbt_rp()

        # Linear coupling: friction_static is gamma
        if self.use_linear_coupling:
            gamma = float(self.friction_static)
            if gamma < 0.0:
                raise ValueError("friction_static (gamma) must be non-negative for Markovian linear coupling.")

            vnm = self.nm.pnm / self.nm.dynm3
            nmodes, ndof_nm = vnm.shape

            self.nm.pnm += -(gamma * vnm) * pdt

            amp = np.sqrt(2.0 * kbt * gamma * pdt)
            for n in range(nmodes):
                self.nm.pnm[n, :] += amp * self.prng.gvec(ndof_nm)
            return

        info("Markov: before get sigma", verbosity.low)
        # General coupling: requires sigma per bead
        sigma = self._get_sigma()  # (nbeads, nbath, ndof)
        info("Markov: after get sigma", verbosity.low)
        nbeads, nbath, ndof = sigma.shape

        sigma_eff = self.friction_static * sigma
        v = self.beads.p / self.beads.m3

        u = np.einsum("baj,bj->ba", sigma_eff, v)
        self.beads.p += -np.einsum("baj,ba->bj", sigma_eff, u) * pdt

        amp = np.sqrt(2.0 * kbt * pdt)
        for b in range(nbeads):
            g = self.prng.gvec(nbath)
            self.beads.p[b, :] += amp * (sigma_eff[b, :, :].T @ g)

    # ==========================================================================
    # Non-Markovian OU: exact block update
    # ==========================================================================

    def _ou_block_exact_update(
        self,
        y: np.ndarray,
        u: np.ndarray,
        g: float,
        w: float,
        dt: float,
        kbt: float,
    ) -> None:
        """
        Exact discrete-time update for one damped-cosine OU block:

            d/dt y0 = -g y0 + w y1 + u + noise
            d/dt y1 = -w y0 - g y1     + noise

        Uses exact exp(M dt) + forcing integral + exact discrete-time noise covariance:
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

    # Actually will be outsourced to George's code
    # =========================================================================

    def _step_non_markovian(self, pdt: float) -> None:
        if self._ou_params is None:
            raise RuntimeError("non-markovian requested but OU params not initialized.") 
        
        return

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
        f = CubicSpline(omega, Lambda * (wk**2) / (omega**2 + wk**2))
        alpha[i] = (2.0 / np.pi) * f.integrate(0.0, omega[-1])
        info(f"Friction: wk={wk} alpha={alpha[i]}", verbosity.high)

    return alpha

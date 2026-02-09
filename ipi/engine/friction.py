from __future__ import annotations

import numpy as np

from ipi.engine.motion import Motion
from ipi.engine.normalmodes import NormalModes
from ipi.engine.beads import Beads

from ipi.utils.depend import depend_value,  dproperties
from ipi.utils.messages import info, verbosity, warning


class Friction:
    """
    Friction operator for friction-enabled dynamics.
    """

    # -------------------------
    # Input-configured
    # -------------------------
    variable_friction: bool
    bath_mode: str       # "none" | "markovian" | "non-markovian"
    mf_mode: str         # "none" | "linear" | "reconstruct"

    Lambda: np.ndarray  # [omega, J(omega)] for non-markovian OU fit
    alpha_input: np.ndarray      # optional [omega_k, alpha]

    friction_static: float

    # OU fitting (non-markovian)
    ou_nterms: int
    ou_tmax: float
    ou_nt: int
    ou_print: bool
    ou_propagator: str           # "exact" | "euler"

    # Extras parsing
    sigma_key: str
    friction_key: str

    # Optional: if driver returns reduced matrix for a subset of atoms
    friction_atoms: np.ndarray | list | None  # 1-based atom indices

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
        bath_mode: str = "non-markovian",
        mf_mode: str = "reconstruct",
        Lambda=np.zeros((0, 2), float),
        alpha_input=np.zeros((0, 2), float),
        friction_static: float = 1.0,
        ou_nterms: int = 4,
        ou_tmax: float = 200.0,
        ou_nt: int = 2000,
        ou_print: bool = True,
        ou_propagator: str = "exact",
        sigma_key: str = "sigma",
        friction_key: str = "friction",
        friction_atoms = np.zeros(0, dtype=int),
    ):
        """Initialises the friction object.
        Args:
            Lambda: Cosine transform of the time-dependent factor in the friction kernel,
                divided by frequncy. Supplied as a 2d array of two columns containing frequency and
                spectral density, respectively.
                Defaults to np.zeros(0, float).
            alpha_input: Normal-mode coefficients in expression for the frictional mean-field
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
        self.mf_mode = str(mf_mode)

        # Kernel shape
        self.Lambda = np.asanyarray(Lambda, dtype=float).copy()
        self.alpha_input = np.asanyarray(alpha_input, dtype=float).copy()

        self._sigma = depend_value(name="sigma", func=self._get_sigma)
    
    #     # Friction coupling: F(q), such that Σ{i,α} = ∂F(q) / ∂q{i,α}
        self._friction_coupling_nm = depend_value(
            name="friction_coupling_nm",
            func=self.get_friction_coupling_nm,
            dependencies=[self._sigma],
        )
        # Frictional mean-field force
        self._fmf_nm = depend_value(
            name="fmf_nm",
            func=self.get_fmf_nm,
            dependencies=[self._friction_coupling_nm],
        )
        self._fmf = depend_value(
            name="fmf", func=self.get_fmf, dependencies=[self._fmf_nm]
        )

        #Conserved mean-field potential
        self._emf = depend_value(name="emf", value=0.0)

        self.friction_static = float(friction_static)

        self.ou_nterms = int(ou_nterms)
        self.ou_tmax = float(ou_tmax)
        self.ou_nt = int(ou_nt)
        self.ou_print = bool(ou_print)
        self.ou_propagator = str(ou_propagator)

        self.sigma_key = str(sigma_key)
        self.friction_key = str(friction_key)

        self.friction_atoms = None if friction_atoms is None else np.asarray(friction_atoms)

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



        # Non-Markovian OU fit requires spectral density
        if self.bath_mode == "non-markovian":
            if self.Lambda.size == 0:
                raise ValueError(
                    "non-markovian requires Lambda to fit OU embedding "
                    "(or provide alpha_input + pre-set OU params)."
                )




            # self._fit_ou_params()

        # Setup alpha for MF (may be zeros if MF disabled)
        self.alpha = self._setup_alpha()

        # init reconstruction memory of positions
        self._q_prev = np.array(self.beads.q, copy=True)

        info(
            "Friction.bind:\n"
            f"  variable_friction = {self.variable_friction}\n"
            f"  bath_mode           = {self.bath_mode}\n"
            f"  mf_mode             = {self.mf_mode}\n"
            f"  friction_static     = {self.friction_static}\n"
            f"  ou_propagator       = {self.ou_propagator}\n"
            f"  sigma_key           = '{self.sigma_key}'\n"
            f"  friction_key        = '{self.friction_key}'\n",
            verbosity.low,
        )

        if (self.variable_friction) and self.mf_mode == "linear":
            warning(
                "mf_mode='linear' is only meaningful with variable_friction=False. Setting mf_mode='none'.",
                verbosity.low,
            )
            self.mf_mode = "none"

        if self.bath_mode == "non-markovian" and self.ou_propagator == "euler":
            warning(
                "Friction: ou_propagator='euler' is legacy/debug and can be unstable for stiff gamma*dt >> 1.",
                verbosity.low,
            )

        # Dependencies
        self._sigma.add_dependency(self.forces._extras)
        self._friction_coupling_nm.add_dependency(self.beads._q)
        self._emf.add_dependency(self._friction_coupling_nm)
        self._emf._func = self.get_emf

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

        warning(
            "Friction: MF requested but no alpha_input (and no Lambda). Disabling MF.",
            verbosity.low,
        )
        self.mf_mode = "none"
        return np.zeros(nmodes, dtype=float)

    # ==========================================================================
    # Parse Sigma
    # ==========================================================================

    def _embed_if_needed(self, M2d: np.ndarray, natoms: int) -> np.ndarray:
        # Placehold for embedded partial sigma delivered by model into full n_atoms version
        return

    # AKA get_diffusion_coefficient
    def _get_sigma(self) -> np.ndarray:
        """ 

        Reads i-PI combined extras:
            extras[key] is expected to be a 
            each entry being the per-bead payload (a dict).

        Returns:
            sigma : (shape)

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

        # sigma_eff = self.friction_static * sigma
        f_beads = -np.einsum("baj,ba->bj", sigma, g_beads)

        # self.emf = 0.5 * float(np.sum(self.alpha[:, None] * (F_nm * F_nm)))
        return f_beads.reshape(nbeads, ndof)

    # def meanfield_forces(self) -> np.ndarray:
    #     """
    #     MF force in bead representation.

    #     - mf_mode='none' : zero
    #     - variable_friction = false i.e linear coupling + mf_mode='linear' : f_nm = -alpha q_nm, transformed to beads
    #     - variable_friction true i.e separable coupling + mf_mode='reconstruct' : reconstructed-F MF (requires sigma)
    #     """
    #     nbeads = int(self.beads.nbeads)
    #     ndof = 3 * int(self.beads.natoms)

    #     if self.mf_mode == "none":
    #         self.emf = 0.0
    #         return np.zeros((nbeads, ndof), float)

    #     # linear coupling MF
    #     if not self.variable_friction:

    #         if self.mf_mode == "reconstruct":
    #             # reconstruct requires sigma from extras. We probably dont ask for that in the linear case. Might be a good validation of reconstruct method in future.
    #             raise ValueError("if variable_friction is false, mf_mode must be linear or none. reconstruct is currently not supported")


    #         if self.mf_mode != "linear":
    #             # either none or reconstruct
    #             self.emf = 0.0
    #             return np.zeros((nbeads, ndof), float)

    #         # analytical mean-field force when the coupling is linear
    #         qnm = self.nm.qnm
    #         fnm = -self.alpha[:, None] * qnm
    #         f_beads = self.nm.transform.nm2b(fnm)
    #         self.emf = 0.5 * np.einsum("n,nm,nm->", self.alpha, qnm, qnm)
    #         return f_beads

    #     # general coupling MF
    #     if self.mf_mode == "reconstruct":
    #         sigma = self._get_sigma()  # (nbeads, nbath, ndof)
    #         self._update_reconstructed_F(sigma)
    #         return self._mf_forces_reconstruct(sigma)

    #     self.emf = 0.0
    #     return np.zeros((nbeads, ndof), float)
    

    def get_friction_coupling_nm(self):
        """Compute the friction coupling for each normal-mode index"""
        if self.variable_friction:
            raise NotImplementedError(
                "The calculation of friction coupling for position-dependent diffusion coefficients is not implemented."
            )
        else:
            # Here we assume that the interaction potential, F(Q) in https://doi.org/10.1103/PhysRevLett.134.226201,
            # is of the form F(q) = SUM[ c{i,α} q{i,α}, {{i,0,n_atom-1}, {α,0,2}} ] where α indexes Cartesian components
            # The diffusion coefficients for bead index l returned by the driver are expected to be packed as
            # Σ{i,α} = ∂F(q) / ∂q{i,α} = diffusion_coeff[l, 3*i+α].
            return np.sum(self.sigma * self.nm.qnm, axis=-1)

    def get_emf(self):
        """Compute the frictional potential of mean field, Eq. (S19) of https://doi.org/10.1103/PhysRevLett.134.226201"""

        if self.mf_mode == "none":
            return 0.0
        
        return np.sum(self.alpha * self.friction_coupling_nm**2) / 2
        
    def get_fmf_nm(self):
        """Negative derivative of the frictional potential of mean field with respect to normal modes"""
        return -(self.alpha * self.friction_coupling_nm)[:, np.newaxis] * self.sigma

    def get_fmf(self):
        """Negative derivative of the frictional potential of mean field with respect to bead positions"""
        return self.nm.transform.nm2b(self.fmf_nm)

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

        kbt = self._kbt()

        # Linear coupling: friction does not depend on instaneous atomic configuration.
        # friction_static is gamma
        if not self.variable_friction:
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
        # variable friction / seperable coupling: requires sigma per bead
        sigma = self._get_sigma()  # (nbeads, nbath, ndof)
        info("Markov: after get sigma", verbosity.low)
        nbeads, nbath, ndof = sigma.shape
        info("nbeads from sigma"+str(nbeads), verbosity.low)

        v = self.beads.p / self.beads.m3

        u = np.einsum("baj,bj->ba", sigma, v)
        self.beads.p += -np.einsum("baj,ba->bj", sigma, u) * pdt

        amp = np.sqrt(2.0 * kbt * pdt)
        for b in range(nbeads):
            g = self.prng.gvec(nbath)
            self.beads.p[b, :] += amp * (sigma[b, :, :].T @ g)

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
    # Non-Markovian OU: step

    # Fit will be outsourced to George's code
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
          - bath kick (markovian/non-markovian) - dissipative and random forces
        """

        # MF
        if self.mf_mode != "none":
            self.beads.p += self.fmf * pdt

        # Bath
        if self.bath_mode == "none":
            return
        if self.bath_mode == "markovian":
            self._step_markovian(pdt)
            return
        if self.bath_mode == "non-markovian":
            self._step_non_markovian(pdt)
            return

        raise RuntimeError("bath_mode must be one of none, markovian or non-markovian")



dproperties(Friction, ["sigma", "friction_coupling_nm", "emf", "fmf_nm", "fmf"])



def get_alpha_numeric(Lambda: np.ndarray, omega: np.ndarray, omegak: np.ndarray) -> np.ndarray:
    """Numerically compute alpha^(n) from Lambda(omega)."""
    try:
        from scipy.interpolate import CubicSpline
    except ModuleNotFoundError as e:
        raise ModuleNotFoundError(
            "Friction: scipy is required to compute alpha from Lambda. "
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

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

    sigma_static: float

    # OU fitting (non-markovian)
    ou_nterms: int
    ou_tmax: float
    ou_nt: int
    ou_print: bool
    ou_propagator: str           # "exact" | "euler"

    # Extras parsing
    sigma_key: str
    friction_key: str

    # Optional: subset of atoms (0-based) that participate in electronic friction
    friction_atoms: np.ndarray | list | None

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
        sigma_static: float = 1.0,
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

        self.sigma_static = float(sigma_static)

        self.ou_nterms = int(ou_nterms)
        self.ou_tmax = float(ou_tmax)
        self.ou_nt = int(ou_nt)
        self.ou_print = bool(ou_print)
        self.ou_propagator = str(ou_propagator)

        self.sigma_key = str(sigma_key)
        self.friction_key = str(friction_key)

        if friction_atoms is None:
            self.friction_atoms = None
        else:
            f_atoms = np.asarray(friction_atoms, dtype=int).flatten()
            # Empty means "not specified" -> all atoms participate (resolved at bind time).
            self.friction_atoms = None if f_atoms.size == 0 else f_atoms
        self._friction_atoms_idx: np.ndarray | None = None
        self._friction_dof_idx: np.ndarray | None = None

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

         # bookkeeping: cumulative energy exchange with bath (Markovian)
        self._ediss = depend_value(name="ediss", value=0.0)     # positive = system -> bath via friction
        self._erand = depend_value(name="erand", value=0.0)  # positive = bath -> system via random force



    # ==========================================================================
    # bind
    # ==========================================================================

    def bind(self, motion: Motion) -> None:
        self.beads = motion.beads
        self.nm = motion.nm
        self.ensemble = motion.ensemble
        self.forces = motion.ensemble.forces
        self.prng = motion.prng
        self._setup_friction_atoms()

        if self.bath_mode not in ("none", "markovian", "non-markovian"):
            raise ValueError("bath_mode must be one of: 'none', 'markovian', 'non-markovian'.")

        if self.mf_mode not in ("none", "linear", "reconstruct"):
            raise ValueError("mf_mode must be one of: 'none', 'linear', 'reconstruct'.")

        if self.ou_propagator not in ("exact", "euler"):
            raise ValueError("ou_propagator must be one of: 'exact', 'euler'.")


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
            f"  sigma_static     = {self.sigma_static}\n"
            f"  ou_propagator       = {self.ou_propagator}\n"
            f"  sigma_key           = '{self.sigma_key}'\n"
            f"  friction_key        = '{self.friction_key}'\n",
            verbosity.low,
        )

        if (self.variable_friction) and self.mf_mode == "linear":
            raise ValueError(
                "mf_mode='linear' is only meaningful with variable_friction=False.")


        # Dependencies
        self._sigma.add_dependency(self.forces._extras)
        self._friction_coupling_nm.add_dependency(self.beads._q)
        self._emf.add_dependency(self._friction_coupling_nm)
        self._emf._func = self.get_emf

    def _setup_friction_atoms(self) -> None:
        """Resolves and validates the 0-based list of atoms affected by friction."""
        natoms = int(self.beads.natoms)
        if self.friction_atoms is None:
            atoms = np.arange(natoms, dtype=int)
        else:
            atoms = np.asarray(self.friction_atoms, dtype=int).flatten()

        if np.any(atoms < 0) or np.any(atoms >= natoms):
            raise ValueError(
                f"friction_atoms must be 0-based indices in [0, {natoms - 1}], got {atoms}."
            )
        if np.unique(atoms).size != atoms.size:
            raise ValueError(f"friction_atoms contains duplicate indices: {atoms}")

        self._friction_atoms_idx = atoms
        self._friction_dof_idx = np.concatenate(
            [np.arange(3 * a, 3 * a + 3, dtype=int) for a in atoms]
        )

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

        raise ValueError(
            "Friction: MF requested but no alpha_input (and no Lambda)."
        )

    # ==========================================================================
    # Parse Sigma
    # ==========================================================================

    def _embed_if_needed(self, sigma3d: np.ndarray, ndof: int) -> np.ndarray:
        """Embeds reduced sigma payload into full Cartesian space with zero padding."""
        ndof_reduced = int(sigma3d.shape[2])
        if ndof_reduced == ndof:
            return sigma3d

        if self._friction_dof_idx is None:
            raise RuntimeError("Friction atoms have not been initialized. Call bind() first.")

        if ndof_reduced != self._friction_dof_idx.size:
            raise ValueError(
                f"Reduced {self.sigma_key} ndof={ndof_reduced} does not match "
                f"3*len(friction_atoms)={self._friction_dof_idx.size}."
            )

        sigma_full = np.zeros(
            (sigma3d.shape[0], sigma3d.shape[1], ndof), dtype=sigma3d.dtype
        )
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
        info("Inside get_sigma", verbosity.low) 


        if (not self.variable_friction):
            return float(self.sigma_static)

        nbeads = int(self.beads.nbeads)
        ndof = 3 * int(self.beads.natoms)

        sigma = self.forces.extras.get(self.sigma_key)


        if sigma is None:
            raise KeyError(
                f"Did not find 'sigma' among the force extras = {self.forces.extras}"
            )
        else:
            info(str(sigma), verbosity.low)
            sigma = np.asarray(sigma,dtype=float)

        
        if sigma.ndim != 3:
            raise ValueError(f"{self.sigma_key} must have ndim=3 (nbeads, nbath, ndof). Got shape {sigma.shape}.")
        if sigma.shape[0] != nbeads:
            raise ValueError(
                f"{self.sigma_key} shape {sigma.shape} incompatible with nbeads={nbeads}."
            )
        sigma = self._embed_if_needed(sigma, ndof=ndof)
        return sigma
    

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
            if self._friction_dof_idx is None:
                raise RuntimeError("Friction atoms have not been initialized. Call bind() first.")
            q_active = self.nm.qnm[:, self._friction_dof_idx]
            return np.sum(self.sigma * q_active, axis=-1)

    def get_emf(self):
        """Compute the frictional potential of mean field, Eq. (S19) of https://doi.org/10.1103/PhysRevLett.134.226201"""

        if self.mf_mode == "none":
            return 0.0

        info("alpha "+str(self.alpha), verbosity.low)
        info("coupling "+str(self.friction_coupling_nm**2))
        info("EMF "+str(np.sum(self.alpha * self.friction_coupling_nm**2)), verbosity.low) 
        return np.sum(self.alpha * self.friction_coupling_nm**2) / 2
        
    def get_fmf_nm(self):
        """Negative derivative of the frictional potential of mean field with respect to normal modes"""
        if self.variable_friction:
            return -(self.alpha * self.friction_coupling_nm)[:, np.newaxis] * self.sigma

        if self._friction_dof_idx is None:
            raise RuntimeError("Friction atoms have not been initialized. Call bind() first.")

        fmf_nm = np.zeros_like(self.nm.qnm)
        fmf_nm[:, self._friction_dof_idx] = (
            -(self.alpha * self.friction_coupling_nm)[:, np.newaxis] * self.sigma
        )
        return fmf_nm

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
          gamma = sigma_static * sigma_static
          dp_nm = -gamma v_nm dt + sqrt(2 kBT gamma dt) N

        Variable friction:
          sigma is position - dependent
          dp_drift = -dt sigma^T (sigma v)
          dp_noise = sqrt(2 kBT dt) sigma^T g, g~N(0,I)
        """
        from ipi.utils.mathtools import matrix_exp, root_herm

        info("Inside Markovian step:", verbosity.low)
        
        if pdt <= 0.0:
            return

        kbt = self._kbt_rp()

        # Linear coupling: friction does not depend on instaneous atomic configuration.
        # Exact OU update for: dp = -gamma v dt + sqrt(2 kBT gamma) dW,  v = p/m
        if not self.variable_friction:
            sigma = float(self.sigma_static)
            gamma = sigma * sigma
            if gamma < 0.0:
                raise ValueError("gamma must be non-negative for Markovian linear coupling.")
            if gamma == 0.0:
                return  # no coupling

            # masses and momenta in NM space
            m = self.nm.dynm3.copy()   # (nmodes, ndof_nm)
            p = self.nm.pnm.copy()     # (nmodes, ndof_nm)
            sm = np.sqrt(m)
            if self._friction_dof_idx is None:
                raise RuntimeError("Friction atoms have not been initialized. Call bind() first.")

            info("before: " + str(self.beads.p), verbosity.low)

            # --- ThermoLangevin coefficients but elementwise (because tau = m/gamma) ---
            # tau_ij = m_ij / gamma  =>  T_ij = exp(-pdt/tau_ij) = exp(-(gamma/m_ij)*pdt)
            gamma_nm = np.zeros_like(m)
            gamma_nm[:, self._friction_dof_idx] = gamma
            T = np.exp(-(gamma_nm / m) * pdt)

            # ThermoLangevin: S = sqrt(kB*Tsys * (1 - T^2))
            # NOTE: this S is the coefficient in mass-scaled space (same as i-PI)
            S = np.sqrt(kbt * (1.0 - T * T))

            # --- goes in a single step to mass scaled coordinates and applies damping ---
            # i-PI: p = dstrip(self.p) * dstrip(self.T_on_sm)
            # with T_on_sm = T / sm
            p = p * (T / sm)

            # --- ThermoLangevin deltah bookkeeping (thermostat work) ---
            # ThermoLangevin: deltah = noddot(p,p)/(T^2)  # must correct for the "pre-damping"
            # Here T is array, so do elementwise divide before summing.
            deltah = np.sum((p * p) / (T * T))

            # ThermoLangevin: p += S * prng.gvec(...)
            p += S * self.prng.gvec(p.shape)

            # ThermoLangevin: deltah -= noddot(p,p)
            deltah -= np.sum(p * p)

            # ThermoLangevin: self.p[:] = p * sm
            p = p * sm
            self.nm.pnm[:] = p
            self.beads.p = self.nm.transform.nm2b(p)

            # THermoLangvin: self.ethermo += deltah * 0.5.
            # The *net* bath exchange is exactly 0.5*deltah, so track it consistently:
            # Convention: + means energy to bath (same sign spirit as i-PI's ethermo).
            self.ediss += 0.5 * deltah

            info("after: " + str(self.beads.p), verbosity.low)
            return

        # =========================
        # Variable friction: exact OU with frozen sigma over the substep
        # =========================

        # For ring-polymer momenta, i-PI thermostats typically use kB * (P*T).
        # If you explicitly want physical T instead, swap to self._kbt().
        # inside _step_markovian, variable_friction branch
        kbt = self._kbt_rp()  # usually correct for RP momenta

        sigma = self._get_sigma()  # (nbeads, nbath, ndof)
        nbeads, nbath, ndof = sigma.shape

        p = self.beads.p
        m = self.beads.m3
        sm = np.sqrt(m)

        # optional i-PI-esque net bath work bookkeeping (ThermoCL/ThermoGLE style)
        et = 0.0

        for b in range(nbeads):
            # build Lambda = sigma^T sigma  (ndof, ndof)
            sig = sigma[b, :, :]               # (nbath, ndof)
            Lam = sig.T @ sig                  # friction tensor in p-space: dp = -Lam v dt + ...

            inv_sm = 1.0 / sm[b, :]            # (ndof,)
            # mass-weighted A = M^{-1/2} Lam M^{-1/2}
            A = (inv_sm[:, None] * Lam) * inv_sm[None, :]

            # mass-scaled momentum
            s = p[b, :] * inv_sm

            # i-PI-style work accumulation: et += 1/2 s^2 (before), et -= 1/2 s^2 (after)
            et += 0.5 * float(np.dot(s, s))

            # diagonalize symmetric PSD A
            # (use eigh; A should be symmetric up to roundoff)
            A = 0.5 * (A + A.T)
            a, V = np.linalg.eigh(A)
            a = np.clip(a, 0.0, None)          # guard tiny negatives

            # OU coefficients along eigenmodes
            c = np.exp(-a * pdt)               # exact OU drift factors
            s2 = np.sqrt(1.0 - c * c)          # exact OU noise factors

            # transform, damp, add noise, transform back
            y = V.T @ s
            y = c * y + np.sqrt(kbt) * s2 * self.prng.gvec(ndof)
            s = V @ y

            et -= 0.5 * float(np.dot(s, s))

            # back to physical momentum
            p[b, :] = s * sm[b, :]

        self.beads.p[:] = p
        self.ediss += et   # this is really "ethermo-like" net exchange

        return




        # # # Euler-Maruyama (Kloeden pg 339)
        # info("Markov: before get sigma", verbosity.low)
        # # variable friction / seperable coupling: requires sigma per bead
        # sigma = self._get_sigma()  # (nbeads, nbath, ndof)
        # info("Markov: after get sigma", verbosity.low)

        # Gamma = sigma[0,:,:].T @ sigma[0,:,:]
        # info("Gamma shape "+str(np.shape(Gamma)), verbosity.low)
        # info("Gamma "+str(Gamma), verbosity.low)

        # nbeads, nbath,  ndof = sigma.shape

        # # Record energy before dissipative kick
        # p0 = self.beads.p.copy()
        # K0 = self._kinetic(p0, self.beads.m3)
        
        # # Dissipative force
        # v = self.beads.p / self.beads.m3
        # u = np.einsum("baj,bj->ba", sigma, v)
        # self.beads.p += -np.einsum("baj,ba->bj", sigma, u) * pdt

        # # Calculate addional dissipation energy
        # p1 = self.beads.p
        # K1 = self._kinetic(p1, self.beads.m3)
        # self.ediss += -(K1 - K0)

        # # Random force
        # amp = np.sqrt(2.0 * kbt * pdt)
        # for b in range(nbeads):
        #     g = self.prng.gvec(nbath)
        #     self.beads.p[b, :] += amp * (sigma[b, :, :].T @ g)

        # p2 = self.beads.p
        # K2 = self._kinetic(p2, self.beads.m3)
        # self.erand += (K2 - K1)
        return


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



dproperties(Friction, ["sigma", "friction_coupling_nm", "emf", "ediss","erand","fmf_nm", "fmf"])



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

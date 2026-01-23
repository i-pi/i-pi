"""Creates objects to deal with friction."""

# This file is part of i-PI.
# i-PI Copyright (C) 2025 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np

from ipi.engine.motion import Motion
from ipi.engine.normalmodes import NormalModes
from ipi.engine.beads import Beads
from ipi.utils.depend import depend_value, dproperties
from ipi.utils.messages import info, verbosity


class Friction:
    spectral_density: np.ndarray  # (n, 2)
    """Input spectral density of omega and J(omega) value pairs"""

    beads: Beads
    """Reference to the beads"""
    nm: NormalModes
    """Reference to the normal modes"""

    def __init__(
        self,
        spectral_density=np.zeros(0, float),
        alpha_input=np.zeros(0, float),
        position_dependent: bool = False,
    ):
        """Initialises the friction object.

        Args:
            spectral_density: Cosine transform of the time-dependent factor in the friction kernel,
                divided by frequncy. Supplied as a 2d array of two columns containing frequency and
                spectral density, respectively.
                Defaults to np.zeros(0, float).
            alpha_input: Normal-mode coefficients in expression for the frictional mean-field
                potential [Eq. (8b) in https://doi.org/10.1103/PhysRevLett.134.226201].
                Defaults to np.zeros(0, float).
            position_dependent (bool, optional): True if the gradient of the friction coupling F(q)
                [introduced in Eq. (5) of https://doi.org/10.1103/PhysRevLett.134.226201]
                depends on position.
                Defaults to False.

        """
        # TODO: make these depend objects?
        self.spectral_density = np.asanyarray(spectral_density, dtype=float).copy()
        self.alpha_input = np.asanyarray(alpha_input, dtype=float).copy()
        self.position_dependent = position_dependent
        # Diffusion coefficient
        self._Sigma = depend_value(name="Sigma", func=self.get_diffusion_coefficient)
        # Friction coupling: F(q), such that Σ{i,α} = ∂F(q) / ∂q{i,α}
        self._friction_coupling_nm = depend_value(
            name="friction_coupling_nm",
            func=self.get_friction_coupling_nm,
            dependencies=[self._Sigma],
        )
        # Frictional mean-field force
        self._ffric_nm = depend_value(
            name="ffric_nm",
            func=self.get_ffric_nm,
            dependencies=[self._friction_coupling_nm],
        )
        self._ffric = depend_value(
            name="ffric", func=self.get_ffric, dependencies=[self._ffric_nm]
        )
        # Frictional mean-field potential
        self._efric = depend_value(name="efric", value=0.0)

        if self.position_dependent:
            raise NotImplementedError(
                "Position dependent friction not implemented yet."
            )

    def bind(self, motion: Motion) -> None:
        self.beads = motion.beads
        self.nm = motion.nm
        self.forces = motion.ensemble.forces
        if self.alpha_input.shape[0] == self.nm.omegak.size:
            # if self.alpha is already provided as a file, use it
            self.alpha = self.alpha_input[:, 1]
            omegak_input = self.alpha_input[:, 0]
            if not np.allclose(omegak_input, self.nm.omegak):
                raise ValueError(
                    "The provided alpha values do not correspond to the current normal mode frequencies."
                )
        else:
            # otherwise, compute alpha numerically
            assert self.spectral_density.ndim == 2
            # TODO: should the user provide J or Λ? The integral over Λ starts at ω = 0, for
            # which the line below is numerically problematic.
            Lambda = self.spectral_density[:, 1] / self.spectral_density[:, 0]
            omega = self.spectral_density[:, 0]
            self.alpha = get_alpha_numeric(
                Lambda=Lambda,
                omega=omega,
                omegak=self.nm.omegak,
            )  # (nmodes,)
            info("compute alpha using get_alpha_numeric().")
        self._Sigma.add_dependency(self.forces._extras)
        self._friction_coupling_nm.add_dependency(self.beads._q)
        self._efric.add_dependency(self._friction_coupling_nm)
        self._efric._func = self.get_efric

    def get_diffusion_coefficient(self):
        Sigma = self.forces.extras.get("diffusion_coefficient")
        if Sigma is None:
            raise KeyError(f"Did not find 'diffusion_coefficient' among the force extras = {self.forces.extras}")
        else:
            return np.asarray(Sigma)

    def get_friction_coupling_nm(self):
        """Compute the friction coupling for each normal-mode index"""
        if self.position_dependent:
            raise NotImplementedError(
                "The calculation of friction coupling for position-dependent diffusion coefficients is not implemented."
            )
        else:
            # Here we assume that the interaction potential, F(Q) in https://doi.org/10.1103/PhysRevLett.134.226201,
            # is of the form F(q) = SUM[ c{i,α} q{i,α}, {{i,0,n_atom-1}, {α,0,2}} ] where α indexes Cartesian components
            # The diffusion coefficients for bead index l returned by the driver are expected to be packed as
            # Σ{i,α} = ∂F(q) / ∂q{i,α} = diffusion_coeff[l, 3*i+α].
            return np.sum(self.Sigma * self.nm.qnm, axis=-1)

    def get_efric(self):
        """Compute the frictional potential of mean field, Eq. (S19) of https://doi.org/10.1103/PhysRevLett.134.226201"""
        return np.sum(self.alpha * self.friction_coupling_nm**2) / 2

    def get_ffric_nm(self):
        """Negative derivative of the frictional potential of mean field with respect to normal modes"""
        return -(self.alpha * self.friction_coupling_nm)[:, np.newaxis] * self.Sigma

    def get_ffric(self):
        """Negative derivative of the frictional potential of mean field with respect to bead positions"""
        return self.nm.transform.nm2b(self.ffric_nm)

    def step(self, pdt: float) -> None:
        self.beads.p += self.ffric * pdt


dproperties(Friction, ["Sigma", "friction_coupling_nm", "efric", "ffric_nm", "ffric"])


def get_alpha_numeric(
    Lambda: np.ndarray, omega: np.ndarray, omegak: np.ndarray
) -> np.ndarray:
    try:
        from scipy.interpolate import CubicSpline
    except ModuleNotFoundError as e:
        raise ModuleNotFoundError(
            "Friction class requires scipy to work, please install scipy"
        ) from e

    alpha = np.zeros(omegak.shape)
    for idx, omegak in enumerate(omegak):
        # TODO: what if omega[0] > 0?
        f = CubicSpline(omega, Lambda * omegak**2 / (omega**2 + omegak**2))
        alpha[idx] = 2 / np.pi * f.integrate(0, omega[-1])
        info(
            f"for normal mode {omegak} alpha is {alpha[idx]}",
            verbosity.high,
        )
    return alpha

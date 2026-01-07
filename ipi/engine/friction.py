"""Creates objects to deal with friction."""

# This file is part of i-PI.
# i-PI Copyright (C) 2025 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np

from ipi.engine.motion import Motion
from ipi.engine.normalmodes import NormalModes
from ipi.engine.beads import Beads
from ipi.utils.depend import depend_value
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
        efric=0.0,
        position_dependent: bool = False,
    ):
        """Initialises friction

        Args:
        spectral density:
        alpha:
        efric: The initial friction energy.
            Default to 0.0. It will be non-zero if the friction class is initialised from a checkpoint file.
        """
        self.spectral_density = np.asanyarray(spectral_density, dtype=float)
        self.alpha_input = np.asanyarray(alpha_input, dtype=float)
        self._efric = depend_value(name="efric", value=efric)
        self.position_dependent = position_dependent

        if self.position_dependent:
            raise NotImplementedError(
                "Position dependent friction not implemented yet."
            )

    def bind(self, motion: Motion) -> None:
        self.beads = motion.beads
        self.nm = motion.nm
        self.forces = motion.ensemble.forces

        # if self.alpha is already provided as a file, use it
        # if self.alpha.size ==self.nm.omegak.size
        if self.alpha_input.shape[0] == self.nm.omegak.size:
            self.alpha = self.alpha_input[:, 1]
            omegak_input = self.alpha_input[:, 0]
            if not np.allclose(omegak_input, self.nm.omegak):
                raise ValueError(
                    "The provided alpha values do not correspond to the current normal mode frequencies."
                )
        else:
            assert self.spectral_density.ndim == 2
            Lambda = self.spectral_density[:, 1] / self.spectral_density[:, 0]
            omega = self.spectral_density[:, 0]

            # otherwise, compute alpha numerically

            self.alpha = get_alpha_numeric(
                Lambda=Lambda,
                omega=omega,
                omegak=self.nm.omegak,
            )  # (nmodes,)
            info("compute alpha using get_alpha_numeric().")

    def fric_forces(self) -> np.ndarray:
        fnm = -self.alpha[:, np.newaxis] * self.nm.qnm  # (nmodes, 3 * natoms)
        forces = self.nm.transform.nm2b(fnm)  # (nbeads, 3 * natoms)
        eta0 = np.asarray(self.forces.extras["eta0"])
        if self.position_dependent:
            ...  # To be implemented
        return forces * eta0[:, np.newaxis]

    def step(self, pdt: float) -> None:
        fric_forces = self.fric_forces()
        self.beads.p += fric_forces * pdt
        self.efric = 0.5 * np.einsum("n,nm,nm->", self.alpha, self.nm.qnm, self.nm.qnm)


# dproperties(Friction, ["efric"])


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
        f = CubicSpline(omega, Lambda * omegak**2 / (omega**2 + omegak**2))
        alpha[idx] = 2 / np.pi * f.integrate(0, omega[-1])
        info(
            f"for normal mode {omegak} alpha is {alpha[idx]}",
            verbosity.high,
        )
    return alpha

"""Creates objects to deal with friction."""

# This file is part of i-PI.
# i-PI Copyright (C) 2025 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np
from typing import Protocol

from ipi.engine.motion import Motion
from ipi.engine.normalmodes import NormalModes
from ipi.engine.beads import Beads


class FrictionProtocol(Protocol):
    def bind(self, motion: Motion) -> None:
        ...

    def forces(self) -> np.ndarray:
        ...

class Friction(FrictionProtocol):
    spectral_density: np.ndarray  # (n, 2)
    """Input spectral density of omega and J(omega) value pairs"""

    beads: Beads
    """Reference to the beads"""
    nm: NormalModes
    """Reference to the normal modes"""

    def __init__(
        self,
        spectral_density=np.zeros(0, float),
    ):
        self.spectral_density = np.asanyarray(spectral_density, dtype=float)

    def bind(self, motion: Motion) -> None:
        self.beads = motion.beads
        self.nm = motion.nm
        assert self.spectral_density.ndim == 2
        Lambda = self.spectral_density[:, 1] / self.spectral_density[:, 0]
        omega = self.spectral_density[:, 0]
        self.alpha = get_alpha_numeric(
            Lambda=Lambda,
            omega=omega,
            nm=self.nm,
        )  # (nmodes,)

    def forces(self) -> np.ndarray:

        fnm = self.alpha[:, np.newaxis] * self.nm.qnm  # (nmodes, 3 * natoms)
        forces = self.nm.transform.nm2b(fnm)  # (nbeads, 3 * natoms)
        return forces


def get_alpha_numeric(
    Lambda: np.ndarray, omega: np.ndarray, nm: NormalModes
) -> np.ndarray:
    try:
        from scipy.interpolate import CubicSpline
        from scipy.integrate import quad
    except ModuleNotFoundError as e:
        raise ModuleNotFoundError(
            "Friction class requires scipy to work, please install scipy"
        ) from e

    alpha = np.zeros(nm.omegak.shape)
    for idx, omegak in enumerate(nm.omegak):
        f = CubicSpline(omega, Lambda * omegak**2 / (omega**2 + omegak**2))
        alpha[idx] = 2 / np.pi * quad(f, 0, omega[-1])[0]
    return alpha

# def get_eta(beads: Beads, forces: Forces) -> np.ndarray:
#    """
#    Get the friction matrix from the forces object.

#    Returns
#    -------
#    np.ndarray
#        The friction matrix, shape (nbeads, 3 * natoms, 3 * natoms)
#    """
#    shape = (beads.nbeads, 3 * beads.natoms, 3 * beads.natoms)
#    return np.array(forces.extras["friction"]).reshape(shape)

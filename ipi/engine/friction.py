"""Creates objects to deal with friction."""

# This file is part of i-PI.
# i-PI Copyright (C) 2025 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np

# from scipy.linalg import sqrtm
# from scipy.interpolate import interp1d
from scipy.integrate import quad

from ipi.engine.motion import Motion
from ipi.engine.normalmodes import NormalModes
from ipi.engine.beads import Beads


class Friction:
    beads: Beads
    nm: NormalModes

    def __init__(
        self,
        frictionSD: bool = True,
        fric_spec_dens=np.zeros(0, float),
        fric_spec_dens_ener=0.0,
        bath_mass: float = 1.0,
        eta0: float = 0.0,
        omega_cutoff: float = 0.0,
        *args,
        **kwargs,
    ):
        self.has_numerical_spectral_density = frictionSD
        self.Lambda = fric_spec_dens / fric_spec_dens_ener
        self.omega = fric_spec_dens_ener
        self.eta0 = eta0
        self.omega_cutoff = omega_cutoff

    def bind(self, motion: Motion) -> None:
        self.beads = motion.beads
        self.nm = motion.nm

    def forces(self) -> np.ndarray:
        if self.has_numerical_spectral_density:
            alpha = get_alpha_numeric(
                Lambda=self.Lambda,
                omega=self.omega,
                nm=self.nm,
            )  # (nmodes,)
        else:
            alpha = get_alpha(
                eta0=self.eta0,
                omega_cutoff=self.omega_cutoff,
                nm=self.nm,
            )
        fnm = alpha[:, np.newaxis] * self.nm.qnm  # (nmodes, 3 * natoms)
        forces = self.nm.transform.nm2b(fnm)  # (nbeads, 3 * natoms)
        return forces


def get_alpha_numeric(
    Lambda: np.ndarray, omega: np.ndarray, nm: NormalModes
) -> np.ndarray:
    Lambda = Lambda[:, np.newaxis]
    omega2 = omega[:, np.newaxis] ** 2
    omegak2 = nm.omegak[np.newaxis, :] ** 2
    alpha = 2 / np.pi * np.sum(Lambda * omegak2 / (omega2 + omegak2), axis=1)
    return alpha


def get_alpha(eta0: float, omega_cutoff: float, nm: NormalModes) -> np.ndarray:

    def Lambda(omega: float) -> float:
        return eta0 * np.exp(-omega / omega_cutoff)

    alpha = np.zeros_like(nm.omegak)
    for idx, omegak in enumerate(nm.omegak):

        def Value(omega: float) -> float:
            return Lambda(omega) * omegak**2 / (omega**2 + omegak**2)

        # Integrate from 0 to infinity
        alpha[idx] = 2 / np.pi * quad(Value, 0, np.inf)[0]
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

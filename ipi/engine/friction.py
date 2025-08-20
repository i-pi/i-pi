"""Creates objects to deal with friction."""

# This file is part of i-PI.
# i-PI Copyright (C) 2025 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np

from scipy.interpolate import CubicSpline
from scipy.integrate import quad

from ipi.engine.motion import Motion
from ipi.engine.normalmodes import NormalModes
from ipi.engine.beads import Beads


class Friction:
    beads: Beads
    nm: NormalModes

    def __init__(
        self,
        fric_spec_dens=np.zeros(0, float),
        fric_spec_dens_ener=0.0,
        *args,
        **kwargs,
    ):
        self.Lambda = fric_spec_dens / fric_spec_dens_ener
        self.omega = fric_spec_dens_ener

    def bind(self, motion: Motion) -> None:
        self.beads = motion.beads
        self.nm = motion.nm
        self.alpha = get_alpha_numeric(
            Lambda=self.Lambda,
            omega=self.omega,
            nm=self.nm,
        )  # (nmodes,)

    def forces(self, time: float) -> np.ndarray:
        del time  # unused for now

        fnm = self.alpha[:, np.newaxis] * self.nm.qnm  # (nmodes, 3 * natoms)
        forces = self.nm.transform.nm2b(fnm)  # (nbeads, 3 * natoms)
        return forces


def get_alpha_numeric(
    Lambda: np.ndarray, omega: np.ndarray, nm: NormalModes
) -> np.ndarray:
    alpha = np.zeros_like(nm.omegak)
    for idx, omegak in enumerate(nm.omegak):
        f = CubicSpline(omega, Lambda * omegak**2 / (omega**2 + omegak**2))
        alpha[idx] = 2 / np.pi * quad(f, 0, np.inf)[0]
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

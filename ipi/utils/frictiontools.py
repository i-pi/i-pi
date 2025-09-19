from typing import Union
from ipi.engine.friction import FrictionProtocol

import numpy as np
from scipy.integrate import quad
from scipy.special import sici
from ipi.engine.normalmodes import NormalModes
from ipi.engine.motion import Motion


class Fricition_eq133(FrictionProtocol):
    def __init__(self, omega_cutoff: float, eta: float):
        self.omega_cutoff = omega_cutoff
        self.eta = eta

    def bind(self, motion: Motion) -> None:
        self.beads = motion.beads
        self.nm = motion.nm

        alpha = np.zeros(self.nm.omegak.shape)
        for idx, omegak in enumerate(self.nm.omegak):
            z = omegak / self.omega_cutoff
            si, ci = sici(z)
            alpha[idx] = (
                2
                / np.pi
                * self.eta
                * omegak
                * (ci * np.sin(z) - (si - np.pi * 0.5) * np.cos(z))
            )
        self.alpha = alpha

    def forces(self) -> np.ndarray:

        fnm = self.alpha[:, np.newaxis] * self.nm.qnm  # (nmodes, 3 * natoms)
        forces = self.nm.transform.nm2b(fnm)  # (nbeads, 3 * natoms)
        return forces


class Fricition_eq134(FrictionProtocol):
    def __init__(self, omega_cutoff: float, eta: float):
        self.omega_cutoff = omega_cutoff
        self.eta = eta

    def bind(self, motion: Motion) -> None:
        self.beads = motion.beads
        self.nm = motion.nm

        alpha = np.zeros(self.nm.omegak.shape)
        for idx, omegak in enumerate(self.nm.omegak):
            z = omegak / self.omega_cutoff
            alpha[idx] = (
                2
                / np.pi
                * self.eta
                * omegak
                * (z * (np.euler_gamma + np.log(z)) - (z - np.pi * 0.5))
            )
        self.alpha = alpha

    def forces(self) -> np.ndarray:

        fnm = self.alpha[:, np.newaxis] * self.nm.qnm  # (nmodes, 3 * natoms)
        forces = self.nm.transform.nm2b(fnm)  # (nbeads, 3 * natoms)
        return forces


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


def expohmic_J(
    omega: Union[float, np.ndarray], eta, omega_cut
) -> Union[float, np.ndarray]:
    """Spectral density at frequency omega."""
    return eta * np.abs(omega) * np.exp(-np.abs(omega) / omega_cut)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="tool to compute the spectral densities for the electronic friction"
    )
    parser.add_argument("--omega-cut", required=True, type=float)
    parser.add_argument("--eta", required=True, type=float)
    parser.add_argument("output", type=str)
    args = parser.parse_args()

    omega_cut = args.omega_cut

    eta = args.eta
    omega = np.arange(omega_cut)

    J = expohmic_J(omega, eta, omega_cut)

    print(np.stack((omega, J)))
    with open(args.output, "w") as fd:
        for i in range(J.size):
            fd.write(f"{omega[i]} {J[i]}\n")

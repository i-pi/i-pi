import numpy as np
from scipy.integrate import quad

from ipi.engine.normalmodes import NormalModes


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

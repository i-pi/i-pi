from typing import Union

import numpy as np
from scipy.special import sici


def get_alpha_ohmic(omegak, omega_cutoff, eta):
    r"""Analytical expression of alpha based on the exponentially damped ohmic spectral density J.

    .. math::

       \eta \omega_n^2 \int_0^\infty \frac{e^{-\omega/\omega_c}}{\omega^2 + \omega_n^2} \, d\omega
       = \eta \omega_n \big[ \text{Ci}(z_n)\sin(z_n) - (\text{Si}(z_n) - \tfrac{\pi}{2}) \cos(z_n) \big]
    """

    alpha = np.zeros(omegak.shape)
    for idx, omegak in enumerate(omegak):
        z = omegak / omega_cutoff
        si, ci = sici(z)
        alpha[idx] = (
            2 / np.pi * eta * omegak * (ci * np.sin(z) - (si - np.pi * 0.5) * np.cos(z))
        )
    return alpha


def expohmic_J(
    omega: Union[float, np.ndarray], eta, omega_cut
) -> Union[float, np.ndarray]:
    r"""Spectral density at frequency omega.
    Within this function the exp ohmic J can be calculated by providing eta and omegacut

    .. math::

       J(\omega) = \eta \, \omega \, e^{-\omega / \omega_c}

    attributes:
      eta:  is a coupling strength (related to the friction coefficient at low frequency)
      omega_cut:   is a cutoff frequency describing how quickly the coupling decays for high-frequency modes


    python frictiontools.py --omega-cut OMEGA_CUT --eta ETA output
    """
    return np.abs(omega) * expohmic_Lambda(omega, eta, omega_cut)


def expohmic_Lambda(
    omega: Union[float, np.ndarray], eta, omega_cut
) -> Union[float, np.ndarray]:
    r"""Spectrum at frequency omega.
    Within this function the exp ohmic Lambda can be calculated by providing eta and omegacut

    .. math::

       Lambda(\omega) = \eta \,  e^{-\omega / \omega_c}

    attributes:
      eta:  is a coupling strength (related to the friction coefficient at low frequency)
      omega_cut:   is a cutoff frequency describing how quickly the coupling decays for high-frequency modes


    python frictiontools.py --omega-cut OMEGA_CUT --eta ETA output
    """
    return eta * np.exp(-np.abs(omega) / omega_cut)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="tool to compute the spectral densities for the electronic friction"
    )
    parser.add_argument("--omega-cut", required=True, type=float)
    parser.add_argument("--eta", required=True, type=float)
    parser.add_argument("output", type=str)
    parser.add_argument("--points", type=int, default=999)
    args = parser.parse_args()

    omega_cut = args.omega_cut

    eta = args.eta
    omega = np.linspace(0, omega_cut, args.points + 1)[1:]

    J = expohmic_J(omega, eta, omega_cut)

    print(np.stack((omega, J)))
    with open(args.output, "w") as fd:
        for i in range(J.size):
            fd.write(f"{omega[i]} {J[i]}\n")

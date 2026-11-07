from typing import Union

import numpy as np
from scipy.special import sici


def get_alpha_ohmic(omegak, omega_cutoff, eta):
    r"""Analytical expression of alpha based on the exponentially damped ohmic spectral density J.

    .. math::

       \eta \omega_n^2 \int_0^\infty \frac{e^{-\omega/\omega_c}}{\omega^2 + \omega_n^2} \, d\omega
       = \eta \omega_n \big[ \text{Ci}(z_n)\sin(z_n) - (\text{Si}(z_n) - \tfrac{\pi}{2}) \cos(z_n) \big]
    """

    omegak = np.asarray(omegak, dtype=float)
    if (
        not np.isfinite(omega_cutoff)
        or omega_cutoff <= 0
        or not np.isfinite(eta)
        or eta < 0
        or not np.all(np.isfinite(omegak))
        or np.any(omegak < 0)
    ):
        raise ValueError(
            "Require positive cutoff and finite nonnegative frequencies and eta."
        )
    alpha = np.zeros(omegak.shape)
    positive = omegak > 0
    frequency = omegak[positive]
    z = frequency / omega_cutoff
    si, ci = sici(z)
    alpha[positive] = (
        2 / np.pi * eta * frequency * (ci * np.sin(z) - (si - np.pi * 0.5) * np.cos(z))
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


    python frictiontools.py --omega-cut OMEGA_CUT --eta ETA --quantity J output
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
        description="Generate an exponential Ohmic spectrum in atomic units."
    )
    parser.add_argument(
        "--omega-cut",
        required=True,
        type=float,
        help="Exponential decay frequency, not the grid endpoint.",
    )
    parser.add_argument("--eta", required=True, type=float)
    parser.add_argument(
        "--omega-max",
        type=float,
        help="Grid endpoint; default 20 times omega-cut. Check cutoff convergence.",
    )
    parser.add_argument(
        "--quantity",
        choices=["Lambda", "J"],
        default="Lambda",
        help="Lambda is the spectrum expected by the friction input.",
    )
    parser.add_argument("output", type=str)
    parser.add_argument("--points", type=int, default=10001)
    args = parser.parse_args()
    maximum = args.omega_max if args.omega_max is not None else 20 * args.omega_cut
    if (
        not np.isfinite(args.omega_cut)
        or args.omega_cut <= 0
        or not np.isfinite(args.eta)
        or args.eta < 0
        or not np.isfinite(maximum)
        or maximum <= 0
        or args.points < 2
    ):
        parser.error(
            "Require finite positive frequencies, nonnegative eta and at least two points."
        )
    omega = np.linspace(0, maximum, args.points)
    spectrum = (expohmic_Lambda if args.quantity == "Lambda" else expohmic_J)(
        omega, args.eta, args.omega_cut
    )
    np.savetxt(
        args.output,
        np.column_stack((omega, spectrum)),
        header=f"omega, {args.quantity}(omega); atomic units",
    )

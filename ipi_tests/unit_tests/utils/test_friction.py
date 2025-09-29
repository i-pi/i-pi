"""Test electron-friction implementations. """

import numpy as np
from ipi.utils.frictiontools import get_alpha_eq133


def test_analytical_alpha():
    r"""
    In this test we are calculating the alpha based on the analytical ohmic spectral density J. Eq 1.33 from George draft.

    .. math::
       \eta \omega_n^2 \int_0^\infty \frac{e^{-\omega/\omega_c}}{\omega^2 + \omega_n^2} \, d\omega
       = \eta \omega_n \big[ \text{Ci}(z_n)\sin(z_n) - (\text{Si}(z_n) - \tfrac{\pi}{2}) \cos(z_n) \big]
    """
    omegak = np.array(
        [
            0.0,
            0.005700267359999999,
            0.009873152684246512,
            0.01140053472,
            0.009873152684246512,
            0.005700267359999999,
        ]
    )
    omega_cutoff = 2.0
    eta = 3.0

    alpha = get_alpha_eq133(omegak, omega_cutoff, eta)
    assert np.allclose(
        alpha,
        np.array(
            [
                0.0,
                5.480081735587132e-05,
                0.00016439167542196503,
                0.00021918171683906422,
                0.00016439167542196503,
                5.480081735587132e-05,
            ]
        ),
    )

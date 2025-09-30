"""Test electron-friction implementations. """

import pytest
import numpy as np
from ipi.utils.frictiontools import get_alpha_eq133
from ipi.utils.frictiontools import get_alpha_eq134


@pytest.fixture
def omega_cutoff() -> float:
    return 2.0


@pytest.fixture
def eta() -> float:
    return 0.5


def test_analytical_alphaeq133(omega_cutoff: float, eta: float):
    r"""
    In this test we are calculating the alpha based on the analytical ohmic spectral density J. Eq 1.33 from George draft.

    .. math::
       \eta \omega_n^2 \int_0^\infty \frac{e^{-\omega/\omega_c}}{\omega^2 + \omega_n^2} \, d\omega
       = \eta \omega_n \big[ \text{Ci}(z_n)\sin(z_n) - (\text{Si}(z_n) - \tfrac{\pi}{2}) \cos(z_n) \big]
    """
    omegak = np.array(
        [
            0.005700267359999999,
            0.009873152684246512,
            0.01140053472,
            0.009873152684246512,
            0.005700267359999999,
        ]
    )

    alpha = get_alpha_eq133(omegak, omega_cutoff, eta)
    assert np.allclose(
        alpha,
        np.array([0.00281764, 0.00484762, 0.00558463, 0.00484762, 0.00281764]),
    )


def test_analytical_alphaeq134(omega_cutoff: float, eta: float):
    r"""
    In this test we are calculating alpha based on the analytical solution of the eq1.33 at small values of zn (zn=ωn/ωc).
    """
    omegak = np.array(
        [
            0.005700267359999999,
            0.009873152684246512,
            0.01140053472,
            0.009873152684246512,
            0.005700267359999999,
        ]
    )
    alpha = get_alpha_eq134(omegak, omega_cutoff, eta)
    assert np.allclose(
        alpha,
        np.array([0.00281764, 0.00484762, 0.00558463, 0.00484762, 0.00281764]),
    )

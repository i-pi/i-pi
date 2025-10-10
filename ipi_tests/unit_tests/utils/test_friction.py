"""Test electron-friction implementations.


This test can be executed using the following command:

pytest -v test_friction.py

or 

pytest -s test_friction.py

where test_friction.py is located in "i-pi/ipi_tests/unit_tests/utils".

In this test, the analytical and numerical evaluations of the alpha parameter within the friction class and frictiontool.py are performed.
The model system consists of a hydrogen atom represented by six beads, where the normal mode frequencies (omega_k, as defined in i-PI) are obtained by running i-PI (i-pi input.xml) in the directory "ipi_tests/regression_tests/tests/NVE/NVE_with_friction/double_well".

These normal mode frequencies are then used as input for Eqs. (1.33) and (1.34):

 - Eq. (1.33) represents the analytical expression for alpha, assuming an exponential Ohmic spectral density J(ω). The zero-frequency term is omitted to avoid undefined (NaN) values.

 - Eq. (1.34) provides the asymptotic analytical approximation of Eq. (1.33).

The reference values correspond to the results of Eq. (1.34).
During the unit test, alpha is computed using:

 - The analytical formulations defined in frictiontool.py ("i-pi/ipi/utils/frictiontools.py"), and

 - The numerical implementation in friction.py ("i-pi/ipi/engine/friction.py").

The computed results are then compared against the reference values.
For the analytical case (Eq. 1.33), the deviation is within 1e-8, while for the numerical implementation, the deviation remains within 1e-5.


"""

import pytest
import numpy as np
from ipi.utils.frictiontools import get_alpha_eq133
from ipi.utils.frictiontools import get_alpha_eq134
from ipi.engine.friction import get_alpha_numeric
from ipi.utils.frictiontools import expohmic_J


@pytest.fixture
def omega_cutoff() -> float:
    return 2.0


@pytest.fixture
def eta() -> float:
    return 0.5


OMEGAK = np.array(
    [
        0.005700267359999999,
        0.009873152684246512,
        0.01140053472,
        0.009873152684246512,
        0.005700267359999999,
    ]
)

ALPHAK = np.array([0.00281764, 0.00484762, 0.00558463, 0.00484762, 0.00281764])


def test_analytical_alphaeq133(omega_cutoff: float, eta: float):
    r"""
    In this test we are calculating the alpha based on the analytical ohmic spectral density J. (Eq 1.33 from George draft.)

    .. math::
       \eta \omega_n^2 \int_0^\infty \frac{e^{-\omega/\omega_c}}{\omega^2 + \omega_n^2} \, d\omega
       = \eta \omega_n \big[ \text{Ci}(z_n)\sin(z_n) - (\text{Si}(z_n) - \tfrac{\pi}{2}) \cos(z_n) \big]
    """

    alpha = get_alpha_eq133(OMEGAK, omega_cutoff, eta)
    assert np.allclose(
        alpha,
        ALPHAK,
        atol=1e-7,
    )


def test_analytical_alphaeq134(omega_cutoff: float, eta: float):
    r"""
    In this test we are calculating alpha based on the analytical solution of the eq1.33 at small values of zn (zn=ωn/ωc).Eq 1.34 given as
    .. math::
    \sim \eta \omega_n \big[ z_n (\gamma + \ln z_n) - z_n - \tfrac{\pi}{2} \big]
     at small values of zn tis the asymptotic answer to Eq1.33.
    """

    alpha = get_alpha_eq134(OMEGAK, omega_cutoff, eta)
    assert np.allclose(
        alpha,
        ALPHAK,
        atol=1e-7,
    )


def test_numerical_alpha(omega_cutoff: float, eta: float):
    r"""In this test, the numerical evaluation of alpha, as implemented in ipi/engine/friction, is carried out.
    The computed numerical values of alpha are subsequently compared with the corresponding analytical expression given in Eq. 1.34.
    """

    omega = np.arange(0.0001, omega_cutoff, 0.0001)
    J = expohmic_J(omega, eta, omega_cutoff)
    Lambda = J / omega
    print(Lambda)
    alpha = get_alpha_numeric(Lambda, omega, OMEGAK)
    assert np.allclose(
        alpha,
        ALPHAK,
        atol=1e-5,  # Note: this is not a very great accuracy. Can this be made better?
    )
    print(alpha, ALPHAK)

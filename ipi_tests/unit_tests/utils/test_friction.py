"""Test electron-friction implementations.


This test can be executed using the following command:

pytest -v test_friction.py

where test_friction.py is located in "i-pi/ipi_tests/unit_tests/utils".

In this test, the analytical and numerical evaluations of the alpha parameter within the friction class and frictiontool.py are performed.
The model system consists of a hydrogen atom represented by six beads, where the normal mode frequencies (omega_k, as defined in i-PI) are obtained by running i-PI.

These normal mode frequencies are then used as input for two different expressions for the ohmic spectral density

 - `get_alpha_ohmic` represents the analytical expression for alpha, assuming an exponential Ohmic spectral density J(ω). The zero-frequency term is omitted to avoid undefined (NaN) values.

The reference values correspond to the results of the asymptotic version.
During the unit test, alpha is computed using:

 - The analytical formulations defined in frictiontool.py ("i-pi/ipi/utils/frictiontools.py"), and

 - The numerical implementation in friction.py ("i-pi/ipi/engine/friction.py").

The computed results are then compared against the reference values.
For the analytical case, the deviation is within 1e-8, while for the numerical implementation, the deviation remains within 1e-5.


"""

import pytest
import numpy as np
from ipi.utils.frictiontools import expohmic_Lambda, get_alpha_ohmic

from ipi.engine.friction import get_alpha_numeric


@pytest.fixture
def omega_cutoff() -> float:
    return 2.0


@pytest.fixture
def eta() -> float:
    return 0.5


@pytest.fixture
def omegak() -> np.ndarray:
    return np.array(
        [
            0.005700267359999999,
            0.009873152684246512,
            0.01140053472,
            0.009873152684246512,
            0.005700267359999999,
        ]
    )


@pytest.fixture
def alphak(omega_cutoff: float, eta: float, omegak: np.ndarray) -> np.ndarray:
    """Reference alpha_k values for the test cases computed with analytical formula."""
    return get_alpha_ohmic(omegak, omega_cutoff, eta)


def test_numerical_alpha(
    omega_cutoff: float, eta: float, omegak: np.ndarray, alphak: np.ndarray
) -> None:
    r"""In this test, the numerical evaluation of alpha, as implemented in ipi/engine/friction, is carried out.
    The computed numerical values of alpha are subsequently compared with the corresponding analytical expression.
    """

    omega = np.linspace(0, omega_cutoff, 1000000)
    Lambda = expohmic_Lambda(omega, eta, omega_cutoff)
    # print(Lambda)
    alpha = get_alpha_numeric(Lambda, omega, omegak)
    assert np.allclose(
        alpha,
        alphak,
        atol=1e-5,
    )
    print(alpha, alphak)

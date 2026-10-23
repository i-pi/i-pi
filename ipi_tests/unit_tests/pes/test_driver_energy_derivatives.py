"""Potential forces checked by displacing coordinates, including explicit baths."""

import numpy as np
import pytest
from ipi.pes.doublewell import DoubleWell_driver
from ipi.pes.doublewell_with_bath import DoubleWell_with_explicit_bath_driver
from ipi.pes.doublewell_with_friction import DoubleWell_with_friction_driver
from ipi.pes.bath import Harmonic_Bath_explicit


@pytest.mark.parametrize(
    "kind",
    [
        "doublewell",
        pytest.param(
            "bath",
            marks=pytest.mark.xfail(
                strict=True,
                raises=TypeError,
                reason="Double-well bath driver calls a bath object with no __call__ method",
            ),
        ),
        "friction",
    ],
)
def test_driver_force_is_negative_energy_derivative(kind):
    common = dict(w_b=500.0, v0=2000.0, m=1837.0, delta=0.0)
    if kind == "doublewell":
        driver = DoubleWell_driver(**common)
    elif kind == "bath":
        driver = DoubleWell_with_explicit_bath_driver(
            **common, eta0=0.2, eps1=0.3, eps2=0.1, deltaQ=0.7, w_c=500.0
        )
    else:
        driver = DoubleWell_with_friction_driver(
            **common, eta0=0.2, eps1=0.3, eps2=0.1, deltaQ=0.7
        )
    q = np.array([[0.3, -0.2, 0.1]])
    cell = np.eye(3) * 10
    energy, force, _, _ = driver(cell, q)
    assert np.isfinite(energy).all()
    for i in range(q.size):
        delta = np.zeros_like(q)
        delta.flat[i] = 1e-5
        gradient = (driver(cell, q + delta)[0] - driver(cell, q - delta)[0]) / (2e-5)
        np.testing.assert_allclose(force.flat[i], -gradient, rtol=1e-6, atol=1e-9)


@pytest.mark.parametrize("cutoff", [0.0, 0.7])
def test_explicit_bath_coupling_derivative_and_forces(cutoff):
    bath = Harmonic_Bath_explicit(3, 2.0, 0.4, 0.3, 0.2, 0.8, cutoff)
    q = np.array([0.2])
    x = np.array([[0.1], [-0.3], [0.4]])
    eps = 1e-6
    np.testing.assert_allclose(
        bath.dSD_dq(q), (bath.S(q + eps) - bath.S(q - eps)) / (2 * eps), atol=1e-10
    )
    energy, fq, fx = bath.compute_structure(q, x)
    np.testing.assert_allclose(
        fq,
        -(bath.compute_structure(q + eps, x)[0] - bath.compute_structure(q - eps, x)[0])
        / (2 * eps),
        atol=1e-10,
    )
    for i in range(3):
        delta = np.zeros_like(x)
        delta[i] = eps
        np.testing.assert_allclose(
            fx[i],
            -(
                bath.compute_structure(q, x + delta)[0]
                - bath.compute_structure(q, x - delta)[0]
            )
            / (2 * eps),
            atol=1e-10,
        )
    assert np.all(energy >= 0)

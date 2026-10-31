"""Potential forces checked by displacing coordinates, including explicit baths."""

import numpy as np
import pytest
from ipi.pes.doubledoublewell import DDW_with_explicit_bath_driver
from ipi.pes.doublewell import DoubleWell_driver
from ipi.pes.doublewell_with_bath import DoubleWell_with_explicit_bath_driver
from ipi.pes.doublewell_with_friction import DoubleWell_with_friction_driver
from ipi.pes.bath import Harmonic_Bath_explicit


@pytest.mark.parametrize("kind", ["doublewell", "bath", "friction", "double-bath"])
def test_driver_force_is_negative_energy_derivative(kind):
    common = dict(w_b=500.0, v0=2000.0, m=1837.0, delta=0.0)
    bath = dict(eta0=0.2, eps1=0.3, eps2=0.1, deltaQ=0.7)
    q = np.array([[0.3, -0.2, 0.1]])
    if kind == "doublewell":
        driver = DoubleWell_driver(**common)
    elif kind == "bath":
        driver = DoubleWell_with_explicit_bath_driver(**common, **bath, w_c=500.0)
    elif kind == "friction":
        driver = DoubleWell_with_friction_driver(**common, **bath)
    else:
        driver = DDW_with_explicit_bath_driver(
            wb1=500.0,
            v1=2000.0,
            wb2=400.0,
            v2=1500.0,
            C=1e-5,
            m=1837.0,
            delta=0.0,
            **bath,
            w_c=500.0
        )
        q = np.array([[0.3, -0.2, 0.1, 0.2, -0.1, 0.05]])
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

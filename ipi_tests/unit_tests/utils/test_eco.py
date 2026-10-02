"""Tests the economised (eco) path integral frequencies."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import pytest
import numpy as np

from ipi.utils import nmtransform


def rms_r2_error(eva, nbeads, xmax, npts=2000):
    """RMS fractional error in the radius of gyration of harmonic oscillators
    with 0 <= beta*hbar*omega <= xmax, for dimensionless mode eigenvalues eva."""

    t, w = np.polynomial.legendre.leggauss(npts)
    x = 0.5 * xmax * (t + 1.0)
    y = eva[1:] * nbeads  # y_k = beta*hbar*omega_k
    f = nmtransform._eco_f(x)
    r = f * (1.0 / (x[:, np.newaxis] ** 2 + y**2)).sum(axis=1) - 1.0
    return np.sqrt(0.5 * (w * r**2).sum())


@pytest.mark.parametrize("nbeads", [2, 3, 4, 8, 16, 32, 33])
@pytest.mark.parametrize("xmax", [2.0, 10.0, 50.0])
def test_eco_eva_structure(nbeads, xmax):
    """Checks centroid, symmetry and positivity of the eco eigenvalues."""

    eva = nmtransform.eco_eva(nbeads, xmax)
    assert eva.shape == (nbeads,)
    assert eva[0] == 0.0
    assert np.all(eva[1:] > 0)
    assert np.all(np.isfinite(eva))
    # internal modes k and nbeads-k must be degenerate
    np.testing.assert_allclose(eva[1:], eva[1:][::-1])


@pytest.mark.parametrize("nbeads", [8, 16, 32, 33, 64])
@pytest.mark.parametrize("xmax", [2.0, 10.0, 50.0])
def test_eco_eva_accuracy(nbeads, xmax):
    """Checks that eco frequencies reproduce harmonic radii of gyration
    better than the Trotter ones."""

    err_eco = rms_r2_error(nmtransform.eco_eva(nbeads, xmax), nbeads, xmax)
    err_trotter = rms_r2_error(nmtransform.nm_eva(nbeads), nbeads, xmax)
    assert err_eco <= err_trotter
    # well-converged regime: eco error should be dramatically smaller
    if nbeads >= 4 * xmax:
        assert err_eco < 1e-2 * err_trotter


def test_eco_f_small_x():
    """Checks the small-x expansion of the objective kernel is smooth."""

    x = np.array([1e-8, 0.1, 0.4999, 0.5001, 1.0])
    f = nmtransform._eco_f(x)
    assert np.all(np.isfinite(f))
    np.testing.assert_allclose(f[0], 12.0, rtol=1e-10)
    # continuity across the series/direct switchover
    assert abs(f[2] - f[3]) < 1e-4


def test_eco_eva_invalid_xmax():
    """Checks that non-positive maximum frequencies are rejected."""

    with pytest.raises(ValueError):
        nmtransform.eco_eva(8, 0.0)
    with pytest.raises(ValueError):
        nmtransform.eco_eva(8, -1.0)


def test_eco_eva_classical_limit():
    """A single bead has no springs regardless of the fit."""

    assert np.all(nmtransform.eco_eva(1, 10.0) == 0.0)


def test_eco_eva_warm_start():
    """A fit warm-started from a nearby solution must match a cold fit."""

    nbeads = 16
    cold = nmtransform.eco_eva(nbeads, 20.0)
    y0 = cold[1 : nbeads // 2 + 1] * nbeads  # previous dimensionless solution
    warm = nmtransform.eco_eva(nbeads, 20.5, y0)
    ref = nmtransform.eco_eva(nbeads, 20.5)
    np.testing.assert_allclose(warm, ref, rtol=1e-6)


@pytest.mark.parametrize(
    "y0",
    [
        np.array([1.0, 2.0, 3.0]),  # wrong length
        np.zeros(8),  # not strictly positive
        np.linspace(50.0, 1.0, 8),  # descending
    ],
)
def test_eco_eva_bad_guess_falls_back(y0):
    """Invalid initial guesses are ignored, falling back to the Matsubara start."""

    ref = nmtransform.eco_eva(16, 20.0)
    np.testing.assert_allclose(nmtransform.eco_eva(16, 20.0, y0), ref, rtol=1e-8)


@pytest.mark.parametrize("nbeads", [2, 3, 8, 16, 33])
def test_spring_energy_parseval(nbeads):
    """Checks that the normal-mode spring energy with Trotter eigenvalues
    matches the bead-difference formula (the identity used to compute the
    primitive kinetic energy estimator in the normal-mode representation)."""

    rng = np.random.RandomState(31415)
    q = rng.normal(size=(nbeads, 6))
    cmat = nmtransform.mk_nm_matrix(nbeads)
    qnm = cmat @ q
    eva = nmtransform.nm_eva(nbeads)
    vnm = 0.5 * ((eva**2) @ (qnm**2))
    vbead = 0.5 * ((q - np.roll(q, 1, axis=0)) ** 2).sum(axis=0)
    np.testing.assert_allclose(vnm, vbead, rtol=1e-10)

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


def rms_open_errors(eva, nbeads, xmax, npts=2000):
    """RMS fractional errors of the open-path end-to-end variance (including
    the endpoint-kernel term 1/P, over its representable window) and radius
    of gyration, for dimensionless open-chain mode eigenvalues eva."""

    P = nbeads
    x = (np.arange(npts) + 0.5) * (xmax / npts)
    y = eva[1:] * P
    k = np.arange(1, P)
    wa = np.where(k % 2 == 1, 8.0 * np.cos(k * np.pi / (2 * P)) ** 2, 0.0)
    t = nmtransform._eco_open_t(x)
    g = nmtransform._eco_open_g(x)
    d = 1.0 / (y[np.newaxis, :] ** 2 + x[:, np.newaxis] ** 2)
    in_a = t >= 2.0 / P
    ra = (((d * wa).sum(axis=1) + 1.0 / P) / t - 1.0)[in_a]
    rb = d.sum(axis=1) / g - 1.0
    rms_a = np.sqrt((ra**2).mean()) if in_a.any() else np.nan
    return rms_a, np.sqrt((rb**2).mean())


@pytest.mark.parametrize("nbeads", [2, 3, 4, 8, 16, 32, 33])
@pytest.mark.parametrize("xmax", [2.0, 10.0, 50.0])
def test_eco_o_eva_structure(nbeads, xmax):
    """Checks centroid, positivity and boundedness of the open eco eigenvalues."""

    eva = nmtransform.eco_o_eva(nbeads, xmax)
    assert eva.shape == (nbeads,)
    assert eva[0] == 0.0
    assert np.all(eva[1:] > 0)
    assert np.all(np.isfinite(eva))
    # the fit caps the frequencies at twice the stiffest Trotter mode scale
    assert np.all(eva[1:] * nbeads <= 4.0 * nbeads + 1e-10)


@pytest.mark.parametrize("nbeads", [8, 16, 32, 33, 64])
@pytest.mark.parametrize("xmax", [2.0, 10.0])
def test_eco_o_eva_accuracy(nbeads, xmax):
    """Checks that open eco frequencies reproduce both the end-to-end variance
    and the open-path gyration better than the Trotter ones."""

    ea_eco, eb_eco = rms_open_errors(nmtransform.eco_o_eva(nbeads, xmax), nbeads, xmax)
    ea_tr, eb_tr = rms_open_errors(nmtransform.o_nm_eva(nbeads), nbeads, xmax)
    assert ea_eco <= ea_tr
    assert eb_eco <= eb_tr
    # well-converged regime: both errors should be dramatically smaller (at
    # the smallest bead numbers the fixed endpoint-kernel term and the exact
    # Trotter free-particle limit leave less room for improvement)
    if nbeads >= 4 * xmax and nbeads >= 16:
        assert ea_eco < ea_tr / 30.0
        assert eb_eco < eb_tr / 30.0


def test_eco_open_targets_small_x():
    """Checks the small-x expansions of the open-path targets are smooth."""

    x = np.array([1e-8, 0.1, 0.4999, 0.5001, 1.0])
    t = nmtransform._eco_open_t(x)
    g = nmtransform._eco_open_g(x)
    assert np.all(np.isfinite(t)) and np.all(np.isfinite(g))
    np.testing.assert_allclose(t[0], 1.0, rtol=1e-10)
    np.testing.assert_allclose(g[0], 1.0 / 6.0, rtol=1e-10)
    # continuity across the series/direct switchover
    assert abs(t[2] - t[3]) < 1e-4
    assert abs(g[2] - g[3]) < 1e-4


def test_eco_o_eva_free_particle_limit():
    """The Trotter open chain is exact for the free particle: the fitted
    frequencies must preserve the end-to-end variance at x -> 0 to within
    the fit accuracy."""

    nbeads, xmax = 32, 5.0
    eva = nmtransform.eco_o_eva(nbeads, xmax)
    y = eva[1:] * nbeads
    k = np.arange(1, nbeads)
    wa = np.where(k % 2 == 1, 8.0 * np.cos(k * np.pi / (2 * nbeads)) ** 2, 0.0)
    d2_0 = (wa / y**2).sum() + 1.0 / nbeads
    np.testing.assert_allclose(d2_0, 1.0, rtol=1e-3)


def test_eco_o_eva_classical_limit():
    """A single bead has no springs regardless of the fit."""

    assert np.all(nmtransform.eco_o_eva(1, 10.0) == 0.0)


def test_eco_o_eva_invalid_xmax():
    """Checks that non-positive maximum frequencies are rejected."""

    with pytest.raises(ValueError):
        nmtransform.eco_o_eva(8, 0.0)
    with pytest.raises(ValueError):
        nmtransform.eco_o_eva(8, -1.0)


def test_eco_o_eva_warm_start():
    """A fit warm-started from a nearby solution must have the same quality
    as a cold fit (the joint objective has flat valleys, so the frequencies
    themselves may differ slightly between equally good optima)."""

    nbeads = 16
    y0 = nmtransform.eco_o_eva(nbeads, 5.0)[1:] * nbeads  # previous solution
    warm = nmtransform.eco_o_eva(nbeads, 5.1, y0)
    cold = nmtransform.eco_o_eva(nbeads, 5.1)
    ea_w, eb_w = rms_open_errors(warm, nbeads, 5.1)
    ea_c, eb_c = rms_open_errors(cold, nbeads, 5.1)
    assert ea_w < 2.0 * ea_c
    assert eb_w < 2.0 * eb_c


@pytest.mark.parametrize(
    "y0",
    [
        np.array([1.0, 2.0, 3.0]),  # wrong length
        np.zeros(15),  # not strictly positive
        np.linspace(50.0, 1.0, 15),  # descending within parity classes
    ],
)
def test_eco_o_eva_bad_guess_falls_back(y0):
    """Invalid initial guesses are ignored, falling back to the two-stage fit."""

    ref = nmtransform.eco_o_eva(16, 5.0)
    np.testing.assert_allclose(nmtransform.eco_o_eva(16, 5.0, y0), ref, rtol=1e-8)


@pytest.mark.parametrize("nbeads", [2, 3, 8, 16, 33])
def test_open_spring_energy_parseval(nbeads):
    """Checks that the open normal-mode spring energy with Trotter eigenvalues
    matches the bead-difference formula of the open chain (nbeads-1 springs,
    no ring closure)."""

    rng = np.random.RandomState(27182)
    q = rng.normal(size=(nbeads, 6))
    cmat = nmtransform.mk_o_nm_matrix(nbeads)
    qnm = cmat @ q
    eva = nmtransform.o_nm_eva(nbeads)
    vnm = 0.5 * ((eva**2) @ (qnm**2))
    vbead = 0.5 * (np.diff(q, axis=0) ** 2).sum(axis=0)
    np.testing.assert_allclose(vnm, vbead, rtol=1e-10)

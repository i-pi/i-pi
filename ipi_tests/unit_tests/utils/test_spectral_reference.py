"""Spectral-analysis checks against oscillator and discrete-correlation identities."""

import numpy as np
import pytest
from ipi.utils.tools.gle import Aqp, Dqp, Cqp, gle_frequency_kernel, isra_deconvolute
from ipi.utils.tools.acf_xyz import compute_acf_xyz


@pytest.mark.parametrize("frequency", [1.0, 0.2, 3.0])
@pytest.mark.parametrize("auxiliary", [False, True])
def test_gle_covariance_solves_lyapunov_equation(frequency, auxiliary):
    ap = np.array([[1.0, 0.3], [-0.3, 2.0]]) if auxiliary else np.array([[0.7]])
    dp = ap + ap.T
    a = Aqp(frequency, ap)
    d = Dqp(frequency, dp)
    c = Cqp(frequency, a, d)
    np.testing.assert_allclose(a @ c + c @ a.T, d, atol=1e-12)
    np.testing.assert_allclose(
        c, np.diag([1 / frequency**2] + [1.0] * len(ap)), atol=1e-12
    )


def test_langevin_spectrum_matches_damped_oscillator_response():
    omega = np.linspace(0.05, 5, 100)
    # Kept off the grid's critical damping, omega = gamma/2, where the
    # oscillator's drift matrix cannot be diagonalized.
    gamma = 0.61
    dw = omega[1] - omega[0]
    kernel = gle_frequency_kernel(omega, np.array([[gamma]]), np.array([[2 * gamma]]))
    w = omega[:, None]
    w0 = omega[None, :]
    expected = 2 / np.pi * dw * gamma * w**2 / ((w0**2 - w**2) ** 2 + gamma**2 * w**2)
    np.testing.assert_allclose(kernel, expected, rtol=1e-10, atol=1e-12)


def test_deconvolution_recovers_identity_spectrum():
    omega = np.arange(20.0)
    y = np.exp(-(((omega - 10) / 3) ** 2)) + 0.01
    result, spectra, errors, curvature = isra_deconvolute(
        omega, y, np.eye(20), steps=20, stride=2
    )
    np.testing.assert_allclose(result, y, atol=1e-14)
    np.testing.assert_allclose(spectra, np.tile(y, (9, 1)), atol=1e-14)
    np.testing.assert_allclose(errors, 0, atol=1e-25)
    assert np.all(curvature >= 0)


@pytest.mark.parametrize(
    "window",
    [
        "none",
        "cosine-hanning",
        "cosine-hamming",
        "cosine-blackman",
        "triangle-bartlett",
    ],
)
@pytest.mark.parametrize("derivative", [False, True])
def test_xyz_correlation_matches_direct_circular_sum(tmp_path, window, derivative):
    data = np.column_stack(
        [
            np.sin(np.arange(16) * 2 * np.pi / 16),
            np.cos(np.arange(16) * 4 * np.pi / 16),
            np.zeros(16),
        ]
    )
    path = tmp_path / "series.xyz"
    with path.open("w") as f:
        for xyz in data:
            f.write("2\nframe\nH " + " ".join(map(str, xyz)) + "\nO 9 9 9\n")
    time, acf, err, omega, spectrum, specerr = compute_acf_xyz(
        str(path),
        4,
        block_length=16,
        labels=["H"],
        timestep=0.5,
        spectral_windowing=window,
        compute_derivative=derivative,
    )
    if derivative:
        data = np.gradient(data, axis=0) / 0.5
    expected = (
        np.array(
            [np.mean(np.sum(data * np.roll(data, k, axis=0), axis=1)) for k in range(5)]
        )
        * 0.5
        / (2 * np.pi)
    )
    np.testing.assert_allclose(acf, expected, atol=1e-12)
    np.testing.assert_allclose(time, np.arange(5) * 0.5)
    np.testing.assert_allclose(err, 0, atol=1e-14)
    np.testing.assert_allclose(specerr, 0, atol=1e-14)
    assert spectrum.shape == omega.shape

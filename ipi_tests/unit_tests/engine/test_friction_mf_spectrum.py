"""Mean force follows the physical spectrum independently of dynamical fits."""

from types import SimpleNamespace

import numpy as np
import pytest

from ipi.engine.friction import Friction


def test_mf_uses_spectrum_even_with_different_mode_fits():
    omega = np.r_[0.0, np.geomspace(1e-8, 1e5, 12000)]
    # Lorentzian spectrum: alpha(w) = theta**2*w/(rate+w).
    rate, theta = 1.5, 0.7
    spectrum = np.c_[omega, theta**2 * rate / (rate**2 + omega**2)]
    wk = np.array([0.0, 2.0, 3.0])
    friction = Friction(
        bath_mode="non-markovian",
        debug_mf_mode="on",
        Lambda=spectrum,
        Ap=np.array([[[0.0, c], [-c, 4.0]] for c in (1.0, 2.0, 3.0)]),
    )
    friction.nm = SimpleNamespace(omegak=wk)
    expected = theta**2 * wk / (rate + wk)
    np.testing.assert_allclose(friction._setup_alpha(), expected, rtol=1e-9)
    friction.Ap *= 2.0  # Including the centroid fit must not alter MF.
    np.testing.assert_allclose(friction._setup_alpha(), expected, rtol=1e-9)


def test_mf_does_not_fall_back_to_Ap():
    friction = Friction(
        bath_mode="non-markovian", Ap=np.array([[0.0, 1.0], [-1.0, 2.0]])
    )
    friction.nm = SimpleNamespace(omegak=np.array([0.0, 2.0]))
    with pytest.raises(ValueError, match="MF requires Lambda"):
        friction._setup_alpha()
    friction.debug_mf_mode = "off"
    np.testing.assert_array_equal(friction._setup_alpha(), [0.0, 0.0])


@pytest.mark.parametrize(
    "spectrum",
    [
        [0.0, 1.0],
        [[0.0, 1.0]],
        [[0.0, 1.0], [0.0, 2.0]],
        [[-1.0, 1.0], [1.0, 1.0]],
        [[0.0, -1.0], [1.0, 1.0]],
        [[0.0, 1.0], [1.0, np.nan]],
    ],
)
def test_invalid_mf_spectrum_is_rejected(spectrum):
    friction = Friction(Lambda=np.array(spectrum))
    friction.nm = SimpleNamespace(omegak=np.array([0.0, 2.0]))
    with pytest.raises(ValueError, match="Lambda"):
        friction._setup_alpha()


@pytest.mark.parametrize(
    "values", [[float("nan"), 1.0], [0.0, -1.0], [1.0, 1.0], [0.0, float("inf")]]
)
def test_invalid_explicit_alpha_rejected(values):
    friction = Friction(bath_mode="none", debug_alpha_input=np.c_[[0.0, 1.0], values])
    friction.nm = SimpleNamespace(omegak=np.array([0.0, 1.0]))
    with pytest.raises(ValueError, match="debug_alpha_input"):
        friction._setup_alpha()


def test_positive_spectrum_has_positive_alpha_without_boundary_extrapolation():
    from ipi.engine.friction import get_alpha_numeric

    # The old cubic extrapolation gave a negative alpha for this positive input.
    with pytest.raises(ValueError, match="starting at zero"):
        get_alpha_numeric([0.0, 1.0, 0.0, 0.0], [1.0, 2.0, 3.0, 4.0], [0.0, 1.0])
    alpha = get_alpha_numeric(
        [0.0, 0.0, 1.0, 0.0, 0.0], [0.0, 1.0, 2.0, 3.0, 4.0], [0.0, 0.1, 1.0, 10.0]
    )
    assert alpha[0] == 0
    assert np.all(alpha[1:] > 0)


def test_ohmic_centroid_limit():
    from ipi.utils.frictiontools import get_alpha_ohmic
    from scipy.integrate import quad

    frequencies = np.array([0.0, 0.01, 1.0, 10.0])
    alpha = get_alpha_ohmic(frequencies, 1.0, 2.0)
    assert alpha[0] == 0
    for frequency, actual in zip(frequencies[1:], alpha[1:]):
        expected = (
            2
            / np.pi
            * quad(
                lambda w: 2 * np.exp(-w) * frequency**2 / (w * w + frequency**2),
                0,
                np.inf,
            )[0]
        )
        np.testing.assert_allclose(actual, expected, rtol=1e-10)


@pytest.mark.parametrize("quantity", ["Lambda", "J"])
def test_spectrum_cli_convention(tmp_path, quantity):
    import subprocess
    import sys
    from ipi.utils import frictiontools

    path = tmp_path / "spectrum.dat"
    args = [
        sys.executable,
        frictiontools.__file__,
        "--omega-cut",
        "2",
        "--eta",
        "3",
        "--points",
        "51",
        str(path),
    ]
    if quantity == "J":
        args += ["--quantity", "J"]
    subprocess.run(args, check=True, capture_output=True)
    omega, values = np.loadtxt(path).T
    assert omega[0] == 0 and omega[-1] == 40
    expected = 3 * np.exp(-omega / 2)
    if quantity == "J":
        expected *= omega
    np.testing.assert_allclose(values, expected)

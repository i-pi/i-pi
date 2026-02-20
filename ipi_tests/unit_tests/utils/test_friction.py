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

import json
import pytest
import numpy as np
from ipi.utils.frictiontools import expohmic_Lambda, get_alpha_ohmic

from ipi.engine.friction import Friction, get_alpha_numeric
from ipi.pes.doublewell_with_friction import DoubleWell_with_friction_driver


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

    omega = np.linspace(0, 10 * omega_cutoff, 100000)
    Lambda = expohmic_Lambda(omega, eta, omega_cutoff)
    # print(Lambda)
    alpha = get_alpha_numeric(Lambda, omega, omegak)
    assert np.allclose(
        alpha,
        alphak,
        rtol=1e-6,
    )
    print(alpha, alphak)


class _IdentityTransform:
    @staticmethod
    def nm2b(arr: np.ndarray) -> np.ndarray:
        return np.array(arr, copy=True)

    @staticmethod
    def b2nm(arr: np.ndarray) -> np.ndarray:
        return np.array(arr, copy=True)


class _DummyPRNG:
    def gvec(self, shape):
        # Deterministic thermostat/friction noise for reproducible tests.
        return np.zeros(shape, dtype=float)


class _DummyForces:
    def __init__(self, extras: dict):
        self.extras = extras


class _DummyEnsemble:
    def __init__(self, temp: float, forces: _DummyForces):
        self.temp = temp
        self.forces = forces


class _DummyBeads:
    def __init__(self, q: np.ndarray, p: np.ndarray, m3: np.ndarray):
        self.q = np.array(q, dtype=float, copy=True)
        self.p = np.array(p, dtype=float, copy=True)
        self.m3 = np.array(m3, dtype=float, copy=True)
        self.nbeads = int(self.q.shape[0])
        self.natoms = int(self.q.shape[1] // 3)


class _DummyNM:
    def __init__(self, beads: _DummyBeads):
        self.pnm = np.array(beads.p, dtype=float, copy=True)
        self.qnm = np.array(beads.q, dtype=float, copy=True)
        self.dynm3 = np.array(beads.m3, dtype=float, copy=True)
        self.transform = _IdentityTransform()


def _build_doublewell_diffusion(qx: float) -> np.ndarray:
    driver = DoubleWell_with_friction_driver(
        w_b=500.0,
        v0=2085.0,
        m=1837.36223469,
        eta0=1.5,
        eps1=-1.0,
        eps2=0.0,
        delta=0.0,
        deltaQ=1.0,
    )
    _, _, _, extras_json = driver.compute_structure(
        np.zeros((3, 3), dtype=float),
        np.array([[qx, 0.0, 0.0]], dtype=float),
    )
    extras = json.loads(extras_json)
    # Canonical shape expected by Friction: (nbeads, nbath, ndof_subset/full)
    return np.asarray(extras["diffusion_coefficient"], dtype=float)[:, np.newaxis, :]


@pytest.mark.parametrize("bath_mode", ["none", "markovian"])
@pytest.mark.parametrize("variable_friction", [False, True])
def test_doublewell_friction_bath_modes(
    bath_mode: str, variable_friction: bool
) -> None:
    """Covers friction step for DW friction model with/without variable friction.

    Cases:
      - bath_mode='none'      -> momentum unchanged
      - bath_mode='markovian' -> momentum updated
    """
    q = np.array([[0.2, 0.0, 0.0]], dtype=float)
    p0 = np.array([[1.0, 0.5, -0.3]], dtype=float)
    m3 = np.ones_like(p0) * 1837.36223469

    beads = _DummyBeads(q=q, p=p0, m3=m3)
    nm = _DummyNM(beads)

    sigma = _build_doublewell_diffusion(qx=float(q[0, 0]))
    forces = _DummyForces(extras={"diffusion_coefficient": sigma})
    ensemble = _DummyEnsemble(temp=300.0, forces=forces)

    friction = Friction(
        variable_friction=variable_friction,
        bath_mode=bath_mode,
        mf_mode="none",
        sigma_static=1.2,
        sigma_key="diffusion_coefficient",
        friction_atoms=np.array([0], dtype=int),
    )
    friction.beads = beads
    friction.nm = nm
    friction.ensemble = ensemble
    friction.forces = forces
    friction.prng = _DummyPRNG()
    friction._setup_friction_atoms()

    p_before = np.array(beads.p, copy=True)
    friction.step(0.1)
    p_after = np.array(beads.p, copy=True)

    if bath_mode == "none":
        np.testing.assert_allclose(p_after, p_before)
        return

    # Markovian: must change at least one momentum component.
    assert np.linalg.norm(p_after - p_before) > 0.0

    if variable_friction:
        # For DW friction driver, diffusion vector has only x component.
        np.testing.assert_allclose(p_after[:, 1:], p_before[:, 1:])
        assert np.linalg.norm(p_after[:, 0] - p_before[:, 0]) > 0.0


def test_nested_sigma_payload_concatenates_channels() -> None:
    q = np.array([[0.0, 0.0, 0.0]], dtype=float)
    p = np.zeros_like(q)
    m3 = np.ones_like(q)
    beads = _DummyBeads(q=q, p=p, m3=m3)
    nm = _DummyNM(beads)

    # Two independent sigma channels represented as 1x3 matrices.
    # They must be concatenated along bath rows, not summed.
    extras_sigma = {
        "equ": {"1": [[1.0, 2.0, 3.0]]},
        "inv": {"1": [[4.0, 5.0, 6.0]]},
    }
    forces = _DummyForces(extras={"sigma": extras_sigma})
    ensemble = _DummyEnsemble(temp=300.0, forces=forces)

    friction = Friction(
        variable_friction=True,
        bath_mode="markovian",
        mf_mode="none",
        sigma_key="sigma",
        friction_atoms=np.array([0], dtype=int),
    )
    friction.beads = beads
    friction.nm = nm
    friction.ensemble = ensemble
    friction.forces = forces
    friction.prng = _DummyPRNG()
    friction._setup_friction_atoms()

    sigma = friction._get_sigma()
    expected = np.array([[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]], dtype=float)
    np.testing.assert_allclose(sigma, expected)


def test_nested_sigma_single_dict_raises_for_multiple_beads() -> None:
    q = np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]], dtype=float)
    p = np.zeros_like(q)
    m3 = np.ones_like(q)
    beads = _DummyBeads(q=q, p=p, m3=m3)
    nm = _DummyNM(beads)

    extras_sigma = {
        "equ": {"1": [[1.0, 2.0, 3.0]]},
        "inv": {"1": [[4.0, 5.0, 6.0]]},
    }
    forces = _DummyForces(extras={"sigma": extras_sigma})
    ensemble = _DummyEnsemble(temp=300.0, forces=forces)

    friction = Friction(
        variable_friction=True,
        bath_mode="markovian",
        mf_mode="none",
        sigma_key="sigma",
        friction_atoms=np.array([0], dtype=int),
    )
    friction.beads = beads
    friction.nm = nm
    friction.ensemble = ensemble
    friction.forces = forces
    friction.prng = _DummyPRNG()
    friction._setup_friction_atoms()

    with pytest.raises(ValueError, match="single dict payload"):
        friction._get_sigma()


def test_nested_sigma_list_of_dicts_length_mismatch_raises() -> None:
    q = np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]], dtype=float)
    p = np.zeros_like(q)
    m3 = np.ones_like(q)
    beads = _DummyBeads(q=q, p=p, m3=m3)
    nm = _DummyNM(beads)

    # nbeads=2 but provide only one bead payload.
    extras_sigma = [
        {
            "equ": {"1": [[1.0, 2.0, 3.0]]},
            "inv": {"1": [[4.0, 5.0, 6.0]]},
        }
    ]
    forces = _DummyForces(extras={"sigma": extras_sigma})
    ensemble = _DummyEnsemble(temp=300.0, forces=forces)

    friction = Friction(
        variable_friction=True,
        bath_mode="markovian",
        mf_mode="none",
        sigma_key="sigma",
        friction_atoms=np.array([0], dtype=int),
    )
    friction.beads = beads
    friction.nm = nm
    friction.ensemble = ensemble
    friction.forces = forces
    friction.prng = _DummyPRNG()
    friction._setup_friction_atoms()

    with pytest.raises(ValueError, match="list-of-dicts length"):
        friction._get_sigma()


def test_gamma_from_sigma_variable_is_sigma_t_sigma() -> None:
    q = np.array([[0.0, 0.0, 0.0]], dtype=float)
    p = np.zeros_like(q)
    m3 = np.ones_like(q)
    beads = _DummyBeads(q=q, p=p, m3=m3)
    nm = _DummyNM(beads)

    # (nbeads=1, nbath=2, ndof=3)
    sigma = np.array([[[1.0, 2.0, 0.0], [0.0, 1.0, 3.0]]], dtype=float)
    forces = _DummyForces(extras={"sigma": sigma})
    ensemble = _DummyEnsemble(temp=300.0, forces=forces)

    friction = Friction(
        variable_friction=True,
        bath_mode="markovian",
        mf_mode="none",
        sigma_key="sigma",
        friction_atoms=np.array([0], dtype=int),
    )
    friction.beads = beads
    friction.nm = nm
    friction.ensemble = ensemble
    friction.forces = forces
    friction.prng = _DummyPRNG()
    friction._setup_friction_atoms()

    gamma = friction.gamma
    expected = np.array([sigma[0].T @ sigma[0]], dtype=float)
    np.testing.assert_allclose(gamma, expected)


def test_gamma_from_sigma_static_is_scalar_square() -> None:
    friction = Friction(
        variable_friction=False,
        bath_mode="markovian",
        mf_mode="none",
        sigma_static=1.7,
    )
    np.testing.assert_allclose(friction.gamma, 1.7 * 1.7)

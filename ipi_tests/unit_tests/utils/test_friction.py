"""Electronic-friction algebra, canonical driver payloads and spectral quadrature."""

import json
import pytest
import numpy as np
from ipi.utils.depend import depend_value
from ipi.utils.frictiontools import (
    expohmic_Lambda,
    get_alpha_ohmic,
)

from ipi.engine.friction import Friction, FrictionGLE, get_alpha_numeric
from ipi.inputs.friction import InputFriction
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
    alpha = get_alpha_numeric(Lambda, omega, omegak)
    assert np.allclose(
        alpha,
        alphak,
        rtol=1e-6,
    )


def test_debug_mf_mode_defaults_on_and_maps_legacy_none_to_off() -> None:
    assert Friction().debug_mf_mode == "on"
    assert Friction(debug_mf_mode="off").debug_mf_mode == "off"
    assert Friction(debug_mf_mode="none").debug_mf_mode == "off"
    with pytest.raises(ValueError, match="debug_mf_mode"):
        Friction(debug_mf_mode="linear")


def test_canonical_friction_input_keys_are_the_defaults() -> None:
    friction = InputFriction().fetch()
    assert friction.coupling_jacobian_key == "friction_coupling_jacobian"
    assert friction.gamma_key == "friction_gamma"
    assert friction.friction_meta_key == "friction_meta"


def test_mode_resolved_Ap_input_round_trip() -> None:
    Ap = np.array(
        [
            [[0.0, 0.5], [-0.5, 1.5]],
            [[0.0, 0.8], [-0.8, 2.0]],
        ],
        dtype=float,
    )
    source = Friction(Ap=Ap, bath_mode="non-markovian")
    input_friction = InputFriction()
    input_friction.store(source)
    restored = input_friction.fetch()

    np.testing.assert_allclose(restored.Ap, Ap)
    assert restored.Ap.shape == (2, 2, 2)


def test_debug_alpha_overrides_Ap_and_Lambda() -> None:
    q = np.zeros((3, 3), dtype=float)
    p = np.zeros_like(q)
    m3 = np.ones_like(q)
    beads = _DummyBeads(q=q, p=p, m3=m3)
    nm = _DummyNM(beads)
    nm.omegak = np.array([0.0, 2.0, 3.0], dtype=float)
    forces = _DummyForces(extras={})
    ensemble = _DummyEnsemble(temp=300.0, forces=forces)

    friction = Friction(
        variable_friction=False,
        bath_mode="non-markovian",
        debug_mf_mode="on",
        Ap=np.array([[0.0, 0.5], [-0.5, 1.5]], dtype=float),
        # Debug coefficients explicitly override the spectrum; Ap is dynamics only.
        Lambda=np.array([[0.0, 1.0], [10.0, 0.1]]),
        debug_alpha_input=np.array([[0.0, 0.0], [2.0, 9.0], [3.0, 9.0]], dtype=float),
    )
    friction.beads = beads
    friction.nm = nm
    friction.ensemble = ensemble
    friction.forces = forces
    friction.prng = _DummyPRNG()

    alpha = friction._setup_alpha()
    expected = np.array([0.0, 9.0, 9.0], dtype=float)
    np.testing.assert_allclose(alpha, expected)


def test_debug_mf_mode_off_disables_mf_kick() -> None:
    q = np.array([[1.0, 0.0, 0.0]], dtype=float)
    p0 = np.array([[0.2, -0.1, 0.3]], dtype=float)
    m3 = np.ones_like(q)
    beads = _DummyBeads(q=q, p=p0, m3=m3)
    nm = _DummyNM(beads)
    nm.omegak = np.array([1.0], dtype=float)
    forces = _DummyForces(extras={})
    ensemble = _DummyEnsemble(temp=300.0, forces=forces)

    friction = Friction(
        variable_friction=False,
        bath_mode="none",
        debug_mf_mode="off",
        sigma_static=2.0,
        debug_alpha_input=np.array([[1.0, 5.0]], dtype=float),
    )
    friction.beads = beads
    friction.nm = nm
    friction.ensemble = ensemble
    friction.forces = forces
    friction.prng = _DummyPRNG()
    friction.alpha = friction._setup_alpha()

    assert friction.get_energy_mf() == 0.0
    friction.step(0.25)
    np.testing.assert_allclose(beads.p, p0)


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
    @property
    def pnm_(self):
        return self.pnm

    def __init__(self, beads: _DummyBeads):
        self.pnm = np.array(beads.p, dtype=float, copy=True)
        self.qnm = np.array(beads.q, dtype=float, copy=True)
        self.dynm3 = np.array(beads.m3, dtype=float, copy=True)
        self.omegak = np.zeros(beads.nbeads, dtype=float)
        self.transform = _IdentityTransform()


def _build_doublewell_diffusion(qx: float) -> np.ndarray:
    driver = DoubleWell_with_friction_driver(
        w_b=500.0,
        v0=2085.0,
        m=1837.36223469,
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


def _canonical_extras(
    jacobian: np.ndarray,
    active_atoms: list[int] | None = None,
    labels: list[str] | None = None,
) -> dict:
    """Build the model-independent extras contract used by friction tests."""
    jacobian = np.asarray(jacobian, dtype=float)
    if active_atoms is None:
        active_atoms = list(range(jacobian.shape[-1] // 3))
    if labels is None:
        labels = [f"channel_{i + 1}" for i in range(jacobian.shape[1])]
    gamma = np.einsum("bci,bcj->bij", jacobian, jacobian)
    meta = {
        "schema": "ipi_friction_v1",
        "active_atoms": active_atoms,
        "channel_labels": labels,
        "jacobian_units": "ipi_internal_mass_scaled",
        "gamma_units": "ipi_internal_mass_scaled",
    }
    return {
        "friction_coupling_jacobian": jacobian,
        "friction_gamma": gamma,
        "friction_meta": [dict(meta) for _ in range(jacobian.shape[0])],
    }


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

    jacobian = _build_doublewell_diffusion(qx=float(q[0, 0]))
    forces = _DummyForces(extras=_canonical_extras(jacobian))
    ensemble = _DummyEnsemble(temp=300.0, forces=forces)

    friction = Friction(
        variable_friction=variable_friction,
        bath_mode=bath_mode,
        debug_mf_mode="none",
        sigma_static=1.2,
    )
    friction.beads = beads
    friction.nm = nm
    friction.ensemble = ensemble
    friction.forces = forces
    friction.prng = _DummyPRNG()
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


def test_canonical_jacobian_preserves_channels_and_embeds_active_dofs() -> None:
    q = np.zeros((1, 6), dtype=float)
    p = np.zeros_like(q)
    beads = _DummyBeads(q=q, p=p, m3=np.ones_like(q))
    reduced = np.array([[[1.0, 2.0, 3.0], [10.0, 20.0, 30.0]]])
    forces = _DummyForces(extras=_canonical_extras(reduced, active_atoms=[1]))
    friction = Friction(
        variable_friction=True, bath_mode="markovian", debug_mf_mode="off"
    )
    friction.beads = beads
    friction.nm = _DummyNM(beads)
    friction.forces = forces
    friction.ensemble = _DummyEnsemble(300.0, forces)
    friction.prng = _DummyPRNG()

    expected = np.zeros((1, 2, 6))
    expected[:, :, 3:] = reduced
    np.testing.assert_allclose(friction._get_coupling_jacobian(), expected)
    assert friction._channel_labels == ["channel_1", "channel_2"]


def test_canonical_metadata_mismatch_fails_loudly() -> None:
    jacobian = np.ones((2, 1, 3))
    extras = _canonical_extras(jacobian)
    extras["friction_meta"][1]["channel_labels"] = ["different"]
    beads = _DummyBeads(q=np.zeros((2, 3)), p=np.zeros((2, 3)), m3=np.ones((2, 3)))
    friction = Friction(
        variable_friction=True, bath_mode="markovian", debug_mf_mode="off"
    )
    friction.beads = beads
    friction.nm = _DummyNM(beads)
    friction.forces = _DummyForces(extras)
    with pytest.raises(ValueError, match="differs across beads"):
        friction._get_coupling_jacobian()


def test_canonical_channel_contract_cannot_change_between_steps() -> None:
    jacobian = np.ones((1, 1, 3))
    extras = _canonical_extras(jacobian, labels=["first"])
    beads = _DummyBeads(q=np.zeros((1, 3)), p=np.zeros((1, 3)), m3=np.ones((1, 3)))
    friction = Friction(
        variable_friction=True, bath_mode="markovian", debug_mf_mode="off"
    )
    friction.beads = beads
    friction.nm = _DummyNM(beads)
    friction.forces = _DummyForces(extras)
    friction._get_coupling_jacobian()
    extras["friction_meta"][0]["channel_labels"] = ["changed"]
    with pytest.raises(ValueError, match="changed during the simulation"):
        friction._get_coupling_jacobian()


def test_driver_gamma_is_validated_against_jacobian() -> None:
    jacobian = np.array([[[1.0, 2.0, 3.0]]])
    extras = _canonical_extras(jacobian)
    extras["friction_gamma"][0, 0, 0] += 1.0
    beads = _DummyBeads(q=np.zeros((1, 3)), p=np.zeros((1, 3)), m3=np.ones((1, 3)))
    friction = Friction(
        variable_friction=True, bath_mode="markovian", debug_mf_mode="off"
    )
    friction.beads = beads
    friction.nm = _DummyNM(beads)
    friction.forces = _DummyForces(extras)
    with pytest.raises(ValueError, match="inconsistent with coupling Jacobian"):
        _ = friction.gamma


def test_gamma_from_static_coupling_is_scalar_square() -> None:
    friction = Friction(
        variable_friction=False,
        bath_mode="markovian",
        debug_mf_mode="off",
        sigma_static=1.7,
    )
    np.testing.assert_allclose(friction.gamma, 1.7 * 1.7)


def test_markovian_bath_is_unified_friction_gle() -> None:
    friction = Friction(
        variable_friction=False,
        bath_mode="markovian",
        debug_mf_mode="none",
    )
    bath = friction._build_bath()
    assert isinstance(bath, FrictionGLE)


def test_step_builds_friction_gle_lazily() -> None:
    q = np.array([[0.2, 0.0, 0.0]], dtype=float)
    p0 = np.array([[1.0, 0.5, -0.3]], dtype=float)
    m3 = np.ones_like(p0) * 1837.36223469
    beads = _DummyBeads(q=q, p=p0, m3=m3)
    nm = _DummyNM(beads)
    forces = _DummyForces(extras={})
    ensemble = _DummyEnsemble(temp=300.0, forces=forces)

    friction = Friction(
        variable_friction=False,
        bath_mode="markovian",
        debug_mf_mode="none",
        sigma_static=1.2,
    )
    friction.beads = beads
    friction.nm = nm
    friction.ensemble = ensemble
    friction.forces = forces
    friction.prng = _DummyPRNG()

    assert friction.bath is None
    friction.step(0.1)
    assert isinstance(friction.bath, FrictionGLE)


def test_non_markovian_friction_gle_has_auxiliary_scaffold() -> None:
    q = np.array([[0.2, 0.0, 0.0]], dtype=float)
    p0 = np.array([[1.0, 0.5, -0.3]], dtype=float)
    m3 = np.ones_like(p0)
    beads = _DummyBeads(q=q, p=p0, m3=m3)
    nm = _DummyNM(beads)
    nm.omegak = np.array([0.0], dtype=float)
    forces = _DummyForces(extras={})
    ensemble = _DummyEnsemble(temp=300.0, forces=forces)

    friction = Friction(
        variable_friction=False,
        bath_mode="non-markovian",
        debug_mf_mode="none",
        Ap=np.array([[0.0, 0.5], [-0.5, 1.5]], dtype=float),
    )
    friction.beads = beads
    friction.nm = nm
    friction.ensemble = ensemble
    friction.forces = forces
    friction.prng = _DummyPRNG()
    friction._ensure_bath_bound()

    assert isinstance(friction.bath, FrictionGLE)
    assert friction.bath.state_shape() == (1, 1, 3)
    assert friction.bath.s.shape == (1, 1, 3)
    np.testing.assert_allclose(friction.bath.theta, np.array([[0.5]]))
    np.testing.assert_allclose(friction.bath.A_aux, np.array([[[1.5]]]))


def test_multibead_non_markovian_rejects_centroid_only_Ap() -> None:
    q = np.zeros((2, 3), dtype=float)
    beads = _DummyBeads(q=q, p=np.zeros_like(q), m3=np.ones_like(q))
    friction = Friction(
        variable_friction=False,
        bath_mode="non-markovian",
        debug_mf_mode="off",
        Ap=np.array([[0.0, 0.5], [-0.5, 1.5]], dtype=float),
    )
    friction.beads = beads

    with pytest.raises(ValueError, match="cannot be broadcast"):
        friction._validate_Ap()


def test_multibead_non_markovian_rejects_wrong_mode_count() -> None:
    q = np.zeros((3, 3), dtype=float)
    beads = _DummyBeads(q=q, p=np.zeros_like(q), m3=np.ones_like(q))
    friction = Friction(
        variable_friction=False,
        bath_mode="non-markovian",
        debug_mf_mode="off",
        Ap=np.array(
            [
                [[0.0, 0.5], [-0.5, 1.5]],
                [[0.0, 0.8], [-0.8, 2.0]],
            ],
            dtype=float,
        ),
    )
    friction.beads = beads

    with pytest.raises(ValueError, match="has 2 modes.*3 beads"):
        friction._validate_Ap()


def test_non_markovian_rejects_nondissipative_auxiliary_block() -> None:
    q = np.zeros((1, 3), dtype=float)
    beads = _DummyBeads(q=q, p=np.zeros_like(q), m3=np.ones_like(q))
    friction = Friction(
        variable_friction=False,
        bath_mode="non-markovian",
        debug_mf_mode="off",
        Ap=np.array(
            [
                [0.0, 0.5, 0.2],
                [-0.5, 1.0, 4.0],
                [-0.2, 0.0, 1.0],
            ],
            dtype=float,
        ),
    )
    friction.beads = beads

    with pytest.raises(ValueError, match="symmetric part.*not positive semidefinite"):
        friction._validate_Ap()


def test_mode_resolved_Ap_applies_distinct_auxiliary_dynamics() -> None:
    q = np.zeros((2, 3), dtype=float)
    p0 = np.array([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]], dtype=float)
    beads = _DummyBeads(q=q, p=p0, m3=np.ones_like(q))
    nm = _DummyNM(beads)
    nm.omegak = np.array([0.0, 2.0], dtype=float)
    Ap = np.array(
        [
            [[0.0, 0.5], [-0.5, 1.5]],
            [[0.0, 0.8], [-0.8, 2.0]],
        ],
        dtype=float,
    )
    forces = _DummyForces(extras={})
    friction = Friction(
        variable_friction=False,
        bath_mode="non-markovian",
        debug_mf_mode="off",
        sigma_static=1.0,
        Ap=Ap,
    )
    friction.beads = beads
    friction.nm = nm
    friction.ensemble = _DummyEnsemble(temp=300.0, forces=forces)
    friction.forces = forces
    friction.prng = _DummyPRNG()
    friction._validate_Ap()
    friction._ensure_bath_bound()

    friction.bath.s[:] = np.array([[[0.2, 0.0, 0.0]], [[0.3, 0.0, 0.0]]])
    friction.bath.bp_f_step(0.1)
    expected_p = p0.copy()
    expected_p[0, 0] -= 0.1 * 0.5 * 0.2
    expected_p[1, 0] -= 0.1 * 0.8 * 0.3
    np.testing.assert_allclose(friction.nm.pnm, expected_p)

    s_before = friction.bath.s.copy()
    friction.bath.bs_step(0.2)
    expected_s = s_before.copy()
    expected_s[0, 0, :] += 0.2 * 0.5 * expected_p[0, :]
    expected_s[1, 0, :] += 0.2 * 0.8 * expected_p[1, :]
    np.testing.assert_allclose(friction.bath.s, expected_s)

    friction.bath.os_step(0.3)
    expected_s[0] *= np.exp(-1.5 * 0.3)
    expected_s[1] *= np.exp(-2.0 * 0.3)
    np.testing.assert_allclose(friction.bath.s, expected_s)


def test_mode_resolved_Ap_uses_matching_variable_coupling_mode() -> None:
    q = np.zeros((2, 3), dtype=float)
    p0 = np.array([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]], dtype=float)
    beads = _DummyBeads(q=q, p=p0, m3=np.ones_like(q))
    nm = _DummyNM(beads)
    nm.omegak = np.array([0.0, 2.0], dtype=float)
    Ap = np.array(
        [
            [[0.0, 0.5], [-0.5, 1.5]],
            [[0.0, 0.8], [-0.8, 2.0]],
        ],
        dtype=float,
    )
    # Equal bead-space Jacobians transform to delta[r,n], so each auxiliary
    # mode must couple only to the normal mode with the same index.
    jacobian = np.array([[[1.0, 0.0, 0.0]], [[1.0, 0.0, 0.0]]])
    forces = _DummyForces(extras=_canonical_extras(jacobian))
    friction = Friction(
        variable_friction=True,
        bath_mode="non-markovian",
        debug_mf_mode="off",
        Ap=Ap,
    )
    friction.beads = beads
    friction.nm = nm
    friction.ensemble = _DummyEnsemble(temp=300.0, forces=forces)
    friction.forces = forces
    friction.prng = _DummyPRNG()
    friction._validate_Ap()
    friction._ensure_bath_bound()
    friction.bath._get_non_markovian_jacobian()

    friction.bath.s[:] = np.array([[[0.2]], [[0.3]]])
    friction.bath.bp_f_step(0.1)
    expected_p = p0.copy()
    expected_p[0, 0] -= 0.1 * 0.5 * 0.2
    expected_p[1, 0] -= 0.1 * 0.8 * 0.3
    np.testing.assert_allclose(friction.nm.pnm, expected_p)

    s_before = friction.bath.s.copy()
    friction.bath.bs_step(0.2)
    expected_s = s_before.copy()
    expected_s[0, 0, 0] += 0.2 * 0.5 * expected_p[0, 0]
    expected_s[1, 0, 0] += 0.2 * 0.8 * expected_p[1, 0]
    np.testing.assert_allclose(friction.bath.s, expected_s)


def test_non_markovian_step_runs_with_explicit_Ap() -> None:
    q = np.array([[0.2, 0.0, 0.0]], dtype=float)
    p0 = np.array([[1.0, 0.5, -0.3]], dtype=float)
    m3 = np.ones_like(p0)
    beads = _DummyBeads(q=q, p=p0, m3=m3)
    nm = _DummyNM(beads)
    nm.omegak = np.array([0.0], dtype=float)
    forces = _DummyForces(extras={})
    ensemble = _DummyEnsemble(temp=300.0, forces=forces)

    friction = Friction(
        variable_friction=False,
        bath_mode="non-markovian",
        debug_mf_mode="none",
        Ap=np.array([[0.0, 0.5], [-0.5, 1.5]], dtype=float),
    )
    friction.beads = beads
    friction.nm = nm
    friction.ensemble = ensemble
    friction.forces = forces
    friction.prng = _DummyPRNG()

    friction._ensure_bath_bound()
    friction.bath.s[:] = np.array([[[0.2, -0.1, 0.3]]])

    friction.step(0.1)

    assert np.all(np.isfinite(friction.bath.s))
    assert np.all(np.isfinite(friction.nm.pnm))
    assert np.linalg.norm(friction.bath.s) > 0.0
    assert np.linalg.norm(friction.nm.pnm - p0) > 0.0


def test_non_markovian_position_independent_coupling_steps_are_deterministic() -> None:
    q = np.array([[0.0, 0.0, 0.0]], dtype=float)
    p0 = np.array([[1.0, 2.0, -1.0]], dtype=float)
    m3 = np.ones_like(p0)
    beads = _DummyBeads(q=q, p=p0, m3=m3)
    nm = _DummyNM(beads)
    nm.omegak = np.array([0.0], dtype=float)
    forces = _DummyForces(extras={})
    ensemble = _DummyEnsemble(temp=300.0, forces=forces)

    friction = Friction(
        variable_friction=False,
        bath_mode="non-markovian",
        debug_mf_mode="none",
        sigma_static=2.0,
        Ap=np.array([[0.0, 0.5], [-0.5, 1.5]], dtype=float),
    )
    friction.beads = beads
    friction.nm = nm
    friction.ensemble = ensemble
    friction.forces = forces
    friction.prng = _DummyPRNG()
    friction._ensure_bath_bound()

    friction.bath.s[:] = np.array([[[0.2, -0.4, 0.6]]])
    friction.bath.bp_f_step(0.1)
    expected_p = p0 - 0.1 * 2.0 * 0.5 * np.array([[0.2, -0.4, 0.6]])
    np.testing.assert_allclose(friction.nm.pnm, expected_p)

    s_before = friction.bath.s.copy()
    friction.bath.bs_step(0.2)
    expected_s = s_before + 0.2 * 0.5 * 2.0 * expected_p[:, np.newaxis, :]
    np.testing.assert_allclose(friction.bath.s, expected_s)

    friction.bath.os_step(0.3)
    np.testing.assert_allclose(
        friction.bath.s,
        expected_s * np.exp(-1.5 * 0.3),
    )


def test_non_markovian_variable_friction_gle_tracks_bath_channels() -> None:
    q = np.array([[0.0, 0.0, 0.0]], dtype=float)
    p0 = np.array([[1.0, 2.0, -1.0]], dtype=float)
    m3 = np.ones_like(p0)
    beads = _DummyBeads(q=q, p=p0, m3=m3)
    nm = _DummyNM(beads)
    nm.omegak = np.array([0.0], dtype=float)
    jacobian = np.array([[[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]]], dtype=float)
    forces = _DummyForces(extras=_canonical_extras(jacobian))
    ensemble = _DummyEnsemble(temp=300.0, forces=forces)

    friction = Friction(
        variable_friction=True,
        bath_mode="non-markovian",
        debug_mf_mode="none",
        Ap=np.array([[0.0, 0.5], [-0.5, 1.5]], dtype=float),
    )
    friction.beads = beads
    friction.nm = nm
    friction.ensemble = ensemble
    friction.forces = forces
    friction.prng = _DummyPRNG()
    friction._ensure_bath_bound()

    assert isinstance(friction.bath, FrictionGLE)
    friction.bath._get_non_markovian_jacobian()
    assert friction.bath.state_shape() == (1, 1, 2)
    assert friction.bath.s.shape == (1, 1, 2)


def test_non_markovian_variable_coupling_steps_are_deterministic() -> None:
    q = np.array([[0.0, 0.0, 0.0]], dtype=float)
    p0 = np.array([[1.0, 2.0, -1.0]], dtype=float)
    m3 = np.ones_like(p0)
    beads = _DummyBeads(q=q, p=p0, m3=m3)
    nm = _DummyNM(beads)
    nm.omegak = np.array([0.0], dtype=float)
    jacobian = np.array([[[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]]], dtype=float)
    forces = _DummyForces(extras=_canonical_extras(jacobian))
    ensemble = _DummyEnsemble(temp=300.0, forces=forces)

    friction = Friction(
        variable_friction=True,
        bath_mode="non-markovian",
        debug_mf_mode="none",
        Ap=np.array([[0.0, 0.5], [-0.5, 1.5]], dtype=float),
    )
    friction.beads = beads
    friction.nm = nm
    friction.ensemble = ensemble
    friction.forces = forces
    friction.prng = _DummyPRNG()
    friction._ensure_bath_bound()
    friction.bath._get_non_markovian_jacobian()

    friction.bath.s[:] = np.array([[[0.2, -0.4]]])
    friction.bath.bp_f_step(0.1)
    expected_p = p0 - 0.1 * np.array([[0.1, -0.4, 0.0]])
    np.testing.assert_allclose(friction.nm.pnm, expected_p)

    s_before = friction.bath.s.copy()
    friction.bath.bs_step(0.2)
    expected_s = s_before + 0.2 * 0.5 * np.array([[[0.99, 4.08]]])
    np.testing.assert_allclose(friction.bath.s, expected_s)

    friction.bath.os_step(0.3)
    np.testing.assert_allclose(
        friction.bath.s,
        expected_s * np.exp(-1.5 * 0.3),
    )


def test_non_markovian_variable_step_runs_with_explicit_Ap() -> None:
    q = np.array([[0.0, 0.0, 0.0]], dtype=float)
    p0 = np.array([[1.0, 2.0, -1.0]], dtype=float)
    m3 = np.ones_like(p0)
    beads = _DummyBeads(q=q, p=p0, m3=m3)
    nm = _DummyNM(beads)
    nm.omegak = np.array([0.0], dtype=float)
    jacobian = np.array([[[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]]], dtype=float)
    forces = _DummyForces(extras=_canonical_extras(jacobian))
    ensemble = _DummyEnsemble(temp=300.0, forces=forces)

    friction = Friction(
        variable_friction=True,
        bath_mode="non-markovian",
        debug_mf_mode="none",
        Ap=np.array([[0.0, 0.5], [-0.5, 1.5]], dtype=float),
    )
    friction.beads = beads
    friction.nm = nm
    friction.ensemble = ensemble
    friction.forces = forces
    friction.prng = _DummyPRNG()
    friction._ensure_bath_bound()
    friction.bath._get_non_markovian_jacobian()
    friction.bath.s[:] = np.array([[[0.2, -0.1]]])

    friction.step(0.1)

    assert np.all(np.isfinite(friction.bath.s))
    assert np.all(np.isfinite(friction.nm.pnm))
    assert np.linalg.norm(friction.bath.s) > 0.0
    assert np.linalg.norm(friction.nm.pnm - p0) > 0.0


@pytest.mark.parametrize("eta0", [1.5, -1.0, float("nan"), float("inf")])
def test_doublewell_rejects_nonzero_eta0(eta0):
    with pytest.raises(ValueError, match="eta0 no longer scales"):
        DoubleWell_with_friction_driver(eta0=eta0)


def test_doublewell_zero_eta0_matches_omitted_argument():
    parameters = dict(
        w_b=500.0,
        v0=2085.0,
        m=1837.36223469,
        eps1=-1.0,
        eps2=0.2,
        delta=0.0,
        deltaQ=1.0,
    )
    pos = np.array([[0.3, 0.0, 0.0]])
    results = [
        DoubleWell_with_friction_driver(**parameters, **extra).compute_structure(
            np.zeros((3, 3)), pos
        )
        for extra in ({}, {"eta0": 0.0})
    ]
    for i in range(3):
        np.testing.assert_array_equal(results[0][i], results[1][i])
    assert json.loads(results[0][3]) == json.loads(results[1][3])

"""Regression tests for active-DOF variable non-Markovian friction algebra."""

from __future__ import annotations

import json

import numpy as np
from ipi.utils.depend import depend_value

from ipi.engine.friction import Friction
from ipi.utils import nmtransform


class _IdentityTransform:
    @staticmethod
    def nm2b(array):
        return np.array(array, copy=True)


class _Beads:
    def __init__(self, nbeads: int, natoms: int):
        ndof = 3 * natoms
        self.q = np.zeros((nbeads, ndof))
        self.p = np.linspace(-0.3, 0.4, nbeads * ndof).reshape(nbeads, ndof)
        self.m3 = np.ones((nbeads, ndof))
        self.nbeads = nbeads
        self.natoms = natoms


class _NM:
    @property
    def pnm_(self):
        return self.pnm

    def __init__(self, beads):
        self.pnm = np.array(beads.p, copy=True)
        self.qnm = np.array(beads.q, copy=True)
        self.dynm3 = np.array(beads.m3, copy=True)
        self.omegak = np.zeros(beads.nbeads)
        self.transform = _IdentityTransform()


class _Forces:
    def __init__(self, extras):
        self.extras = extras


class _Ensemble:
    def __init__(self, forces):
        self.temp = 300.0
        self.forces = forces


class _PRNG:
    def gvec(self, shape):
        return np.zeros(shape)


def _extras(reduced_jacobian):
    gamma = np.einsum("bci,bcj->bij", reduced_jacobian, reduced_jacobian)
    meta = {
        "schema": "ipi_friction_v1",
        "active_atoms": [1],
        "channel_labels": ["a", "b"],
        "jacobian_units": "ipi_internal_mass_scaled",
        "gamma_units": "ipi_internal_mass_scaled",
    }
    return {
        "friction_coupling_jacobian": reduced_jacobian.tolist(),
        "friction_gamma": gamma.tolist(),
        "friction_meta": [json.dumps(meta) for _ in range(reduced_jacobian.shape[0])],
    }


def _friction():
    beads = _Beads(nbeads=3, natoms=2)
    # Driver sends only atom 1's three Cartesian columns.
    reduced = np.array(
        [
            [[1.0, 0.1, -0.2], [0.3, -0.4, 0.2]],
            [[0.2, 0.5, 0.1], [-0.1, 0.4, 0.6]],
            [[0.7, -0.3, 0.2], [0.4, 0.2, -0.5]],
        ]
    )
    forces = _Forces(_extras(reduced))
    friction = Friction(
        variable_friction=True,
        bath_mode="non-markovian",
        debug_mf_mode="off",
        Ap=np.array([[0.0, 0.5], [-0.5, 1.5]]),
    )
    friction.beads = beads
    friction.nm = _NM(beads)
    friction.ensemble = _Ensemble(forces)
    friction.forces = forces
    friction.prng = _PRNG()
    friction._ensure_bath_bound()
    return friction


def test_active_normal_mode_transform_matches_full_zero_padded_reference():
    friction = _friction()
    full_jacobian = np.asarray(friction.coupling_jacobian)
    cmat = np.asarray(nmtransform.mk_nm_matrix(friction.beads.nbeads), dtype=float)
    reference = np.einsum("rb,nb,bci->rnci", cmat, cmat, full_jacobian, optimize=True)

    actual = friction._get_nm_coupling_jacobian()
    active = friction._friction_dof_idx
    np.testing.assert_allclose(actual, reference[..., active], rtol=1e-13, atol=1e-13)


def test_mean_field_force_is_scattered_back_to_the_correct_full_dofs():
    friction = _friction()
    full_jacobian = np.asarray(friction.coupling_jacobian)
    cmat = np.asarray(nmtransform.mk_nm_matrix(friction.beads.nbeads), dtype=float)
    reference_jacobian = np.einsum(
        "rb,nb,bci->rnci", cmat, cmat, full_jacobian, optimize=True
    )
    coupling_nm = np.array([[0.2, -0.1], [0.3, 0.5], [-0.4, 0.6]])
    friction.alpha = np.array([0.7, 0.2, 0.9])
    friction.friction_coupling_nm_ = depend_value(
        name="friction_coupling_nm", value=coupling_nm
    )

    reference_force = -np.einsum(
        "r,rc,rnci->ni", friction.alpha, coupling_nm, reference_jacobian, optimize=True
    )
    actual_force = friction.get_force_mf_nm()
    np.testing.assert_allclose(actual_force, reference_force, rtol=1e-13, atol=1e-13)


def test_non_markovian_step_builds_one_normal_mode_jacobian_per_operator_call(
    monkeypatch,
):
    friction = _friction()
    bath = friction.bath
    calls = 0
    original = bath._get_nm_transform_matrix

    def counted_transform():
        nonlocal calls
        calls += 1
        return original()

    monkeypatch.setattr(bath, "_get_nm_transform_matrix", counted_transform)
    friction.step(0.1)

    # B-P-F, B-s, B-P-F all share one coordinate-dependent transform.
    assert calls == 1

    friction.step(0.1)
    # The cache is scoped to one operator application, not retained across
    # the next outer splitting operation where coordinates could be new.
    assert calls == 2

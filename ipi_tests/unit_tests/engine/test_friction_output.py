"""Focused tests for compact and scalar friction diagnostics."""

from types import SimpleNamespace

import numpy as np
import pytest

from ipi.engine.properties import Properties, Trajectories


def friction_system(friction):
    return SimpleNamespace(
        motion=SimpleNamespace(integrator=SimpleNamespace(friction=friction))
    )


def test_active_gamma_matrix_preserves_declared_dof_order():
    gamma = np.arange(2 * 6 * 6, dtype=float).reshape(2, 6, 6)
    friction = SimpleNamespace(
        gamma=gamma,
        variable_friction=True,
        _friction_dof_idx=np.array([5, 3, 4]),
    )
    trajectories = Trajectories()
    trajectories.system = friction_system(friction)

    compact = trajectories.get_friction_gamma_active_matrix()

    expected = gamma[:, [5, 3, 4], :][:, :, [5, 3, 4]]
    np.testing.assert_array_equal(compact, expected)
    assert compact.shape == (2, 3, 3)


def test_active_gamma_matrix_requires_variable_friction_metadata():
    friction = SimpleNamespace(
        gamma=np.eye(3)[None, ...],
        variable_friction=True,
        _friction_dof_idx=None,
    )
    trajectories = Trajectories()
    trajectories.system = friction_system(friction)

    with pytest.raises(ValueError, match="active Cartesian DOFs"):
        trajectories.get_friction_gamma_active_matrix()


def test_scalar_friction_diagnostics():
    friction = SimpleNamespace(
        gamma=np.eye(3)[None, ...],
        _gamma_reconstruction_error=2.5e-13,
        bath=SimpleNamespace(s=np.array([[3.0, 4.0]])),
    )
    properties = Properties()
    properties.motion = friction_system(friction).motion

    assert properties.get_friction_gamma_reconstruction_error() == 2.5e-13
    assert properties.get_friction_auxiliary_norm() == 5.0

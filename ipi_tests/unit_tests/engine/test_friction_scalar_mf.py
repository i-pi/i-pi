"""The scalar bath couples independently to each Cartesian coordinate."""

from types import SimpleNamespace

import numpy as np
import pytest

from ipi.engine.friction import Friction


@pytest.mark.parametrize("ndof", [3, 6])
@pytest.mark.parametrize("nbeads", [1, 3])
def test_scalar_mf_has_independent_cartesian_channels(ndof, nbeads):
    q = np.zeros((nbeads, ndof))
    q[0] = 1.0  # A centroid displacement must have no MF cost.
    if nbeads > 1:
        q[1, 0] = 1.0  # Must not create forces on other coordinates.
        q[2, 0], q[2, -1] = 1.0, -1.0  # Must not cancel the MF energy.
    alpha = np.arange(nbeads, dtype=float)
    friction = Friction(
        variable_friction=False, sigma_static=2.0,
        debug_mf_mode="on", bath_mode="none",
    )
    friction.nm = SimpleNamespace(
        qnm=q, transform=SimpleNamespace(nm2b=lambda a: a.copy()),
    )
    friction.beads = SimpleNamespace(p=np.zeros_like(q), nbeads=nbeads)
    friction.alpha = alpha

    def energy():
        # Refresh explicitly: this focused fixture has no bound dependency graph.
        friction.friction_coupling_nm_.taint(taintme=True)
        return friction.get_energy_mf()

    expected_force = -4.0 * alpha[:, None] * q
    assert energy() == pytest.approx(0.5 * np.sum(-expected_force * q))
    force = friction.get_force_mf_nm()
    assert force.shape == q.shape
    np.testing.assert_allclose(force, expected_force)

    # Check the force against the energy gradient for every Cartesian DOF.
    delta = 1e-6
    for index in np.ndindex(q.shape):
        original = q[index]
        q[index] = original + delta
        plus = energy()
        q[index] = original - delta
        minus = energy()
        q[index] = original
        assert -(plus - minus) / (2 * delta) == pytest.approx(
            expected_force[index], abs=1e-8,
        )
    energy()

    # Exercise the actual momentum kick, including its force shape.
    friction.step(0.25)
    np.testing.assert_allclose(friction.beads.p, 0.25 * expected_force)

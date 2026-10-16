"""Constant friction gives a quadratic action with an independently known Hessian."""

from types import SimpleNamespace

import numpy as np
import pytest

from ipi.engine.motion.instanton import FrictionMapper
from ipi.utils.units import unit_to_internal
from ipi_tests.unit_tests.engine.test_simulation_invariants import harmonic_input


@pytest.mark.parametrize("count", [4, 7])
def test_constant_friction_action_gradient_and_hessian(simulation, count):
    with simulation(harmonic_input(nbeads=count)) as sim:
        s = sim.syslist[0]
        eta = np.diag(np.arange(1.0, 7.0))
        mapper = FrictionMapper(False, eta)
        mapper.bind(
            SimpleNamespace(
                beads=s.beads,
                cell=s.cell,
                forces=s.forces,
                nm=s.nm,
                rp_factor=1.0,
                fix=SimpleNamespace(nactive=2),
                coef=np.ones((count + 1, 1)),
                options={"max_ms": 0, "max_e": 0},
            )
        )
        # Constant spectral density over every normal-mode frequency.
        conversion = unit_to_internal("frequency", "inversecm", 1.0)
        data = np.column_stack(
            [np.linspace(0, 2 * max(mapper.omegak) / conversion, 8), np.ones(8)]
        )
        mapper.set_fric_spec_dens(data, 0)
        q = np.linspace(-0.2, 0.4, count)[:, None] * np.arange(1.0, 7.0)[None, :]
        mapper.dbeads.q[:] = q
        mapper.initialize(q, mapper.dforces)
        energy, gradient = mapper.compute_friction_terms()
        n = np.arange(count)
        transform = np.sqrt(2 / count) * np.cos(
            np.pi * n[:, None] * (n[None, :] + 0.5) / count
        )
        transform[0] /= np.sqrt(2)
        omega = 2 * s.nm.omegan * np.sin(np.pi * n / (2 * count))
        bead_matrix = transform.T @ np.diag(omega) @ transform
        expected_gradient = bead_matrix @ q @ eta
        assert energy == pytest.approx(np.sum(q * expected_gradient) / 2, rel=1e-7)
        np.testing.assert_allclose(gradient, expected_gradient, rtol=1e-7, atol=1e-12)
        hessian = mapper.get_fric_rp_hessian(
            np.zeros((count, 6, 6, 6)), mapper.eta, False
        )
        np.testing.assert_allclose(
            hessian, np.kron(bead_matrix, eta), rtol=1e-7, atol=1e-12
        )
        assert np.linalg.eigvalsh(hessian).min() > -1e-12
        # A rigid shift of the whole path must not change the friction action.
        mapper.dbeads.q[:] = q + 0.7
        shifted_energy, shifted_gradient = mapper.compute_friction_terms()
        assert shifted_energy == pytest.approx(energy, abs=1e-12)
        np.testing.assert_allclose(shifted_gradient, gradient, atol=1e-12)

"""Analytic checks for fixed-replica transformations and spring propagation."""

import numpy as np
import pytest

from ipi.engine.beads import Beads
from ipi.engine.ensembles import Ensemble
from ipi.engine.motion import Motion
from ipi.engine.normalmodes import NormalModes
from ipi.utils.depend import dstrip
from ipi.utils.nmtransform import mk_fixed_bead_matrix
from ipi.utils.units import Constants


def make_nm(q, p, masses, omegan, dt, fixed):
    """Binds the real spring propagator without a physical forcefield."""
    beads = Beads(nbeads=len(q), natoms=q.shape[1] // 3)
    beads.q[:] = q
    beads.p[:] = p
    beads.m[:] = masses
    motion = Motion(fixbeads=fixed)
    motion.dt = dt
    nm = NormalModes(transform_method="matrix")
    nm.bind(
        Ensemble(temp=omegan * Constants.hbar / (len(q) * Constants.kb)),
        motion,
        beads=beads,
    )
    return nm


def spring_energy(nm):
    return 0.5 * np.sum(dstrip(nm.beads.p) ** 2 / dstrip(nm.beads.m3)) + nm.vspring


def test_single_fixed_bead_equilibrium_is_constant_path():
    free, fixed, _, _, equilibrium_map = mk_fixed_bead_matrix(8, [3])
    fixed_q = np.array([[1.25, -0.4]])
    np.testing.assert_array_equal(fixed, [3])
    np.testing.assert_array_equal(free, [0, 1, 2, 4, 5, 6, 7])
    np.testing.assert_allclose(equilibrium_map @ fixed_q, np.repeat(fixed_q, 7, axis=0))


def test_two_bead_propagator_matches_exact_anchored_oscillator():
    nm = make_nm(
        np.array([[-0.2, 0, 0], [0.7, 0, 0]]),
        np.array([[4.2, 0, 0], [-0.9, 0, 0]]),
        [2.3],
        1.7,
        0.31,
        [0],
    )
    frequency = np.sqrt(2.0) * 1.7
    expected_q = -0.2 + 0.9 * np.cos(frequency * 0.31)
    expected_q -= 0.9 * np.sin(frequency * 0.31) / (2.3 * frequency)
    expected_p = -0.9 * np.cos(frequency * 0.31)
    expected_p -= 2.3 * frequency * 0.9 * np.sin(frequency * 0.31)
    nm.free_qstep_fixed_beads()
    assert nm.beads.q[0, 0] == -0.2
    np.testing.assert_array_equal(nm.beads.p[0], 0.0)
    np.testing.assert_allclose(nm.beads.q[1], [expected_q, 0, 0], rtol=1e-13)
    np.testing.assert_allclose(nm.beads.p[1], [expected_p, 0, 0], rtol=1e-13)


def test_fixed_bead_propagator_conserves_energy_and_is_reversible():
    rng = np.random.default_rng(9234)
    q = rng.normal(size=(7, 6))
    p = rng.normal(size=(7, 6))
    p[[1, 5]] = 0.0
    nm = make_nm(q, p, [1.0, 3.0], 2.4, 0.17, [1, 5])
    energy = spring_energy(nm)
    nm.free_qstep_fixed_beads()
    np.testing.assert_allclose(spring_energy(nm), energy, rtol=2e-14)
    nm.motion.dt *= -1
    nm.free_qstep_fixed_beads()
    np.testing.assert_allclose(nm.beads.q, q, rtol=2e-14, atol=2e-14)
    np.testing.assert_allclose(nm.beads.p, p, rtol=2e-14, atol=2e-14)


def test_fixed_bead_indices_are_validated():
    with pytest.raises(ValueError, match="At least one"):
        mk_fixed_bead_matrix(4, [])
    with pytest.raises(ValueError, match="out of bounds"):
        mk_fixed_bead_matrix(4, [-1])
    with pytest.raises(ValueError, match="out of bounds"):
        mk_fixed_bead_matrix(4, [4])


def test_all_fixed_beads_preserve_positions_and_remove_momenta():
    q = np.arange(12, dtype=float).reshape(4, 3)
    p = np.ones_like(q)
    nm = make_nm(q, p, [1.0], 2.0, 0.1, [3, 2, 1, 0])
    nm.free_qstep_fixed_beads()
    np.testing.assert_array_equal(nm.beads.q, q)
    np.testing.assert_array_equal(nm.beads.p, 0.0)
    np.testing.assert_array_equal(p, 1.0)


def test_two_anchors_give_piecewise_linear_minimum_energy_path():
    free, fixed, _, _, equilibrium_map = mk_fixed_bead_matrix(6, [0, 3])
    np.testing.assert_array_equal(fixed, [0, 3])
    np.testing.assert_array_equal(free, [1, 2, 4, 5])
    np.testing.assert_allclose(
        equilibrium_map @ np.array([[0.0], [3.0]]),
        [[1.0], [2.0], [2.0], [1.0]],
        atol=1e-14,
    )


def test_zero_spring_frequency_reduces_to_free_particle_motion():
    q = np.arange(12, dtype=float).reshape(4, 3)
    nm = make_nm(q, np.full_like(q, 2.0), [3.0], 0.0, 0.3, [1])
    nm.free_qstep_fixed_beads()
    np.testing.assert_allclose(nm.beads.q[[0, 2, 3]], q[[0, 2, 3]] + 0.2, atol=1e-14)
    np.testing.assert_array_equal(nm.beads.q[1], q[1])
    np.testing.assert_allclose(nm.beads.p[[0, 2, 3]], 2.0, atol=1e-14)
    np.testing.assert_array_equal(nm.beads.p[1], 0.0)


@pytest.mark.parametrize("change", ["dt", "temperature", "mass"])
def test_cached_propagator_updates_for_timestep_temperature_and_mass(change):
    """Dependency updates retain the exact two-bead oscillator trajectory."""
    nm = make_nm(
        np.array([[-0.2, 0, 0], [0.7, 0, 0]]),
        np.array([[0, 0, 0], [-0.9, 0, 0]]),
        [2.3],
        1.7,
        0.31,
        [0],
    )
    previous = dstrip(nm.fixed_prop_pq_ms).copy()
    nm.beads.q[1, 1] = 0.1
    nm.beads.p[1, 1] = 0.2
    assert not nm._fixed_prop_pq_ms.tainted()
    if change == "dt":
        nm.motion.dt *= 1.5
    elif change == "temperature":
        nm.ensemble.temp *= 1.3
    else:
        nm.beads.m *= 1.4
    assert nm._fixed_prop_pq_ms.tainted()
    assert not np.allclose(nm.fixed_prop_pq_ms, previous)
    q = dstrip(nm.beads.q).copy()
    p = dstrip(nm.beads.p).copy()
    frequency = np.sqrt(2) * nm.omegan
    phase = frequency * nm.dt
    expected_q = q[0] + (q[1] - q[0]) * np.cos(phase)
    expected_q += p[1] * np.sin(phase) / (nm.beads.m[0] * frequency)
    expected_p = p[1] * np.cos(phase)
    expected_p -= nm.beads.m[0] * frequency * (q[1] - q[0]) * np.sin(phase)
    nm.free_qstep_fixed_beads()
    np.testing.assert_allclose(nm.beads.q[1], expected_q, rtol=1e-13, atol=1e-14)
    np.testing.assert_allclose(nm.beads.p[1], expected_p, rtol=1e-13, atol=1e-14)
    assert not nm._fixed_prop_pq_ms.tainted()

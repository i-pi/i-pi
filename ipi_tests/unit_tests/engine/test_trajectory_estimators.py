"""Per-atom estimators checked against harmonic covariance and free-particle limits."""

import numpy as np
import pytest

from ipi.engine.properties import Trajectories
from ipi_tests.unit_tests.engine.test_simulation_invariants import harmonic_input


def test_trajectory_kinetic_tensor_and_gyration_from_covariance(simulation):
    with simulation(harmonic_input(nbeads=4)) as sim:
        s = sim.syslist[0]
        q = np.arange(24).reshape(4, 6) / 20
        s.beads.q[:] = q
        trajectory = Trajectories()
        trajectory.bind(s)
        dq = q - q.mean(axis=0)
        covariance = (
            np.einsum("bai,baj->aij", dq.reshape(4, 2, 3), dq.reshape(4, 2, 3)) / 4
        )
        np.testing.assert_allclose(
            trajectory.get_akcv(),
            0.5 * s.ensemble.temp
            + 0.25 * np.diagonal(covariance, axis1=1, axis2=2).ravel(),
        )
        np.testing.assert_allclose(
            trajectory.get_rg(),
            np.sqrt(np.diagonal(covariance, axis1=1, axis2=2)).ravel(),
        )
        np.testing.assert_allclose(
            trajectory.get_akcv_od().reshape(2, 3),
            0.25 * covariance[:, [0, 0, 1], [1, 2, 2]],
        )


@pytest.mark.parametrize("alpha", [0.5, 1.0, 2.0])
@pytest.mark.parametrize("atom,selected", [("", [0, 1]), ("H", [0]), ("1", [1])])
def test_per_atom_isotope_weights_match_free_particle_action(
    simulation, alpha, atom, selected
):
    with simulation(harmonic_input(nbeads=4).replace("{k1: 0.5}", "{k1: 0.0}")) as sim:
        s = sim.syslist[0]
        trajectory = Trajectories()
        trajectory.bind(s)
        q = s.beads.q.copy().reshape(4, 2, 3)
        delta = (
            (alpha - 1)
            * 0.5
            * np.array([2.0, 3.0])
            * (4 * s.ensemble.temp) ** 2
            * np.sum((q - np.roll(q, 1, axis=0)) ** 2, axis=(0, 2))
        )
        delta[[i for i in range(2) if i not in selected]] = 0
        expected = np.column_stack(
            [delta, delta**2, np.exp(-delta / (4 * s.ensemble.temp))]
        )
        np.testing.assert_allclose(
            trajectory.get_isotope_zetatd(str(alpha), atom).reshape(2, 3),
            expected,
            atol=1e-14,
        )


@pytest.mark.parametrize("atom,natoms", [("", 2), ("H", 1), ("1", 1)])
def test_fourth_order_estimators_reduce_to_free_particle_limits(
    simulation, atom, natoms
):
    with simulation(harmonic_input(nbeads=4).replace("{k1: 0.5}", "{k1: 0.0}")) as sim:
        s = sim.syslist[0]
        assert s.properties.get_scpottd() == pytest.approx(0)
        assert s.properties.get_sckinop(atom) == pytest.approx(
            1.5 * natoms * s.ensemble.temp
        )
        assert s.properties.get_sckintd(atom) == pytest.approx(
            1.5 * natoms * s.ensemble.temp
        )
        assert s.properties.get_sckinpr() == pytest.approx(
            3 * 4 * s.ensemble.temp - s.nm.vspring / 4
        )


@pytest.mark.parametrize("atom,bead", [("", 0), ("2", 0), ("0", 4)])
def test_atom_vector_rejects_out_of_range_indices(simulation, atom, bead):
    with simulation(harmonic_input(nbeads=4)) as sim:
        s = sim.syslist[0]
        with pytest.raises(IndexError):
            s.properties.get_atom_vec(s.beads.q, atom, str(bead))

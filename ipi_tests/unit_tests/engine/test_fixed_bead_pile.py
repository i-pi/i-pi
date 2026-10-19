"""Constrained PILE-L: exact thermostat law and conditional canonical density."""

import io

import numpy as np
import pytest

from ipi.engine.thermostats import ThermoPILE_G, ThermoPILE_L
from ipi.scripting import InteractiveSimulation
from ipi.utils.depend import dstrip
from ipi.utils.prng import Random
from ipi.utils.units import Constants

from .test_fixed_bead_springs import make_nm
from . import test_path_constraints

simulations = test_path_constraints.simulations


@pytest.mark.parametrize("fixed", [[0], [0, 3]])
def test_fixed_pile_preserves_conditional_spring_density_and_heat(fixed):
    """Independent canonical samples remain canonical at a finite timestep.

    For springs alone the exact position covariance is
    P k_B T (m omega_P^2 L_ff)^-1, and the momentum covariance is
    m P k_B T I. Vary masses across the independent atoms.
    """
    rng = np.random.default_rng(3184)
    nbeads, natoms = 6, 8000
    masses = np.linspace(1.0, 4.0, natoms)
    nm = make_nm(
        np.zeros((nbeads, 3 * natoms)),
        np.zeros((nbeads, 3 * natoms)),
        masses,
        1.7,
        0.31,
        fixed,
    )
    # Construct the cyclic Laplacian independently of the mode transform.
    laplacian = 2 * np.eye(nbeads)
    for j in range(nbeads):
        laplacian[j, (j + 1) % nbeads] -= 1
        laplacian[j, (j - 1) % nbeads] -= 1
    free = np.setdiff1d(np.arange(nbeads), fixed)
    covariance = np.linalg.inv(laplacian[np.ix_(free, free)])
    thermal = nbeads * Constants.kb * nm.ensemble.temp
    mass = np.repeat(masses, 3)
    nm.beads.q[fixed] = rng.normal(size=(len(fixed), 3 * natoms))
    anchors = dstrip(nm.beads.q)[fixed].copy()
    equilibrium = -covariance @ laplacian[np.ix_(free, fixed)] @ anchors
    nm.beads.q[free] = (
        equilibrium
        + (np.linalg.cholesky(covariance) @ rng.normal(size=(len(free), 3 * natoms)))
        * np.sqrt(thermal / mass)
        / nm.omegan
    )
    nm.beads.p[free] = rng.normal(size=(len(free), 3 * natoms)) * np.sqrt(
        thermal * mass
    )
    thermostat = ThermoPILE_L(
        temp=nbeads * nm.ensemble.temp, dt=0.31, scale=0.7, ethermo=2.4
    )
    thermostat.bind(nm=nm, prng=Random(seed=1524))
    for _ in range(12):
        energy = 0.5 * np.sum(dstrip(nm.beads.p) ** 2 / nm.beads.m3)
        heat = thermostat.ethermo
        thermostat.step()
        new_energy = 0.5 * np.sum(dstrip(nm.beads.p) ** 2 / nm.beads.m3)
        assert new_energy + thermostat.ethermo == pytest.approx(
            energy + heat, rel=2e-14
        )
        nm.free_qstep_fixed_beads()
    q = (dstrip(nm.beads.q)[free] - equilibrium) * np.sqrt(mass / thermal) * nm.omegan
    p = dstrip(nm.beads.p)[free] / np.sqrt(mass * thermal)
    # 24000 independent Cartesian samples; tolerances exceed sampling errors.
    np.testing.assert_allclose(q.mean(axis=1), 0, atol=0.025)
    np.testing.assert_allclose(p.mean(axis=1), 0, atol=0.025)
    np.testing.assert_allclose(q @ q.T / q.shape[1], covariance, atol=0.035)
    np.testing.assert_allclose(p @ p.T / p.shape[1], np.eye(len(free)), atol=0.035)
    np.testing.assert_allclose(q @ p.T / p.shape[1], 0, atol=0.035)
    np.testing.assert_array_equal(nm.beads.q[fixed], anchors)
    np.testing.assert_array_equal(nm.beads.p[fixed], 0)


def test_fixed_pile_exact_ou_update_and_dependency_changes():
    """One fixed endpoint has the known sine-mode spectrum, not ring modes."""
    nm = make_nm(
        np.zeros((4, 6)), np.arange(24.0).reshape(4, 6), [1.0, 3.0], 1.7, 0.2, [0]
    )
    thermostat = ThermoPILE_L(temp=10000, dt=0.13, scale=0.6, tau=57)
    thermostat.bind(nm=nm, prng=Random(seed=385), fixdof=0, fixcom=False)
    reference = Random(seed=385)
    # A fixed replica cuts the ring into a chain with two fixed endpoints.
    frequencies = 2 * nm.omegan * np.sin(np.arange(1, 4) * np.pi / 8)
    np.testing.assert_allclose(1 / thermostat.tauk, 2 * 0.6 * frequencies)
    for _ in range(2):
        initial = dstrip(nm.beads.p).copy()
        expected = nm.fixed_transform.T @ initial[nm.freebeads]
        damping = np.exp(-2 * thermostat.pilescale * frequencies * thermostat.dt)
        for k in range(3):
            expected[k] = damping[k] * expected[k] + reference.gvec(6) * np.sqrt(
                nm.beads.m3[0] * Constants.kb * thermostat.temp * (1 - damping[k] ** 2)
            )
        thermostat.step()
        np.testing.assert_allclose(
            nm.beads.p[nm.freebeads], nm.fixed_transform @ expected, atol=2e-14
        )
        np.testing.assert_array_equal(nm.beads.p[0], initial[0])
        # Changing ensemble temperature changes spring frequencies; thermostat
        # temperature, dt, masses and pile_lambda must update the OU law too.
        nm.ensemble.temp *= 1.3
        frequencies *= 1.3
        thermostat.temp *= 1.2
        thermostat.dt *= 0.8
        thermostat.pilescale *= 0.7
        nm.beads.m *= 1.4


@pytest.mark.parametrize("thermostat_class", [ThermoPILE_L, ThermoPILE_G])
def test_ordinary_pile_internal_modes_keep_their_existing_ou_law(thermostat_class):
    nm = make_nm(np.zeros((4, 3)), np.arange(12.0).reshape(4, 3), [2.0], 1.7, 0.2, [])
    thermostat = thermostat_class(temp=10000, dt=0.13, tau=0.7, scale=0.6)
    thermostat.bind(nm=nm, prng=Random(seed=385), fixdof=0, fixcom=False)
    initial = dstrip(nm.pnm).copy()
    # Each non-centroid mode retains its ordinary-ring frequency.
    frequencies = np.array([np.sqrt(2), 2, np.sqrt(2)]) * nm.omegan
    np.testing.assert_allclose(1 / thermostat.tauk, 2 * 0.6 * frequencies)
    # PILE-G replaces only the centroid; use actual child random states to
    # check the internal OU law independently of that centroid's random draws.
    for k in range(1, 4):
        child = thermostat._thermos[k]
        random = Random(state=child.prng.state)
        a = np.exp(-2 * 0.6 * frequencies[k - 1] * thermostat.dt)
        expected = a * initial[k] + random.gvec(3) * np.sqrt(
            nm.dynm3[k] * Constants.kb * thermostat.temp * (1 - a**2)
        )
        child.step()
        np.testing.assert_allclose(nm.pnm[k], expected, atol=2e-14)


@pytest.mark.parametrize("nbeads", [1, 4])
def test_all_fixed_pile_is_a_noop_and_keeps_restart_heat(simulations, nbeads):
    system = simulations(
        constraints=f"<fixbeads>{list(range(nbeads))}</fixbeads>",
        nbeads=nbeads,
        mode="nvt",
        thermostat="pile_l",
    ).syslist[0]
    q = dstrip(system.beads.q).copy()
    system.motion.thermostat.ethermo = 7.2
    for _ in range(5):
        system.motion.step()
    np.testing.assert_array_equal(system.beads.q, q)
    np.testing.assert_array_equal(system.beads.p, 0)
    assert system.motion.thermostat.ethermo == 7.2


@pytest.mark.parametrize("splitting", ["obabo", "baoab"])
def test_fixed_pile_with_fixed_atomic_dofs(simulations, splitting):
    system = simulations(
        constraints="<fixbeads>[0,2]</fixbeads><fixatoms_dof>[0,4]</fixatoms_dof>",
        mode="nvt",
        splitting=splitting,
        thermostat="pile_l",
    ).syslist[0]
    system.beads.q[:, [0, 4]] = system.beads.q[0, [0, 4]]
    q = dstrip(system.beads.q).copy()
    for _ in range(20):
        system.motion.step()
    np.testing.assert_allclose(system.beads.q[:, [0, 4]], q[:, [0, 4]], atol=2e-14)
    np.testing.assert_allclose(system.beads.p[:, [0, 4]], 0, atol=2e-14)
    np.testing.assert_array_equal(system.beads.q[[0, 2]], q[[0, 2]])
    assert np.isfinite(system.ensemble.eens)


def test_fixed_pile_rejects_a_centroid_temperature(simulations):
    xml = test_path_constraints.constraint_xml(
        constraints="<fixbeads>[0]</fixbeads>", mode="nvt", thermostat="pile_l"
    ).replace("<tau>0.5</tau>", "<tau>0.5</tau><pile_centroid_t>100</pile_centroid_t>")
    with pytest.raises(ValueError, match="no centroid temperature"):
        InteractiveSimulation(io.StringIO(xml))


def test_fixed_bead_pile_g_is_rejected_when_bound_directly():
    nm = make_nm(np.zeros((4, 3)), np.zeros((4, 3)), [1.0], 1.7, 0.2, [0])
    with pytest.raises(ValueError, match="PILE-L, not PILE-G"):
        ThermoPILE_G().bind(nm=nm)

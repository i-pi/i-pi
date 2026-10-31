"""Thermostat energy accounting and deterministic damping limits."""

import numpy as np
import pytest

from ipi.engine.atoms import Atoms
from ipi.engine.thermostats import (
    MultiThermo,
    ThermoCL,
    ThermoFFL,
    ThermoGLE,
    ThermoLangevin,
    ThermoSVR,
    _svr_step,
)
from ipi.utils.prng import Random


@pytest.mark.parametrize(
    "kind", ["langevin", "svr", "gle", "soft", "hard", "rescale", "none"]
)
def test_thermostat_accounts_for_exchanged_kinetic_energy(kind):
    """For an isolated thermostat step, kinetic energy plus bath energy is conserved."""
    atoms = Atoms(2)
    atoms.m_[:] = [2, 5]
    atoms.p_[:] = [0.1, -0.2, 0.3, -0.4, 0.5, 0.2]
    options = dict(temp=300, dt=0.1, ethermo=1.25)
    if kind == "langevin":
        thermostat = ThermoLangevin(tau=0.5, **options)
    elif kind == "svr":
        thermostat = ThermoSVR(tau=0.5, **options)
    elif kind == "gle":
        thermostat = ThermoGLE(A=np.array([[2.0, -0.2], [0.2, 1.0]]), **options)
    else:
        thermostat = ThermoFFL(tau=0.5, flip=kind, **options)
    thermostat.bind(atoms=atoms, prng=Random(seed=19), fixdof=0)
    original = atoms.kin + thermostat.ethermo
    for _ in range(5):
        thermostat.step()
        assert atoms.kin + thermostat.ethermo == pytest.approx(
            original, rel=1e-12, abs=1e-12
        )
        assert np.isfinite(atoms.p).all()


def test_zero_temperature_langevin_is_exact_damping():
    """Without thermal noise, each momentum decays exponentially."""
    atoms = Atoms(1)
    atoms.m_[:] = 3
    atoms.p_[:] = [1, -2, 3]
    original = atoms.p.copy()
    thermostat = ThermoLangevin(temp=0, dt=0.2, tau=0.5)
    thermostat.bind(atoms=atoms, prng=Random(seed=29))
    thermostat.step()
    np.testing.assert_allclose(atoms.p, original * np.exp(-0.4), atol=1e-14)
    # Changing a coupling parameter must update cached damping coefficients.
    thermostat.tau = 0.25
    thermostat.step()
    np.testing.assert_allclose(atoms.p, original * np.exp(-1.2), atol=1e-14)


def test_finite_temperature_langevin_matches_published_update(monkeypatch):
    """Check damping and thermal kicks against the published Langevin update."""
    # Ceriotti et al., JCP 133, 124104 (2010), Eqs. (28), (30), (31):
    # https://arxiv.org/pdf/1009.1045
    # For one bead, h is the thermostat step duration, gamma=1/tau, and
    # p' = exp(-h/tau)*p + sqrt(m*kB*T*(1-exp(-2*h/tau)))*xi.
    # All inputs below are in atomic units, where kB=1.
    masses = np.array([2.0, 5.0])
    initial = np.array([1.0, -2.0, 0.0, -0.4, 0.5, 0.2])
    draws = np.array([-1.0, 0.5, 2.0, 1.5, -0.25, 0.0])
    temperature, timestep, tau = 0.2, 0.2, 0.5
    atoms = Atoms(2)
    atoms.m_[:] = masses
    atoms.p_[:] = initial
    prng = Random(seed=37)

    def fixed_gaussian_values(shape):
        assert shape == draws.shape
        return draws.copy()

    monkeypatch.setattr(prng, "gvec", fixed_gaussian_values)
    thermostat = ThermoLangevin(temp=temperature, dt=timestep, tau=tau)
    thermostat.bind(atoms=atoms, prng=prng)
    damping = np.exp(-timestep / tau)
    expected = (
        damping * initial
        + np.sqrt(np.repeat(masses, 3) * temperature * (1 - damping**2)) * draws
    )
    thermostat.step()
    np.testing.assert_allclose(atoms.p, expected, rtol=1e-12, atol=1e-14)


@pytest.mark.parametrize("natoms", [1, 2])
def test_svr_rescaling_matches_published_update_with_fixed_draws(natoms):
    """Check the common momentum scale against the paper's Gaussian construction."""
    # Bussi, Donadio and Parrinello, JCP 126, 014101 (2007), Appendix A:
    # https://arxiv.org/pdf/0803.4060
    # Eqs. (A2)-(A7) construct the new kinetic energy from squared Gaussian
    # variables. This avoids repeating the implementation's expanded alpha**2.
    masses = np.array([2.0, 5.0])[:natoms]
    initial = np.array([1.0, -2.0, 3.0, -0.4, 0.5, 0.2])[: 3 * natoms]
    draws = np.array([0.5, -1.0, 0.25, 0.75, -0.5, 1.25])[: 3 * natoms]
    temperature, timestep, tau = 0.2, 0.2, 0.5
    atoms = Atoms(natoms)
    atoms.m_[:] = masses
    atoms.p_[:] = initial

    class FixedRandom:
        def __init__(self):
            self.values = iter(draws)

        def gvec(self, shape):
            return np.array([next(self.values) for _ in range(shape)])

        def gammavec(self, k, shape):
            # Gamma(k, scale=1) = half a sum of 2*k squared standard
            # Gaussian variables, for the integer values of 2*k used here.
            return np.array(
                [
                    0.5 * sum(next(self.values) ** 2 for _ in range(int(2 * k)))
                    for _ in range(shape)
                ]
            )

    # Atomic units: kB=1. The positive first draw gives a positive scale.
    kinetic = np.sum(initial**2 / (2 * np.repeat(masses, 3)))
    damping = np.exp(-timestep / tau)
    auxiliary = np.sqrt(1 - damping) * draws
    auxiliary[0] += np.sqrt(damping * 2 * kinetic / temperature)
    expected_kinetic = 0.5 * temperature * np.dot(auxiliary, auxiliary)
    expected_scale = np.sqrt(expected_kinetic / kinetic)
    thermostat = ThermoSVR(temp=temperature, dt=timestep, tau=tau)
    thermostat.bind(atoms=atoms, prng=FixedRandom())
    thermostat.step()
    assert atoms.p[0] / initial[0] == pytest.approx(expected_scale, rel=1e-12)
    np.testing.assert_allclose(
        atoms.p, expected_scale * initial, rtol=1e-12, atol=1e-14
    )


@pytest.mark.parametrize("r1", [-1.5, -1.0, -0.5, 0.0, 0.5])
def test_svr_with_one_degree_of_freedom_is_langevin(r1):
    """With one degree of freedom, the SVR step is the Langevin propagator,
    including the sign of the new momentum."""
    # Bussi, Donadio and Parrinello, JCP 126, 014101 (2007), Eq. (A2), with
    # x = p / sqrt(m kT): p' = sqrt(et) p + sqrt((1 - et) m kT) r1.
    # With dt = tau, the momentum is reversed for r1 < -0.76.
    p, m, temperature, et = 1.0, 2.0, 0.5, np.exp(-1.0)
    alpha, _ = _svr_step(
        np.array([p]), np.array([m]), 0.5 * temperature, et, np.array([r1]), np.zeros(1)
    )
    expected = np.sqrt(et) * p + np.sqrt((1 - et) * m * temperature) * r1
    assert alpha * p == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize("intau,idtau", [(0.0, 0.0), (0.7, 0.0), (0.0, 0.7)])
@pytest.mark.parametrize("adaptive", [False, True])
def test_compensating_langevin_accounts_for_heat(intau, idtau, adaptive):

    atoms = Atoms(2)
    atoms.m_[:] = [2, 5]
    atoms.p_[:] = [0.1, -0.2, 0.3, -0.4, 0.5, 0.2]
    bath = ThermoCL(
        temp=0.2,
        dt=0.01,
        tau=0.5,
        intau=intau,
        idtau=idtau,
        apat=1.0 if adaptive else 0.0,
    )
    bath.bind(atoms=atoms, prng=Random(seed=71))
    energy = atoms.kin + bath.ethermo
    for _ in range(4):
        bath.step()
        assert atoms.kin + bath.ethermo == pytest.approx(energy, abs=1e-12)
        assert np.isfinite(atoms.p).all()
    if adaptive and intau:
        assert bath.intau > 0 and bath.intau != intau
    if adaptive and idtau:
        assert bath.idtau > 0 and bath.idtau != idtau


def test_multiple_baths_compose_damping_and_sum_heat():

    atoms = Atoms(1)
    atoms.m_[:] = 2
    atoms.p_[:] = [0.2, -0.3, 0.4]
    initial = atoms.p.copy()
    baths = [ThermoLangevin(tau=0.5), ThermoLangevin(tau=1.0)]
    combined = MultiThermo(temp=0.0, dt=0.2, thermolist=baths)
    combined.bind(atoms=atoms, prng=Random(seed=79))
    energy = atoms.kin
    combined.step()
    np.testing.assert_allclose(atoms.p, initial * np.exp(-0.2 * (2 + 1)), atol=1e-14)
    assert atoms.kin + combined.ethermo == pytest.approx(energy)
    combined.dt = 0.1
    combined.step()
    np.testing.assert_allclose(atoms.p, initial * np.exp(-0.3 * (2 + 1)), atol=1e-14)

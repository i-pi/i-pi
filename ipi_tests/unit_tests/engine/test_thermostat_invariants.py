"""Thermostat energy accounting and deterministic damping limits."""

import numpy as np
import pytest

from ipi.engine.atoms import Atoms
from ipi.engine.thermostats import ThermoFFL, ThermoGLE, ThermoLangevin, ThermoSVR
from ipi.utils.prng import Random


@pytest.mark.parametrize(
    "kind", ["langevin", "svr", "gle", "soft", "hard", "rescale", "none"]
)
def test_thermostat_accounts_for_exchanged_kinetic_energy(kind):
    """For an isolated thermostat step, kinetic energy plus bath energy is conserved."""
    atoms = Atoms(2)
    atoms.m[:] = [2, 5]
    atoms.p[:] = [0.1, -0.2, 0.3, -0.4, 0.5, 0.2]
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
    atoms.m[:] = 3
    atoms.p[:] = [1, -2, 3]
    original = atoms.p.copy()
    thermostat = ThermoLangevin(temp=0, dt=0.2, tau=0.5)
    thermostat.bind(atoms=atoms, prng=Random(seed=29))
    thermostat.step()
    np.testing.assert_allclose(atoms.p, original * np.exp(-0.4), atol=1e-14)
    # Changing a coupling parameter must update cached damping coefficients.
    thermostat.tau = 0.25
    thermostat.step()
    np.testing.assert_allclose(atoms.p, original * np.exp(-1.2), atol=1e-14)


def test_svr_leaves_zero_momenta_unchanged():
    """Velocity rescaling has no direction to rescale at zero kinetic energy."""
    atoms = Atoms(1)
    atoms.m[:] = 2
    thermostat = ThermoSVR(temp=300, ethermo=2.5)
    thermostat.bind(atoms=atoms, prng=Random(seed=31))
    thermostat.step()
    np.testing.assert_array_equal(atoms.p, np.zeros(3))
    assert thermostat.ethermo == 2.5


@pytest.mark.parametrize("intau,idtau", [(0.0, 0.0), (0.7, 0.0), (0.0, 0.7)])
@pytest.mark.parametrize("adaptive", [False, True])
def test_compensating_langevin_accounts_for_heat(intau, idtau, adaptive):
    from ipi.engine.thermostats import ThermoCL

    atoms = Atoms(2)
    atoms.m[:] = [2, 5]
    atoms.p[:] = [0.1, -0.2, 0.3, -0.4, 0.5, 0.2]
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
    from ipi.engine.thermostats import MultiThermo

    atoms = Atoms(1)
    atoms.m[:] = 2
    atoms.p[:] = [0.2, -0.3, 0.4]
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

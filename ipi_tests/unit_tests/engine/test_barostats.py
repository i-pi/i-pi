"""Tests for pressure-control algorithms."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

import math
import xml.etree.ElementTree as ET

from types import SimpleNamespace

import numpy as np
import pytest

from ipi.engine.barostats import BaroBZP, BaroMTK, BaroRGB, BaroSCR, cell_propagators
from ipi.engine.cell import Cell
from ipi.engine.motion.dynamics import Dynamics, NPTIntegrator
from ipi.engine.thermostats import ThermoLangevin
from ipi.inputs.barostats import InputBaro
from ipi.utils.array_backend import to_numpy
from ipi.utils.depend import depend_array, depend_value, dproperties
from ipi.utils.units import Constants, unit_to_internal
from ipi_tests.unit_tests.engine.inputs import pressure_input


class MockBeads:
    """Positions and momenta of a single bead."""

    nbeads = 1

    def __init__(self, q, p):
        self.q_ = depend_array(name="q", value=q)
        self.p_ = depend_array(name="p", value=p)


dproperties(MockBeads, ["q", "p"])


class FixedGaussian:
    """Minimal random-number source returning a fixed Gaussian value."""

    def __init__(self, value):
        self.value = value

    @property
    def g(self):
        return self.value


class OrthorhombicCell:
    """Minimal cell object for testing the SCR coordinate map."""

    def __init__(self, lengths):
        self.h_ = depend_array(name="h", value=np.diag(np.asarray(lengths, float)))

    @property
    def V(self):
        return np.linalg.det(self.h)


dproperties(OrthorhombicCell, ["h"])


def scr_simulation_xml():
    """Returns a small self-contained SCR simulation input."""

    return """<simulation verbosity="quiet" floatformat="%24.16e">
  <total_steps>1</total_steps>
  <prng><seed>31415</seed></prng>
  <ffdirect name="harmonic">
    <pes>harmonic</pes>
    <parameters>{k1: 0.1}</parameters>
  </ffdirect>
  <system>
    <forces><force forcefield="harmonic"/></forces>
    <ensemble>
      <temperature units="kelvin">300</temperature>
      <pressure units="bar">1</pressure>
    </ensemble>
    <motion mode="dynamics">
      <fixcom>False</fixcom>
      <dynamics mode="npt">
        <timestep units="femtosecond">0.5</timestep>
        <barostat mode="stochastic-rescaling">
          <tau units="femtosecond">1000</tau>
          <bulk_modulus units="gigapascal">2.2</bulk_modulus>
        </barostat>
        <thermostat mode="langevin">
          <tau units="femtosecond">100</tau>
        </thermostat>
      </dynamics>
    </motion>
    <beads natoms="1" nbeads="1">
      <q shape="(1,3)">[0.1,0,0]</q>
      <p shape="(1,3)">[1,0.2,-0.1]</p>
      <m shape="(1)">[1837.36223469]</m>
      <names shape="(1)">[H]</names>
    </beads>
    <cell shape="(3,3)">
      [18.897261,0,0,0,18.897261,0,0,0,18.897261]
    </cell>
  </system>
</simulation>"""


def simulation_state(sim):
    """Returns the state needed to compare continuous and restarted runs."""

    system = sim.syslist[0]
    return (
        sim.step,
        to_numpy(system.beads.q).copy(),
        to_numpy(system.beads.p).copy(),
        to_numpy(system.cell.h).copy(),
        float(system.motion.barostat.ebaro),
    )


def scr_two_particle_xml(external_pressure, pair_potential=None):
    """Specify a volume of 8, kinetic energy of 3, and kBT of 2."""
    root = ET.fromstring(scr_simulation_xml())
    provider = root.find("ffdirect")
    provider.set("name", "reference")
    provider.find("pes").text = "dummy" if pair_potential is None else "custom"
    provider.find("parameters").text = "{}"
    if pair_potential is not None:
        ET.SubElement(provider, "pes_path").text = str(pair_potential)
    root.find("system/forces/force").set("forcefield", "reference")

    # These inputs use atomic units, including temperature in energy units.
    for xpath, value in [
        ("system/ensemble/temperature", "2.0"),
        ("system/ensemble/pressure", str(external_pressure)),
        ("system/motion/dynamics/timestep", "0.2"),
        ("system/motion/dynamics/barostat/tau", "4.0"),
        ("system/motion/dynamics/barostat/bulk_modulus", "3.0"),
    ]:
        element = root.find(xpath)
        element.attrib.clear()
        element.text = value

    beads = root.find("system/beads")
    beads.set("natoms", "2")
    for name, shape, values in [
        ("q", "(1,6)", "[0.5,0,0,-0.5,0,0]"),
        ("p", "(1,6)", "[1,1,1,-1,-1,-1]"),
        ("m", "(2)", "[1,1]"),
        ("names", "(2)", "[H,H]"),
    ]:
        beads.find(name).set("shape", shape)
        beads.find(name).text = values
    root.find("system/cell").text = "[2,0,0,0,2,0,0,0,2]"
    return ET.tostring(root, encoding="unicode")


@pytest.fixture
def harmonic_pair_potential(tmp_path):
    """Provide the analytic interaction U = (2/2) |q1-q2|^2."""
    potential = tmp_path / "scr_harmonic_pair.py"
    potential.write_text("""import numpy as np
from ipi.pes.dummy import Dummy_driver

__DRIVER_CLASS__ = "HarmonicPair"

class HarmonicPair(Dummy_driver):
    def compute_structure(self, cell, pos):
        separation = pos[0] - pos[1]
        energy = np.dot(separation, separation)
        forces = np.array([-2.0 * separation, 2.0 * separation])
        virial = -2.0 * np.outer(separation, separation)
        return energy, forces, virial, ""
""")
    return potential


@pytest.mark.parametrize(
    "external_pressure,force_coefficient", [(0.1, 0.55), (0.375, 0.0), (0.7, -0.65)]
)
def test_scr_pressure_force_for_noninteracting_particles(
    simulation, external_pressure, force_coefficient
):
    # Bernetti and Bussi (2020), Eq. 5 defines Pint; Eq. 7 gives the drift:
    # https://arxiv.org/pdf/2006.09250
    # K = 6 * (1^2 / 2) = 3, V = 8, so Pint = 2K/(3V) = 1/4.
    # F_lambda = 2 sqrt(V) [Pint - Pext + kBT/(2V)], with kBT/(2V) = 1/8.
    with simulation(scr_two_particle_xml(external_pressure)) as sim:
        system = sim.syslist[0]
        barostat = system.motion.barostat
        assert float(system.nm.kin) == pytest.approx(3.0, rel=1e-12)
        np.testing.assert_allclose(
            to_numpy(system.forces.vir), np.zeros((3, 3)), atol=1e-12
        )
        stress = to_numpy(barostat.stress_mts(0))
        assert np.trace(stress) / 3.0 == pytest.approx(0.25, rel=1e-12)
        assert float(barostat.get_lambda_force()) == pytest.approx(
            force_coefficient * np.sqrt(8.0), rel=1e-12, abs=1e-12
        )


@pytest.mark.parametrize(
    "external_pressure,force_coefficient",
    [(0.1, 23.0 / 60.0), (7.0 / 24.0, 0.0), (0.7, -49.0 / 60.0)],
)
def test_scr_pressure_force_for_harmonic_pair(
    simulation, harmonic_pair_potential, external_pressure, force_coefficient
):
    # At fixed scaled positions, r scales as V^(1/3), so U scales as V^(2/3).
    # Thus Pint = 2K/(3V) - dU/dV = (2K - k*r^2)/(3V).
    # Here K = 3, k = 2, r = 1, V = 8: U = 1, Pint = 1/6.
    # Eq. 7 of Bernetti and Bussi then gives the tabulated forces above.
    with simulation(
        scr_two_particle_xml(external_pressure, harmonic_pair_potential)
    ) as sim:
        system = sim.syslist[0]
        barostat = system.motion.barostat
        assert float(system.nm.kin) == pytest.approx(3.0, rel=1e-12)
        assert float(system.forces.pot) == pytest.approx(1.0, rel=1e-12)
        np.testing.assert_allclose(
            to_numpy(system.forces.vir),
            np.diag([-2.0, 0.0, 0.0]),
            rtol=1e-12,
            atol=1e-12,
        )
        stress = to_numpy(barostat.stress_mts(0))
        assert np.trace(stress) / 3.0 == pytest.approx(1.0 / 6.0, rel=1e-12)
        assert float(barostat.get_lambda_force()) == pytest.approx(
            force_coefficient * np.sqrt(8.0), rel=1e-12, abs=1e-12
        )


def test_scr_input_roundtrip_preserves_restart_state():
    """SCR parameters and accumulated work survive checkpoint storage."""

    bulk_modulus = unit_to_internal("pressure", "gigapascal", 2.2)
    original = BaroSCR(
        tau=100.0,
        bulk_modulus=bulk_modulus,
        ebaro=3.25,
    )

    stored = InputBaro()
    stored.store(original)
    restored = stored.fetch()

    assert type(restored) is BaroSCR
    assert restored.tau == pytest.approx(original.tau)
    assert restored.bulk_modulus == pytest.approx(original.bulk_modulus)
    assert restored.ebaro == pytest.approx(original.ebaro)


def test_scr_restart_matches_continuous_trajectory(restart_states):
    """Checkpointing preserves SCR work, phase, and random-number state."""

    continuous_state, restarted_state = restart_states(
        scr_simulation_xml(), simulation_state, steps=(6, 4)
    )

    assert restarted_state[0] == continuous_state[0]
    for restarted_value, continuous_value in zip(
        restarted_state[1:4], continuous_state[1:4]
    ):
        np.testing.assert_allclose(
            restarted_value, continuous_value, rtol=1e-12, atol=1e-12
        )
    assert restarted_state[4] == pytest.approx(
        continuous_state[4], rel=1e-12, abs=1e-12
    )


def test_scr_requires_explicit_positive_bulk_modulus():
    """SCR cannot silently guess a material bulk modulus."""

    input_baro = InputBaro()
    input_baro._explicit = True
    input_baro.mode.store("stochastic-rescaling")
    input_baro.tau.store(100.0)

    with pytest.raises(ValueError, match="must be specified explicitly"):
        input_baro.fetch()

    input_baro.bulk_modulus.store(-1.0)
    with pytest.raises(ValueError, match="bulk modulus must be positive"):
        input_baro.fetch()


def test_scr_rejects_a_separate_cell_thermostat():
    """SCR has no piston degree of freedom to thermostat."""

    input_baro = InputBaro()
    input_baro._explicit = True
    input_baro.mode.store("stochastic-rescaling")
    input_baro.tau.store(100.0)
    input_baro.bulk_modulus.store(1.0)
    input_baro.thermostat.store(ThermoLangevin(tau=10.0))

    with pytest.raises(ValueError, match="does not use a separate cell thermostat"):
        input_baro.fetch()


def test_scr_rejects_path_integral_dynamics(simulation):
    """The classical SCR equations reject path-integral dynamics."""

    xml = (
        scr_simulation_xml()
        .replace('nbeads="1"', 'nbeads="2"')
        .replace(
            '<q shape="(1,3)">[0.1,0,0]</q>',
            '<q shape="(2,3)">[0.1,0,0,0.1,0,0]</q>',
        )
        .replace(
            '<p shape="(1,3)">[1,0.2,-0.1]</p>',
            '<p shape="(2,3)">[1,0.2,-0.1,1,0.2,-0.1]</p>',
        )
    )

    with pytest.raises(ValueError, match="classical dynamics only"):
        with simulation(xml):
            pass


def test_scr_trotter_map_and_effective_energy():
    """Check the S7/S12 map and energy correction from move probabilities."""

    duration = 0.2
    kbt = 2.0  # Atomic units use kB = 1.
    relaxation_time = 4.0
    bulk_modulus = 3.0
    external_pressure = 0.7
    initial_ebaro = 1.25
    gaussian_draw = 0.4
    barostat = BaroSCR(
        dt=duration,
        temp=kbt,
        tau=relaxation_time,
        bulk_modulus=bulk_modulus,
        pext=external_pressure,
        ebaro=initial_ebaro,
    )
    barostat.prng = FixedGaussian(gaussian_draw)
    barostat.cell = OrthorhombicCell([2.0, 2.0, 2.0])
    positions = np.asarray([[0.2, -0.3, 0.4, 0.5, -0.1, 0.7]])
    momenta = np.asarray([[0.8, -0.4, 0.2, 0.3, 0.5, -0.6]])
    masses = np.asarray([[1.0, 1.0, 1.0, 2.0, 2.0, 2.0]])
    barostat.beads = MockBeads(positions.copy(), momenta.copy())
    barostat.nm = SimpleNamespace(dynm3=masses)
    barostat.qdt_ = depend_value(name="qdt", value=barostat.dt / 2.0)

    old_volume = barostat.cell.V
    old_lambda = np.sqrt(old_volume)
    old_force = 0.6
    new_force = -0.2
    forces = iter((old_force, new_force))
    barostat.get_lambda_force = lambda: next(forces)

    coupling_dt = 2.0 * barostat.qdt
    diffusion = barostat.temp / (4.0 * barostat.tau * barostat.bulk_modulus)
    delta_lambda = (
        diffusion * old_force * coupling_dt / barostat.temp
        + np.sqrt(2.0 * diffusion * coupling_dt) * barostat.prng.g
    )
    scale = ((old_lambda + delta_lambda) / old_lambda) ** (2.0 / 3.0)
    inverse_scale = 1.0 / scale
    drift_scale = 0.5 * (scale + inverse_scale)

    barostat.qcstep()
    assert barostat._scr_state is not None
    barostat.qcstep()

    np.testing.assert_allclose(
        barostat.beads.q[0],
        scale * positions[0] + drift_scale * momenta[0] * barostat.dt / masses[0],
    )
    np.testing.assert_allclose(barostat.beads.p[0], inverse_scale * momenta[0])
    assert barostat.cell.V == pytest.approx(old_volume * scale**3)

    # Bernetti and Bussi (2020), Supplementary Methods, Eqs. S7--S9:
    # https://arxiv.org/src/2006.09250v4/anc/si.pdf
    # Calculate Gaussian transition densities from S7, without expanding their
    # ratio into the force sums and force squares used by finalize(). All
    # reference inputs below are prescribed, including the initial volume 8.
    variance = kbt * duration / (2.0 * relaxation_time * bulk_modulus)
    drift_per_force = duration / (4.0 * relaxation_time * bulk_modulus)
    start_lambda = np.sqrt(8.0)
    forward_mean = start_lambda + drift_per_force * old_force
    end_lambda = forward_mean + np.sqrt(variance) * gaussian_draw
    reverse_mean = end_lambda + drift_per_force * new_force

    def log_gaussian_density(value, mean):
        return -0.5 * (np.log(2.0 * np.pi * variance) + (value - mean) ** 2 / variance)

    def log_volume_weight(lam):
        # S8/S9: the volume factor is lambda * exp(-Pext*lambda^2/kBT).
        # Kinetic and potential energies are accounted for separately in the
        # total effective energy, so they are not accumulated in ebaro.
        return np.log(lam) - external_pressure * lam**2 / kbt

    log_forward = log_gaussian_density(end_lambda, forward_mean)
    log_reverse = log_gaussian_density(start_lambda, reverse_mean)
    expected_from_probabilities = -kbt * (
        log_volume_weight(end_lambda)
        - log_volume_weight(start_lambda)
        + log_reverse
        - log_forward
    )
    assert barostat.ebaro == pytest.approx(
        initial_ebaro + expected_from_probabilities, rel=1e-12, abs=1e-12
    )

    # Expanded form (S11). Expanding the Gaussian ratio above gives
    # D*duration/(4*kBT) multiplying the force-square difference; the last
    # coefficient printed in S11 differs from this consequence of S7.
    expected_correction = barostat.pext * (barostat.cell.V - old_volume)
    expected_correction -= barostat.temp * np.log(np.sqrt(barostat.cell.V) / old_lambda)
    expected_correction += 0.5 * delta_lambda * (old_force + new_force)
    expected_correction += (
        diffusion * coupling_dt / (4.0 * barostat.temp) * (new_force**2 - old_force**2)
    )

    assert barostat.ebaro == pytest.approx(initial_ebaro + expected_correction)
    assert barostat._scr_state is None


def test_scr_rejects_nonpositive_proposed_volume():
    """Invalid stochastic moves fail loudly instead of being clipped."""

    barostat = BaroSCR(
        dt=1.0,
        temp=1.0,
        tau=1.0,
        bulk_modulus=1.0,
        pext=0.0,
    )
    barostat.cell = OrthorhombicCell([1.0, 1.0, 1.0])
    barostat.prng = FixedGaussian(-10.0)
    barostat.get_lambda_force = lambda: 0.0
    barostat.qdt_ = depend_value(name="qdt", value=barostat.dt / 2.0)

    with pytest.raises(ValueError, match="non-positive or non-finite sqrt"):
        barostat.prepare()


def test_dynamics_uses_standard_npt_integrator_for_scr():
    """SCR uses the same NPT integrator as other isotropic barostats."""

    dynamics = Dynamics(
        timestep=0.5,
        mode="npt",
        barostat=BaroSCR(tau=100.0, bulk_modulus=1.0),
    )
    assert type(dynamics.integrator) is NPTIntegrator


@pytest.mark.parametrize("splitting", ["obabo", "baoab"])
def test_scr_uses_standard_npt_splittings(simulation, splitting):
    """Both standard NPT thermostat splittings complete an SCR move."""

    xml = scr_simulation_xml().replace(
        '<dynamics mode="npt">',
        f'<dynamics mode="npt" splitting="{splitting}">',
    )

    with simulation(xml) as sim:
        sim.run(steps=2, write_outputs=False)
        assert type(sim.syslist[0].motion.integrator) is NPTIntegrator
        assert sim.syslist[0].motion.barostat._scr_state is None


def test_scr_uses_standard_multiple_time_stepping(simulation):
    """The shared NPT integrator can apply SCR with multiple force levels."""

    xml = scr_simulation_xml().replace(
        '<force forcefield="harmonic"/>',
        '<force forcefield="harmonic"><mts_weights>[1,1]</mts_weights></force>',
    )
    xml = xml.replace(
        '<timestep units="femtosecond">0.5</timestep>',
        '<timestep units="femtosecond">0.5</timestep><nmts>[1,2]</nmts>',
    )

    with simulation(xml) as sim:
        sim.run(steps=2, write_outputs=False)
        assert sim.syslist[0].motion.barostat._scr_state is None


def test_flexible_barostat_volume_constraint(simulation):
    with simulation(pressure_input("flexible", constant_volume=True)) as sim:
        initial = float(sim.syslist[0].cell.V)
        sim.run(steps=8, write_outputs=False)
        assert float(sim.syslist[0].cell.V) == pytest.approx(initial, rel=1e-12)


def test_bzp_mass_from_bulk_modulus(simulation):
    """The piston mass is 3 N kB (P T) tau^2, or 3 P V B tau^2 with a bulk modulus:
    the same for that of an ideal gas, N kB T / V, but independent of T."""

    root = ET.fromstring(pressure_input("isotropic", nbeads=2))
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        system = sim.syslist[0]
        barostat = system.motion.barostat
        kbt = Constants.kb * float(system.ensemble.temp)
        mass = 3 * system.beads.natoms * 2 * kbt * barostat.tau**2
        assert float(barostat.m[0]) == pytest.approx(mass, rel=1e-12)
        ideal_gas = system.beads.natoms * kbt / float(system.cell.V)
        system.ensemble.temp = 4 * system.ensemble.temp
        assert float(barostat.m[0]) == pytest.approx(4 * mass, rel=1e-12)

    baro = root.find("system/motion/dynamics/barostat")
    ET.SubElement(baro, "bulk_modulus").text = repr(ideal_gas)
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        system = sim.syslist[0]
        barostat = system.motion.barostat
        assert float(barostat.m[0]) == pytest.approx(mass, rel=1e-12)
        system.ensemble.temp = 4 * system.ensemble.temp
        assert float(barostat.m[0]) == pytest.approx(mass, rel=1e-12)


@pytest.mark.parametrize("mode", ["flexible", "anisotropic"])
def test_flexible_cell_mass_from_bulk_modulus(simulation, mode):
    """The mass of the cell is N kB (P T) tau^2, or P V B tau^2 with a bulk
    modulus, with V the volume of the reference cell for 'anisotropic'."""

    root = ET.fromstring(pressure_input(mode, nbeads=2))
    baro = root.find("system/motion/dynamics/barostat")
    if mode == "anisotropic":
        baro.find("h0").text = str((np.eye(3) * 21).ravel().tolist())
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        system = sim.syslist[0]
        barostat = system.motion.barostat
        kbt = Constants.kb * float(system.ensemble.temp)
        mass = system.beads.natoms * 2 * kbt * barostat.tau**2
        assert float(barostat.m[0]) == pytest.approx(mass, rel=1e-12)

    ET.SubElement(baro, "bulk_modulus").text = "0.5"
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        barostat = sim.syslist[0].motion.barostat
        volume = 21.0**3 if mode == "anisotropic" else 20.0**3
        assert barostat.v0 == pytest.approx(volume, rel=1e-12)
        mass = 2 * volume * 0.5 * barostat.tau**2
        assert float(barostat.m[0]) == pytest.approx(mass, rel=1e-12)


@pytest.mark.parametrize("mode", ["isotropic", "flexible"])
def test_restart_keeps_the_reference_volume(restart_states, mode):
    """The reference volume of the piston mass is the first volume, also after a
    restart, which continues the trajectory of the continuous run."""

    root = ET.fromstring(pressure_input(mode))
    root.find("system/motion/dynamics/timestep").text = "0.5"
    baro = root.find("system/motion/dynamics/barostat")
    ET.SubElement(baro, "bulk_modulus").text = "1e-6"

    def snapshot(sim):
        barostat = sim.syslist[0].motion.barostat
        return simulation_state(sim) + (barostat.v0, float(barostat.m[0]))

    continuous, restarted = restart_states(
        ET.tostring(root, encoding="unicode"), snapshot, steps=(6, 4)
    )

    assert continuous[5] == 8000.0
    assert restarted[5] == continuous[5]
    assert np.linalg.det(continuous[3]) != pytest.approx(8000.0, rel=1e-6)
    assert restarted[6] == continuous[6]
    assert restarted[0] == continuous[0]
    for restarted_value, continuous_value in zip(restarted[1:5], continuous[1:5]):
        np.testing.assert_allclose(
            restarted_value, continuous_value, rtol=1e-12, atol=1e-12
        )


@pytest.mark.parametrize("barostat", [BaroBZP, BaroMTK, BaroRGB])
def test_bulk_modulus_input(barostat):
    """The bulk modulus and the reference volume of the piston mass are stored
    only when the bulk modulus is given (the reference cell sets the volume of
    'anisotropic')."""

    options = {"h0": Cell(np.eye(3) * 2.0)} if barostat is BaroRGB else {"v0": 10.0}
    stored = InputBaro()
    stored.store(barostat(tau=100.0, **options))
    assert not stored.bulk_modulus._explicit and not stored.v0._explicit
    stored = InputBaro()
    stored.store(barostat(tau=100.0, bulk_modulus=0.5, **options))
    restored = stored.fetch()
    assert type(restored) is barostat
    assert restored.bulk_modulus == pytest.approx(0.5)
    if barostat is not BaroRGB:
        assert restored.v0 == pytest.approx(10.0)


def test_bulk_modulus_and_reference_volume_are_rejected_where_unused():
    input_baro = InputBaro()
    input_baro._explicit = True
    input_baro.mode.store("dummy")
    input_baro.bulk_modulus.store(1.0)
    with pytest.raises(ValueError, match="bulk modulus is not used"):
        input_baro.fetch()

    for mode, bulk_modulus in (
        ("isotropic", None),
        ("flexible", None),
        ("stochastic-rescaling", 1.0),
        ("anisotropic", 1.0),
    ):
        input_baro = InputBaro()
        input_baro._explicit = True
        input_baro.mode.store(mode)
        if bulk_modulus is not None:
            input_baro.bulk_modulus.store(bulk_modulus)
        input_baro.v0.store(10.0)
        with pytest.raises(ValueError, match="v0 is only used"):
            input_baro.fetch()


@pytest.mark.parametrize(
    "hfix,expected",
    [
        ([], [3, 1, -1]),
        (["offdiagonal"], [1, 1, 1]),
        (["offdiagonal", "xx", "yy"], [0, 0, 1]),
        (["xz", "yz"], [2, 0, 1]),
        (["yy", "yz"], [3, 0, 0]),
        (["zz"], [3, 1, 0]),
    ],
)
def test_measure_term_of_the_cell(hfix, expected):
    """L_kk is one plus the free elements of row k minus those of column k, for
    the free diagonal elements (zero for the fixed ones)."""

    barostat = BaroMTK(hfix=hfix)
    np.testing.assert_array_equal(np.diag(to_numpy(barostat.L)), expected)


def test_fixing_xz_alone_is_rejected():
    with pytest.raises(ValueError, match="xz requires fixing xy or yz"):
        BaroMTK(hfix=["xz"])


CELL_METRIC = """import numpy as np
from ipi.pes.dummy import Dummy_driver

__DRIVER_CLASS__ = "CellMetric"

class CellMetric(Dummy_driver):
    def compute_structure(self, cell, pos):
        dg = cell.T @ cell - 100.0 * np.eye(3)
        u, f = 0.005 * np.sum(dg**2), np.zeros_like(pos)
        vir = -0.02 * cell @ dg @ cell.T
        # harmonic bonds of length 1.5 between all the pairs of atoms
        for i in range(len(pos)):
            for j in range(i):
                r = pos[i] - pos[j]
                d = np.linalg.norm(r)
                fij = -0.5 * (d - 1.5) * r / d
                u += 0.25 * (d - 1.5) ** 2
                f[i] += fij
                f[j] -= fij
                vir += np.outer(fij, r)
        return u, f, vir, ""
"""


def test_conserved_quantity_with_fixed_cell_elements(simulation, tmp_path):
    """Without thermostats, the conserved quantity of a flexible cell with the
    off-diagonal elements fixed is conserved: its measure term follows hfix.
    One free particle in a cell with U = (k/2) |h^T h - a0^2 I|^2."""

    (tmp_path / "cell_metric.py").write_text(CELL_METRIC)
    xml = f"""<simulation verbosity="quiet">
      <total_steps>1000</total_steps>
      <ffdirect name="metric"><pes>custom</pes><parameters>{{}}</parameters>
        <pes_path>{tmp_path / "cell_metric.py"}</pes_path></ffdirect>
      <system><forces><force forcefield="metric"/></forces>
        <ensemble><temperature>1.0</temperature><pressure>0.0</pressure></ensemble>
        <motion mode="dynamics"><fixcom>False</fixcom>
          <dynamics mode="npt"><timestep>0.05</timestep>
            <barostat mode="flexible"><tau>20</tau><hfix>[offdiagonal]</hfix>
              <p shape="(3,3)">[4, 0, 0, 0, -6, 0, 0, 0, 2]</p></barostat>
          </dynamics></motion>
        <beads natoms="1" nbeads="1"><q shape="(1,3)">[1, 2, 3]</q>
          <p shape="(1,3)">[0.5, -0.3, 0.2]</p><m shape="(1)">[1.0]</m>
          <names shape="(1)">[H]</names></beads>
        <cell shape="(3,3)">[10, 0, 0, 0, 10, 0, 0, 0, 10]</cell>
      </system></simulation>"""

    with simulation(xml) as sim:
        system = sim.syslist[0]
        econs, diagonal = [], []
        for _ in range(200):
            sim.run(steps=1, write_outputs=False)
            econs.append(float(system.ensemble.econs))
            diagonal.append(np.log(np.diag(to_numpy(system.cell.h))))
    assert np.ptp(diagonal, axis=0).min() > 0.01
    assert np.ptp(econs) < 1e-3


@pytest.mark.parametrize(
    "hfix,nbeads", [("[]", 2), ("[offdiagonal, xx, yy]", 1)], ids=["full", "npzt"]
)
def test_conserved_quantity_converges_as_dt_squared(simulation, tmp_path, hfix, nbeads):
    """Without thermostats, the error on the conserved quantity of a flexible
    cell decreases as dt^2: the barostat feels the off-diagonal kinetic stress
    and centroid virial that a shear cell velocity exchanges with the beads,
    and the fixed elements of the cell receive no momentum. Three bonded atoms
    in a cell with U = (k/2) |h^T h - a0^2 I|^2."""

    (tmp_path / "cell_metric.py").write_text(CELL_METRIC)
    rng = np.random.default_rng(3)
    q = np.tile([1.0, 2.0, 3.0, 2.5, 2.0, 3.0, 1.7, 3.3, 3.2], (nbeads, 1))
    if nbeads > 1:
        q += 0.1 * rng.standard_normal(q.shape)
    p = rng.standard_normal(q.shape)

    def spread(dt):
        xml = f"""<simulation verbosity="quiet">
          <total_steps>1000</total_steps>
          <ffdirect name="metric"><pes>custom</pes><parameters>{{}}</parameters>
            <pes_path>{tmp_path / "cell_metric.py"}</pes_path></ffdirect>
          <system><forces><force forcefield="metric"/></forces>
            <ensemble><temperature>1.0</temperature><pressure>0.0</pressure></ensemble>
            <motion mode="dynamics"><fixcom>False</fixcom>
              <dynamics mode="npt"><timestep>{dt}</timestep>
                <barostat mode="flexible"><tau>10</tau><hfix>{hfix}</hfix>
                  <p shape="(3,3)">[2, 10, 10, 0, -3, 10, 0, 0, 1]</p></barostat>
              </dynamics></motion>
            <beads natoms="3" nbeads="{nbeads}">
              <q shape="({nbeads},9)">{q.ravel().tolist()}</q>
              <p shape="({nbeads},9)">{p.ravel().tolist()}</p><m shape="(3)">[1, 1, 1]</m>
              <names shape="(3)">[H, H, H]</names></beads>
            <cell shape="(3,3)">[10, 0, 0, 0, 10, 0, 0, 0, 10]</cell>
          </system></simulation>"""
        with simulation(xml) as sim:
            system = sim.syslist[0]
            econs = []
            for _ in range(round(8.0 / dt)):
                sim.run(steps=1, write_outputs=False)
                econs.append(float(system.ensemble.econs))
        return np.ptp(econs)

    coarse, fine = spread(0.04), spread(0.02)
    assert fine < 5e-3
    assert coarse / fine > 3.0


def test_cell_propagators_match_power_series():
    """exp(+-v dt) and dt sinh(v dt)/(v dt) are correct both for cell velocities
    with zero diagonal (nilpotent) and with non-zero diagonal"""

    dt = 2.0
    shear = np.array([[0.0, 0.3, -0.2], [0.0, 0.0, 0.1], [0.0, 0.0, 0.0]])
    # the scaling and squaring of the matrix exponential loses a few digits
    tolerance = dict(rtol=1e-12, atol=1e-12)
    for v in (shear, shear + np.diag([1e-2, -2e-2, 3e-2])):
        x = v * dt
        powers = [np.linalg.matrix_power(x, k) / math.factorial(k) for k in range(25)]
        expq, expp, sinh = cell_propagators(v, dt)
        assert np.allclose(expq, sum(powers), **tolerance)
        exp_minus = sum(c * (-1) ** k for k, c in enumerate(powers))
        assert np.allclose(expp, exp_minus, **tolerance)
        series = sum(powers[2 * k] / (2 * k + 1) for k in range(12))
        assert np.allclose(sinh, dt * series, **tolerance)

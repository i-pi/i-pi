"""Tests for pressure-control algorithms."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

import math
import xml.etree.ElementTree as ET

import numpy as np
import pytest

from ipi.engine.barostats import BaroSCR, dt_sinch
from ipi.engine.motion.dynamics import Dynamics, NPTIntegrator
from ipi.engine.thermostats import ThermoLangevin
from ipi.inputs.barostats import InputBaro
from ipi.utils.depend import depend_array, depend_value, dproperties
from ipi.utils.units import unit_to_internal
from ipi_tests.unit_tests.engine.inputs import pressure_input


class MockNormalModes:
    """Normal-mode positions and momenta, with the masses of the modes."""

    def __init__(self, qnm, pnm, dynm3):
        self.qnm_ = depend_array(name="qnm", value=qnm)
        self.pnm_ = depend_array(name="pnm", value=pnm)
        self.dynm3 = dynm3


dproperties(MockNormalModes, ["qnm", "pnm"])


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
          <compressibility units="bar^-1">4.5e-5</compressibility>
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
        np.asarray(system.beads.q).copy(),
        np.asarray(system.beads.p).copy(),
        np.asarray(system.cell.h).copy(),
        system.motion.barostat.ebaro,
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
        ("system/motion/dynamics/barostat/compressibility", "0.3"),
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
        assert system.nm.kin == pytest.approx(3.0, rel=1e-12)
        np.testing.assert_allclose(system.forces.vir, np.zeros((3, 3)), atol=1e-12)
        assert np.trace(barostat.stress_mts(0)) / 3.0 == pytest.approx(0.25, rel=1e-12)
        assert barostat.get_lambda_force() == pytest.approx(
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
        assert system.nm.kin == pytest.approx(3.0, rel=1e-12)
        assert system.forces.pot == pytest.approx(1.0, rel=1e-12)
        np.testing.assert_allclose(
            system.forces.vir, np.diag([-2.0, 0.0, 0.0]), rtol=1e-12, atol=1e-12
        )
        assert np.trace(barostat.stress_mts(0)) / 3.0 == pytest.approx(
            1.0 / 6.0, rel=1e-12
        )
        assert barostat.get_lambda_force() == pytest.approx(
            force_coefficient * np.sqrt(8.0), rel=1e-12, abs=1e-12
        )


def test_scr_input_roundtrip_preserves_restart_state():
    """SCR parameters and accumulated work survive checkpoint storage."""

    compressibility = unit_to_internal("inverse-pressure", "bar^-1", 4.5e-5)
    original = BaroSCR(
        tau=100.0,
        compressibility=compressibility,
        ebaro=3.25,
    )

    stored = InputBaro()
    stored.store(original)
    restored = stored.fetch()

    assert type(restored) is BaroSCR
    assert restored.tau == pytest.approx(original.tau)
    assert restored.compressibility == pytest.approx(original.compressibility)
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


def test_scr_requires_explicit_positive_compressibility():
    """SCR cannot silently guess a material compressibility."""

    input_baro = InputBaro()
    input_baro._explicit = True
    input_baro.mode.store("stochastic-rescaling")
    input_baro.tau.store(100.0)

    with pytest.raises(ValueError, match="must be specified explicitly"):
        input_baro.fetch()

    input_baro.compressibility.store(-1.0)
    with pytest.raises(ValueError, match="compressibility must be positive"):
        input_baro.fetch()


def test_scr_rejects_a_separate_cell_thermostat():
    """SCR has no piston degree of freedom to thermostat."""

    input_baro = InputBaro()
    input_baro._explicit = True
    input_baro.mode.store("stochastic-rescaling")
    input_baro.tau.store(100.0)
    input_baro.compressibility.store(1.0)
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
    compressibility = 0.3
    external_pressure = 0.7
    initial_ebaro = 1.25
    gaussian_draw = 0.4
    barostat = BaroSCR(
        dt=duration,
        temp=kbt,
        tau=relaxation_time,
        compressibility=compressibility,
        pext=external_pressure,
        ebaro=initial_ebaro,
    )
    barostat.prng = FixedGaussian(gaussian_draw)
    barostat.cell = OrthorhombicCell([2.0, 2.0, 2.0])
    positions = np.asarray([[0.2, -0.3, 0.4, 0.5, -0.1, 0.7]])
    momenta = np.asarray([[0.8, -0.4, 0.2, 0.3, 0.5, -0.6]])
    masses = np.asarray([[1.0, 1.0, 1.0, 2.0, 2.0, 2.0]])
    barostat.nm = MockNormalModes(positions.copy(), momenta.copy(), masses)
    barostat.qdt_ = depend_value(name="qdt", value=barostat.dt / 2.0)

    old_volume = barostat.cell.V
    old_lambda = np.sqrt(old_volume)
    old_force = 0.6
    new_force = -0.2
    forces = iter((old_force, new_force))
    barostat.get_lambda_force = lambda: next(forces)

    coupling_dt = 2.0 * barostat.qdt
    diffusion = barostat.temp * barostat.compressibility / (4.0 * barostat.tau)
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
        barostat.nm.qnm[0],
        scale * positions[0] + drift_scale * momenta[0] * barostat.dt / masses[0],
    )
    np.testing.assert_allclose(barostat.nm.pnm[0], inverse_scale * momenta[0])
    assert barostat.cell.V == pytest.approx(old_volume * scale**3)

    # Bernetti and Bussi (2020), Supplementary Methods, Eqs. S7--S9:
    # https://arxiv.org/src/2006.09250v4/anc/si.pdf
    # Calculate Gaussian transition densities from S7, without expanding their
    # ratio into the force sums and force squares used by finalize(). All
    # reference inputs below are prescribed, including the initial volume 8.
    variance = kbt * compressibility * duration / (2.0 * relaxation_time)
    drift_per_force = compressibility * duration / (4.0 * relaxation_time)
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
        compressibility=1.0,
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
        barostat=BaroSCR(tau=100.0, compressibility=1.0),
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
        initial = sim.syslist[0].cell.V
        sim.run(steps=8, write_outputs=False)
        assert sim.syslist[0].cell.V == pytest.approx(initial, rel=1e-12)


def test_dt_sinch_matches_power_series():
    """dt sinh(v dt)/(v dt) is correct both for cell velocities with zero
    diagonal (Taylor branch) and with non-zero diagonal (eigendecomposition)"""

    dt = 2.0
    shear = np.array([[0.0, 0.3, -0.2], [0.0, 0.0, 0.1], [0.0, 0.0, 0.0]])
    for v in (shear, shear + np.diag([1e-2, -2e-2, 3e-2])):
        x = v * dt
        series = sum(
            np.linalg.matrix_power(x, 2 * k) / math.factorial(2 * k + 1)
            for k in range(12)
        )
        assert np.allclose(dt_sinch(v, dt), dt * series, rtol=1e-12, atol=1e-14)

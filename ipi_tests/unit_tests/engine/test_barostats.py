"""Tests for pressure-control algorithms."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

import io
import os
from types import SimpleNamespace

import numpy as np
import pytest

from ipi.engine.barostats import BaroSCR
from ipi.engine.motion.dynamics import Dynamics, NPTIntegrator
from ipi.engine.thermostats import ThermoLangevin
from ipi.inputs.barostats import InputBaro
from ipi.scripting import InteractiveSimulation
from ipi.utils.depend import depend_value
from ipi.utils.softexit import softexit
from ipi.utils.units import unit_to_internal


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
        self.h = np.diag(np.asarray(lengths, float))

    @property
    def V(self):
        return np.linalg.det(self.h)


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


def simulation_state(simulation):
    """Returns the state needed to compare continuous and restarted runs."""

    system = simulation.syslist[0]
    return (
        simulation.step,
        np.asarray(system.beads.q).copy(),
        np.asarray(system.beads.p).copy(),
        np.asarray(system.cell.h).copy(),
        system.motion.barostat.ebaro,
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


def test_scr_restart_matches_continuous_trajectory(tmp_path):
    """Checkpointing preserves SCR work, phase, and random-number state."""

    original_directory = os.getcwd()
    continuous_directory = tmp_path / "continuous"
    restarted_directory = tmp_path / "restarted"
    continuous_directory.mkdir()
    restarted_directory.mkdir()

    try:
        os.chdir(continuous_directory)
        continuous = InteractiveSimulation(io.StringIO(scr_simulation_xml()))
        continuous.run(steps=10, write_outputs=False)
        continuous_state = simulation_state(continuous)
        continuous.stop()
        softexit.reset()

        os.chdir(restarted_directory)
        interrupted = InteractiveSimulation(io.StringIO(scr_simulation_xml()))
        interrupted.run(steps=6, write_outputs=False)
        interrupted.stop()
        softexit.reset()

        with open("RESTART") as restart_file:
            restarted = InteractiveSimulation(restart_file)
        restarted.run(steps=10 - restarted.step, write_outputs=False)
        restarted_state = simulation_state(restarted)
        restarted.stop()
        softexit.reset()
    finally:
        os.chdir(original_directory)
        softexit.reset()

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


def test_scr_rejects_path_integral_dynamics():
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

    try:
        with pytest.raises(ValueError, match="classical dynamics only"):
            InteractiveSimulation(io.StringIO(xml))
    finally:
        softexit.reset()


def test_scr_trotter_map_and_effective_energy():
    """The implementation follows Eqs. S7, S12, and S11 of the paper."""

    barostat = BaroSCR(
        dt=0.2,
        temp=2.0,
        tau=4.0,
        compressibility=0.3,
        pext=0.7,
        ebaro=1.25,
    )
    barostat.prng = FixedGaussian(0.4)
    barostat.cell = OrthorhombicCell([2.0, 2.0, 2.0])
    positions = np.asarray([[0.2, -0.3, 0.4, 0.5, -0.1, 0.7]])
    momenta = np.asarray([[0.8, -0.4, 0.2, 0.3, 0.5, -0.6]])
    masses = np.asarray([[1.0, 1.0, 1.0, 2.0, 2.0, 2.0]])
    barostat.nm = SimpleNamespace(
        qnm=positions.copy(),
        pnm=momenta.copy(),
        dynm3=masses,
    )
    barostat._qdt = depend_value(name="qdt", value=barostat.dt / 2.0)

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

    expected_correction = barostat.pext * (barostat.cell.V - old_volume)
    expected_correction -= barostat.temp * np.log(np.sqrt(barostat.cell.V) / old_lambda)
    expected_correction += 0.5 * delta_lambda * (old_force + new_force)
    expected_correction += (
        diffusion * coupling_dt / (4.0 * barostat.temp) * (new_force**2 - old_force**2)
    )

    assert barostat.ebaro == pytest.approx(1.25 + expected_correction)
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
    barostat._qdt = depend_value(name="qdt", value=barostat.dt / 2.0)

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
def test_scr_uses_standard_npt_splittings(splitting):
    """Both standard NPT thermostat splittings complete an SCR move."""

    xml = scr_simulation_xml().replace(
        '<dynamics mode="npt">',
        f'<dynamics mode="npt" splitting="{splitting}">',
    )

    try:
        simulation = InteractiveSimulation(io.StringIO(xml))
        simulation.run(steps=2, write_outputs=False)
        assert type(simulation.syslist[0].motion.integrator) is NPTIntegrator
        assert simulation.syslist[0].motion.barostat._scr_state is None
        simulation.stop()
    finally:
        softexit.reset()


def test_scr_uses_standard_multiple_time_stepping():
    """The shared NPT integrator can apply SCR with multiple force levels."""

    xml = scr_simulation_xml().replace(
        '<force forcefield="harmonic"/>',
        '<force forcefield="harmonic"><mts_weights>[1,1]</mts_weights></force>',
    )
    xml = xml.replace(
        '<timestep units="femtosecond">0.5</timestep>',
        '<timestep units="femtosecond">0.5</timestep><nmts>[1,2]</nmts>',
    )

    try:
        simulation = InteractiveSimulation(io.StringIO(xml))
        simulation.run(steps=2, write_outputs=False)
        assert simulation.syslist[0].motion.barostat._scr_state is None
        simulation.stop()
    finally:
        softexit.reset()

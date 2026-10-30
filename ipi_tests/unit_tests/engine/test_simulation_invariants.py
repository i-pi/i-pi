"""Small simulations checked against analytic motion and initialization rules."""

import copy
import io
import xml.etree.ElementTree as ET

import numpy as np
import pytest

from ipi.engine.simulation import Simulation
from ipi_tests.unit_tests.engine.inputs import harmonic_input


@pytest.mark.parametrize("threading", [False, True])
def test_nve_tracks_exact_harmonic_motion(simulation, threading):
    """Velocity Verlet agrees with the analytic oscillator at a finite time."""
    with simulation(harmonic_input(threading=threading, timestep=0.01)) as sim:
        beads = sim.syslist[0].beads
        q0, p0, masses = beads.q.copy(), beads.p.copy(), beads.m3.copy()
        frequency = np.sqrt(0.5 / masses)
        sim.run(steps=20, write_outputs=False)
        time = 0.2
        expected_q = q0 * np.cos(frequency * time) + p0 / (masses * frequency) * np.sin(
            frequency * time
        )
        expected_p = p0 * np.cos(frequency * time) - masses * frequency * q0 * np.sin(
            frequency * time
        )
        np.testing.assert_allclose(beads.q, expected_q, rtol=0, atol=5e-7)
        np.testing.assert_allclose(beads.p, expected_p, rtol=0, atol=5e-7)
        assert sim.step == 20
        assert sim.syslist[0].ensemble.time == pytest.approx(time)


def test_nve_with_fixed_atoms_tracks_harmonic_motion(simulation):
    """With a fixed atom, the free atom follows the analytic oscillator."""
    xml = harmonic_input(timestep=0.01).replace(
        "<fixcom>False</fixcom>", "<fixcom>False</fixcom><fixatoms>[1]</fixatoms>"
    )
    with simulation(xml) as sim:
        beads = sim.syslist[0].beads
        q0, p0, masses = beads.q.copy(), beads.p.copy(), beads.m3.copy()
        frequency = np.sqrt(0.5 / masses)
        sim.run(steps=20, write_outputs=False)
        time = 0.2
        expected_q = q0 * np.cos(frequency * time) + p0 / (masses * frequency) * np.sin(
            frequency * time
        )
        np.testing.assert_allclose(beads.q[:, :3], expected_q[:, :3], atol=5e-7)
        np.testing.assert_array_equal(beads.q[:, 3:], q0[:, 3:])


def test_integrator_constants_follow_parameters(simulation):
    """The time steps of the integrator follow changes of the time step, of the
    masses and of the multiple time steps."""
    with simulation(harmonic_input()) as sim:
        motion, beads = sim.syslist[0].motion, sim.syslist[0].beads
        integrator = motion.integrator
        qdt = integrator.qdt
        np.testing.assert_allclose(integrator.qdt_on_m, qdt / beads.m3[0])
        motion.dt = 2 * motion.dt
        beads.m = 3 * beads.m
        np.testing.assert_allclose(integrator.qdt_on_m, 2 * qdt / beads.m3[0])
        motion.nmts = [2]
        assert integrator.inmts == 2 and integrator.qdt == pytest.approx(qdt)


@pytest.mark.parametrize("nbeads", [1, 4])
def test_nve_time_reversal(simulation, nbeads):
    """Reversing momenta retraces classical and ring-polymer trajectories."""
    with simulation(harmonic_input(nbeads=nbeads)) as sim:
        beads = sim.syslist[0].beads
        q0, p0 = beads.q.copy(), beads.p.copy()
        sim.run(steps=7, write_outputs=False)
        beads.p_[:] *= -1
        sim.run(steps=7, write_outputs=False)
        np.testing.assert_allclose(beads.q, q0, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(beads.p, -p0, rtol=1e-12, atol=1e-12)


def test_com_removal_from_beads_or_normal_modes(simulation):
    """The motion of the centre of mass is removed in the same way from the
    momenta of the beads and from those of the normal modes."""
    xml = harmonic_input(nbeads=4).replace("<fixcom>False", "<fixcom>True")
    with simulation(xml) as sim:
        s = sim.syslist[0]
        p = np.random.default_rng(0).normal(size=s.beads.p.shape)
        results = []
        for in_normal_modes in (False, True):
            if in_normal_modes:
                s.nm.pnm = s.nm.transform.b2nm(p)
            else:
                s.beads.p = p
            assert s.beads.p_.tainted() == in_normal_modes
            eens = s.ensemble.eens
            s.motion.integrator.pconstraints()
            results.append((s.beads.p.copy(), s.ensemble.eens - eens))
        np.testing.assert_allclose(results[1][0], results[0][0], atol=1e-14)
        assert results[1][1] == pytest.approx(results[0][1], rel=1e-12)
        np.testing.assert_allclose(
            results[0][0].reshape(-1, 3).sum(axis=0), 0, atol=1e-14
        )


@pytest.mark.parametrize("threading", [False, True])
def test_output_stride_matches_completed_steps(simulation, threading):
    """Output rows hold the step and energies of the state at that step."""
    xml = harmonic_input(threading=threading).replace(
        "<system>",
        '<output prefix="sample"><properties filename="properties" stride="2">'
        "[step, potential, kinetic_md]</properties></output><system>",
    )
    expected_energies = []
    with simulation(xml) as sim:
        system = sim.syslist[0]
        for completed_steps in [0, 2, 4]:
            if completed_steps:
                sim.run(steps=2)
            expected_energies.append(
                [
                    0.25 * np.sum(system.beads.q**2),
                    np.sum(system.beads.p**2 / (2 * system.beads.m3)),
                ]
            )
    rows = np.loadtxt("sample.properties")
    np.testing.assert_array_equal(rows[:, 0], [0, 2, 4])
    np.testing.assert_allclose(rows[:, 1:], expected_energies, rtol=1e-12)


@pytest.mark.xfail(
    strict=True,
    raises=ValueError,
    reason="BUG: GLE state initialization uses ensemble.thermostat instead of motion.thermostat",
)
def test_gle_auxiliary_momenta_initialized_from_input(simulation):
    xml = harmonic_input(thermostat="gle").replace(
        "<tau>0.5</tau>", "<A shape='(2,2)'>[1,0,0,1]</A>"
    )
    xml = xml.replace(
        "<system>",
        "<system><initialize nbeads='1'><gle mode='manual'>[1,2,3,4,5,6,7,8,9,10,11,12]</gle></initialize>",
    )
    with simulation(xml) as sim:
        np.testing.assert_array_equal(
            sim.syslist[0].motion.thermostat.s.ravel(), np.arange(1, 13)
        )


@pytest.mark.parametrize("source", ["string", "file"])
def test_read_only_initialization_starts_no_services(tmp_path, monkeypatch, source):
    monkeypatch.chdir(tmp_path)
    xml = harmonic_input()
    observed = []

    def initialize(sim):
        observed.append(sim.chk)
        sim.syslist[0].beads.q_[:] = 2

    sim = Simulation.load_from_xml(
        xml if source == "string" else io.StringIO(xml),
        read_only=True,
        post_init_hook=initialize,
        custom_verbosity="quiet",
        sockets_prefix="/tmp/initialization-test-",
    )
    assert observed == [None]
    assert sim.chk is None
    assert sim.sockets_prefix == "/tmp/initialization-test-"
    assert not list(tmp_path.iterdir())
    assert all(field._thread is None for field in sim.fflist.values())
    np.testing.assert_array_equal(sim.syslist[0].beads.q, 2)
    # Forces are bound to the modified coordinates even in read-only mode.
    assert sim.syslist[0].forces.pot == pytest.approx(6.0)
    np.testing.assert_allclose(sim.syslist[0].forces.f, -1.0)


@pytest.mark.parametrize(
    "xml,message",
    [
        ("not an input file", "neither valid JSON nor valid XML"),
        (harmonic_input().replace("<total_steps>100", "<total_steps>0"), "already run"),
    ],
)
def test_invalid_or_finished_simulation_does_not_start(
    tmp_path, monkeypatch, xml, message
):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError, match=message):
        Simulation.load_from_xml(xml, read_only=True)
    assert not list(tmp_path.iterdir())


def test_multiple_systems_have_independent_state_and_outputs(simulation):
    root = ET.fromstring(harmonic_input())
    first = root.find("system")
    first.set("prefix", "first")
    second = copy.deepcopy(first)
    second.set("prefix", "second")
    second.find("beads/q").text = "[1,1,1,1,1,1]"
    root.append(second)
    output = ET.SubElement(root, "output", prefix="initial")
    ET.SubElement(output, "properties", filename="energy").text = "[potential]"
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        first, second = sim.syslist
        assert first.beads is not second.beads
        assert first.forces.pot == pytest.approx(0.25 * 0.43)
        assert second.forces.pot == pytest.approx(1.5)
        first.beads.q_[:] = 0
        np.testing.assert_array_equal(second.beads.q, 1)
        assert {out.filename for out in sim.outputs} == {
            "first_initial.energy",
            "second_initial.energy",
        }

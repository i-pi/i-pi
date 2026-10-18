"""Small real simulations with analytic and restart consistency checks."""

import numpy as np
import pytest


def harmonic_input(nbeads=1, thermostat=None, threading=False, timestep=0.02):
    """Build a two-atom simulation in atomic units with a known quadratic potential."""
    q = np.tile([0.2, -0.3, 0.1, -0.4, 0.2, 0.3], (nbeads, 1))
    q[:, 0] += np.arange(nbeads) * 0.01
    p = np.tile([0.1, 0.2, -0.1, 0.3, -0.2, 0.1], (nbeads, 1))
    bath = (
        ""
        if thermostat is None
        else f'<thermostat mode="{thermostat}"><tau>0.5</tau></thermostat>'
    )
    mode = "nve" if thermostat is None else "nvt"
    return f"""<simulation verbosity="quiet" threading="{threading}" floatformat="%24.16e">
      <total_steps>100</total_steps><prng><seed>317</seed></prng>
      <ffdirect name="harmonic"><pes>harmonic</pes><parameters>{{k1: 0.5}}</parameters></ffdirect>
      <system><forces><force forcefield="harmonic"/></forces>
        <ensemble><temperature units="kelvin">300</temperature></ensemble>
        <motion mode="dynamics"><fixcom>False</fixcom>
          <dynamics mode="{mode}"><timestep>{timestep}</timestep>{bath}</dynamics>
        </motion>
        <beads natoms="2" nbeads="{nbeads}">
          <q shape="({nbeads},6)">{q.ravel().tolist()}</q>
          <p shape="({nbeads},6)">{p.ravel().tolist()}</p>
          <m shape="(2)">[2, 3]</m><names shape="(2)">[H, O]</names>
        </beads><cell shape="(3,3)">[20,0,0,0,20,0,0,0,20]</cell>
      </system></simulation>"""


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


@pytest.mark.parametrize("nbeads", [1, 4])
@pytest.mark.parametrize("thermostat", [None, "langevin", "svr", "pile_l", "pile_g"])
def test_restart_preserves_next_positions_momenta_and_heat(
    restart_states, nbeads, thermostat
):
    """Continuous and restarted dynamics must use the same state and random draws."""
    xml = harmonic_input(nbeads=nbeads, thermostat=thermostat)

    def snapshot(sim):
        system = sim.syslist[0]
        return dict(
            q=system.beads.q.copy(),
            p=system.beads.p.copy(),
            heat=system.motion.thermostat.ethermo,
            time=system.ensemble.time,
            step=sim.step,
        )

    expected, actual = restart_states(xml, snapshot, steps=(3, 5))
    for key in ["q", "p"]:
        np.testing.assert_allclose(actual[key], expected[key], rtol=1e-12, atol=1e-12)
    assert actual["heat"] == pytest.approx(expected["heat"], rel=1e-12, abs=1e-12)
    assert actual["time"] == pytest.approx(expected["time"])
    assert actual["step"] == 8


@pytest.mark.parametrize("nbeads", [1, 4])
def test_nve_time_reversal_recovers_initial_state(simulation, nbeads):
    """Reversing momenta retraces classical and ring-polymer trajectories."""
    with simulation(harmonic_input(nbeads=nbeads)) as sim:
        beads = sim.syslist[0].beads
        q0, p0 = beads.q.copy(), beads.p.copy()
        sim.run(steps=7, write_outputs=False)
        beads.p[:] *= -1
        sim.run(steps=7, write_outputs=False)
        np.testing.assert_allclose(beads.q, q0, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(beads.p, -p0, rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize("nbeads", [1, 4])
@pytest.mark.parametrize(
    "selection,indices", [("", [0, 1]), ("H", [0]), ("O", [1]), ("0", [0]), ("1", [1])]
)
def test_energy_estimators_respect_atom_selection(
    simulation, nbeads, selection, indices
):
    """Known momenta and a harmonic bead spread give independent kinetic estimates."""
    from ipi.utils.units import Constants, unit_to_internal

    with simulation(harmonic_input(nbeads=nbeads)) as sim:
        props = sim.syslist[0].properties
        # p_H=(.1,.2,-.1), m_H=2; p_O=(.3,-.2,.1), m_O=3.
        kinetic_per_atom = [0.06 / 4, 0.14 / 6]
        expected_kinetic = sum(kinetic_per_atom[index] for index in indices)
        kinetic, count = props.get_kinmd(atom=selection, return_count=True)
        assert count == len(indices)
        assert kinetic == pytest.approx(expected_kinetic)
        assert props.get_kinmd(atom=selection, bead="0") == pytest.approx(
            expected_kinetic
        )
        # Only H is spread: Var(0,.01,...,(P-1)*.01) = .01²(P²-1)/12.
        virial_correction = 0.25 * 0.01**2 * (nbeads**2 - 1) / 12 if 0 in indices else 0
        expected_cv = (
            len(indices)
            * 1.5
            * Constants.kb
            * unit_to_internal("temperature", "kelvin", 300)
            + virial_correction
        )
        assert props.get_kincv(selection) == pytest.approx(expected_cv, abs=1e-12)


def test_properties_reject_out_of_range_and_conflicting_selections(simulation):
    """Invalid atom/bead/mode selectors must not produce a plausible energy."""
    with simulation(harmonic_input(nbeads=4)) as sim:
        props = sim.syslist[0].properties
        with pytest.raises(IndexError, match="atom index"):
            props.get_kinmd(atom="2")
        with pytest.raises(IndexError, match="Bead index"):
            props.get_kinmd(bead="4")
        with pytest.raises(IndexError, match="Normal mode index"):
            props.get_kinmd(nm="4")
        with pytest.raises(ValueError, match="both NM and bead"):
            props.get_kinmd(bead="0", nm="0")
        with pytest.raises(ValueError, match="valid integer"):
            props.get_kinmd(bead="first")


@pytest.mark.parametrize(
    "method", ["sd", "cg", "bfgs", "bfgstrm", "lbfgs", "damped_bfgs"]
)
@pytest.mark.parametrize("selection", [None, "fixatoms", "fixatoms_dof"])
def test_geometry_optimization_finds_known_harmonic_minimum(
    simulation, method, selection, request
):
    """Minimize free coordinates while preserving any fixed atoms or coordinates."""
    import xml.etree.ElementTree as ET

    if method == "damped_bfgs" and selection == "fixatoms_dof":
        request.node.add_marker(
            pytest.mark.xfail(
                strict=True,
                raises=ValueError,
                reason="Damped_BFGS reshapes five free coordinates into triples",
            )
        )
    root = ET.fromstring(harmonic_input())
    motion = root.find("system/motion")
    motion.clear()
    motion.set("mode", "minimize")
    ET.SubElement(motion, "fixcom").text = "False"
    if selection:
        ET.SubElement(motion, selection, shape="(1)").text = "[0]"
    optimizer = ET.SubElement(motion, "optimizer", mode=method)
    ET.SubElement(optimizer, "exit_on_convergence").text = "False"
    fixed = {None: [], "fixatoms": [0, 1, 2], "fixatoms_dof": [0]}[selection]
    free = sorted(set(range(6)) - set(fixed))
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        system = sim.syslist[0]
        initial = system.beads.q.copy()
        sim.run(steps=150 if selection else 20, write_outputs=False)
        np.testing.assert_array_equal(system.beads.q[:, fixed], initial[:, fixed])
        np.testing.assert_allclose(
            system.beads.q[:, free], 0, atol=2e-4 if selection else 1e-5
        )
        if selection:
            assert system.forces.pot == pytest.approx(
                0.25 * np.sum(initial[:, fixed] ** 2), abs=1e-7
            )
        else:
            assert system.forces.pot < 1e-10


@pytest.mark.parametrize("threading", [False, True])
def test_properties_output_stride_matches_completed_steps(simulation, threading):
    """A real output file contains step zero and exactly the requested subsequent steps."""
    xml = harmonic_input(threading=threading).replace(
        "<system>",
        '<output prefix="sample"><properties filename="properties" stride="2">'
        "[step, potential, kinetic_md]</properties></output><system>",
    )
    with simulation(xml) as sim:
        sim.run(steps=4)
        system = sim.syslist[0]
        expected_potential = 0.25 * np.sum(system.beads.q**2)
        expected_kinetic = np.sum(system.beads.p**2 / (2 * system.beads.m3))
    rows = np.loadtxt("sample.properties")
    np.testing.assert_array_equal(rows[:, 0], [0, 2, 4])
    np.testing.assert_allclose(
        rows[-1, 1:], [expected_potential, expected_kinetic], rtol=1e-12
    )


@pytest.mark.parametrize("estimator", ["get_isotope_yama", "get_isotope_thermo"])
@pytest.mark.parametrize(
    "spring",
    [
        0.5,
        pytest.param(
            -0.5,
            marks=pytest.mark.xfail(
                strict=True,
                raises=AssertionError,
                reason="Signed log accumulation gives an incorrect negative isotope kinetic estimator",
            ),
        ),
    ],
)
def test_identical_isotope_weights_recover_mean_kinetic_estimator(
    simulation, estimator, spring
):
    """At mass ratio one, reweighting must recover the ordinary mean, including negative estimates."""
    from ipi.utils.units import unit_to_internal

    xml = harmonic_input(nbeads=4).replace("{k1: 0.5}", f"{{k1: {spring}}}")
    with simulation(xml) as sim:
        system = sim.syslist[0]
        system.beads.q[:, 0] = [-2, -1, 1, 2]
        system.beads.q[:, 3] = [-4, -2, 2, 4]
        original_positions = system.beads.q.copy()
        estimate = getattr(system.properties, estimator)(alpha="1")
        # The atoms' x variances are 2.5 and 10. The mean virial term is k*12.5/4.
        expected = (
            1.5 * unit_to_internal("temperature", "kelvin", 300) + spring * 12.5 / 4
        )
        assert estimate[0] == pytest.approx(0, abs=1e-12)
        assert np.exp(estimate[4]) == pytest.approx(2)
        assert estimate[2] == pytest.approx(expected)
        weighted_mean = estimate[6] * np.exp(estimate[5] - estimate[4])
        assert weighted_mean == pytest.approx(expected)
        np.testing.assert_array_equal(system.beads.q, original_positions)


@pytest.mark.parametrize("statistics", ["cl", "qn"])
@pytest.mark.parametrize(
    "nparallel",
    [
        1,
        pytest.param(
            2,
            marks=pytest.mark.xfail(
                strict=True,
                raises=ValueError,
                reason="SCP parallel sampling tries to clone one bead into two, which Beads.clone rejects",
            ),
        ),
    ],
)
def test_self_consistent_phonon_sampling_of_harmonic_potential(
    simulation, statistics, nparallel
):
    """Trial covariance and sampled forces agree with an independent harmonic solution."""
    import xml.etree.ElementTree as ET
    from ipi.utils.units import unit_to_internal

    root = ET.fromstring(harmonic_input())
    motion = root.find("system/motion")
    motion.clear()
    motion.set("mode", "scp")
    ET.SubElement(motion, "fixcom").text = "False"
    options = ET.SubElement(motion, "scp", mode=statistics)
    ET.SubElement(options, "max_steps").text = "4"
    ET.SubElement(options, "max_iter").text = "2"
    ET.SubElement(options, "nparallel").text = str(nparallel)
    masses = np.repeat([2.0, 3.0], 3)
    dynamical = np.diag(0.5 / masses)
    ET.SubElement(options, "dynmat", shape="(6,6)").text = str(
        dynamical.ravel().tolist()
    )
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        sim.run(steps=1 + 4 // nparallel, write_outputs=False)
        mover = sim.syslist[0].motion
        temperature = unit_to_internal("temperature", "kelvin", 300)
        omega = np.sqrt(0.5 / masses)
        if statistics == "cl":
            variance = np.full(6, temperature / 0.5)
        else:
            variance = 1 / (2 * masses * omega * np.tanh(omega / (2 * temperature)))
        np.testing.assert_allclose(mover.D, np.diag(variance), rtol=1e-12, atol=1e-14)
        np.testing.assert_allclose(mover.K, 0.5 * np.eye(6), atol=1e-14)
        samples = mover.phononator.x[0]
        np.testing.assert_allclose(
            samples.mean(axis=0), sim.syslist[0].beads.q[0], atol=1e-14
        )
        np.testing.assert_allclose(mover.phononator.f[0], -0.5 * samples, atol=1e-14)
        np.testing.assert_allclose(
            mover.phononator.v[0], 0.25 * np.sum(samples**2, axis=1), atol=1e-14
        )


@pytest.mark.xfail(
    strict=True,
    raises=TypeError,
    reason="InputMotion passes fixatoms_dof but ConstrainedDynamics still accepts fixatoms",
)
@pytest.mark.parametrize("ensemble, splitting", [("nve", "obabo")])
def test_constrained_dynamics_preserves_bond_and_tangent_velocities(
    simulation, ensemble, splitting
):
    """The actual integrator preserves the bond while allowing tangential motion."""
    import xml.etree.ElementTree as ET

    root = ET.fromstring(harmonic_input())
    motion = root.find("system/motion")
    motion.clear()
    motion.set("mode", "constrained_dynamics")
    ET.SubElement(motion, "fixcom").text = "False"
    options = ET.SubElement(
        motion, "constrained_dynamics", mode=ensemble, splitting=splitting
    )
    ET.SubElement(options, "timestep").text = "0.01"
    ET.SubElement(options, "nsteps_geo").text = "2"
    ET.SubElement(options, "nsteps_o").text = "2"
    solver = ET.SubElement(options, "csolver")
    ET.SubElement(solver, "tolerance").text = "1e-12"
    constraint = ET.SubElement(options, "constraint", mode="distance")
    ET.SubElement(constraint, "atoms").text = "[0,1]"
    if ensemble == "nvt":
        thermostat = ET.SubElement(options, "thermostat", mode="langevin")
        ET.SubElement(thermostat, "tau").text = "0.5"
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        beads = sim.syslist[0].beads
        initial = beads.q.copy()
        bond_length = np.linalg.norm(initial[0, :3] - initial[0, 3:])
        for _ in range(5):
            sim.run(steps=1, write_outputs=False)
            bond = beads.q[0, :3] - beads.q[0, 3:]
            velocity = beads.p[0] / beads.m3[0]
            assert np.linalg.norm(bond) == pytest.approx(bond_length, abs=1e-11)
            assert np.dot(bond, velocity[:3] - velocity[3:]) == pytest.approx(
                0, abs=1e-11
            )
        assert np.linalg.norm(beads.q - initial) > 1e-4


@pytest.mark.xfail(
    strict=True,
    raises=ValueError,
    reason="GLE initializer looks for the thermostat on Ensemble instead of Motion",
)
def test_gle_auxiliary_momenta_can_be_initialized_from_input(simulation):
    xml = harmonic_input(thermostat="gle").replace(
        "<tau>0.5</tau>", "<A shape='(2,2)'>[1,0,0,1]</A>"
    )
    xml = xml.replace(
        "<system>",
        "<system><initialize nbeads='1'><gle mode='manual'>[1,2,3,4,5,6,7,8,9,10,11,12]</gle></initialize>",
    )
    try:
        with simulation(xml) as sim:
            np.testing.assert_array_equal(
                sim.syslist[0].motion.thermostat.s.ravel(), np.arange(1, 13)
            )
    except ValueError as error:
        assert str(error) == "Ensemble does not have a thermostat to initialize"
        raise


@pytest.mark.parametrize("source", ["string", "file"])
def test_read_only_initialization_binds_state_without_starting_services(
    tmp_path, monkeypatch, source
):
    import io
    from ipi.engine.simulation import Simulation

    monkeypatch.chdir(tmp_path)
    xml = harmonic_input()
    observed = []

    def initialize(sim):
        observed.append(sim.chk)
        sim.syslist[0].beads.q[:] = 2

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
    from ipi.engine.simulation import Simulation

    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError, match=message):
        Simulation.load_from_xml(xml, read_only=True)
    assert not list(tmp_path.iterdir())


def test_multiple_systems_initialize_independent_state_and_outputs(simulation):
    import copy
    import xml.etree.ElementTree as ET

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
        first.beads.q[:] = 0
        np.testing.assert_array_equal(second.beads.q, 1)
        assert {out.filename for out in sim.outputs} == {
            "first_initial.energy",
            "second_initial.energy",
        }

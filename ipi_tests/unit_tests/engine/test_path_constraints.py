"""Scientific and input/restart checks for centroid and replica constraints."""

import io

import numpy as np
import pytest

from ipi.engine.beads import Beads
from ipi.engine.cell import Cell
from ipi.engine.forces import Forces
from ipi.engine.motion.planetary import Planetary
from ipi.engine.thermostats import ThermoPILE_L
from ipi.inputs.beads import InputBeads
from ipi.inputs.cell import InputCell
from ipi.inputs.motion.motion import InputMotion
from ipi.scripting import InteractiveSimulation
from ipi.utils.depend import dstrip
from ipi.utils.softexit import softexit
from ipi.utils.units import Constants


def constraint_xml(
    constraints="",
    mode="nve",
    splitting="obabo",
    nbeads=4,
    nmts=1,
    thermostat="langevin",
    normal_modes="",
    motion_mode="dynamics",
    forces="<force forcefield='harm'/>",
    bias="",
):
    """Two harmonic atoms, with spread replicas and nonzero initial momenta."""
    rng = np.random.default_rng(932)
    beads = Beads(natoms=2, nbeads=nbeads)
    beads.q = rng.normal(size=(nbeads, 6))
    beads.p = rng.normal(size=(nbeads, 6))
    beads.m = np.array([2.0, 2.0])
    beads.names = np.array(["H", "O"])
    ib = InputBeads()
    ib.store(beads)
    ic = InputCell()
    ic.store(Cell(20.0 * np.eye(3)))
    return f"""<simulation verbosity='quiet' floatformat='%24.16e'>
    <total_steps>100</total_steps>
    <prng><seed>9932</seed></prng>
    <ffdirect name='harm'><pes>harmonic</pes><parameters>{{k1: 0.2}}</parameters></ffdirect>
    <ffdirect name='bias'><pes>harmonic</pes><parameters>{{k1: 0.1}}</parameters></ffdirect>
    <system>
      {ib.write('beads')}
      {ic.write('cell')}
      <forces>{forces}</forces>
      <ensemble><temperature units='kelvin'>10000</temperature>{bias}</ensemble>
      <motion mode='{motion_mode}'>
        {constraints}
        <dynamics mode='{mode}' splitting='{splitting}'>
          <timestep>0.05</timestep>
          <nmts>[{nmts}]</nmts>
          <thermostat mode='{thermostat}'><tau>0.5</tau></thermostat>
        </dynamics>
      </motion>
      {normal_modes}
    </system></simulation>"""


@pytest.fixture
def simulations(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    running = []

    def make(**kwargs):
        sim = InteractiveSimulation(io.StringIO(constraint_xml(**kwargs)))
        running.append(sim)
        return sim

    yield make
    for sim in running:
        sim.stop()
    softexit.reset()


@pytest.mark.parametrize(
    "constraints", ["<fixcentroid>True</fixcentroid>", "<fixbeads>[0,2]</fixbeads>"]
)
@pytest.mark.parametrize(
    "mode,splitting", [("nve", "obabo"), ("nvt", "obabo"), ("nvt", "baoab")]
)
@pytest.mark.parametrize("nmts", [1, 2, 3])
def test_constraints_preserved_while_free_coordinates_move(
    simulations, constraints, mode, splitting, nmts
):
    sim = simulations(
        constraints=constraints, mode=mode, splitting=splitting, nmts=nmts
    )
    system = sim.syslist[0]
    initial_q = dstrip(system.beads.q).copy()
    initial_centroid = initial_q.mean(axis=0)
    for _ in range(25):
        system.motion.step()
        if system.motion.fixcentroid:
            np.testing.assert_allclose(
                system.beads.qc, initial_centroid, atol=2e-14, rtol=0
            )
            np.testing.assert_allclose(system.nm.pnm[0], 0.0, atol=2e-14, rtol=0)
        else:
            np.testing.assert_array_equal(system.beads.q[[0, 2]], initial_q[[0, 2]])
            np.testing.assert_allclose(system.beads.p[[0, 2]], 0.0, atol=2e-14)
    assert np.linalg.norm(dstrip(system.beads.q)[1] - initial_q[1]) > 0.01
    assert system.motion.fixcom is False
    assert system.motion.get_fixdof() == (6 if system.motion.fixcentroid else 12)


@pytest.mark.parametrize(
    "constraints,expected",
    [
        ("", 8),
        ("<fixcentroid>True</fixcentroid>", 12),
        ("<fixbeads>[0,2]</fixbeads>", 16),
    ],
)
@pytest.mark.parametrize(
    "mode,splitting", [("nve", "obabo"), ("nvt", "obabo"), ("nvt", "baoab")]
)
def test_fixed_atom_dofs_and_path_constraints_count_overlap_once(
    simulations, constraints, expected, mode, splitting
):
    sim = simulations(
        constraints=constraints + "<fixatoms_dof>[0,4]</fixatoms_dof>",
        mode=mode,
        splitting=splitting,
    )
    system = sim.syslist[0]
    # Identical positions across replicas give zero springs for fixed atoms.
    system.beads.q[:, [0, 4]] = system.beads.q[0, [0, 4]]
    q = dstrip(system.beads.q).copy()
    assert system.motion.get_fixdof() == expected
    assert system.motion.thermostat.ndof == 24 - expected
    for _ in range(20):
        system.motion.step()
    np.testing.assert_allclose(system.beads.q[:, [0, 4]], q[:, [0, 4]], atol=2e-14)
    np.testing.assert_allclose(system.beads.p[:, [0, 4]], 0, atol=2e-14)


def test_fixed_atoms_with_spread_replicas_keep_existing_spring_dynamics(simulations):
    system = simulations(constraints="<fixatoms_dof>[0]</fixatoms_dof>").syslist[0]
    q = dstrip(system.beads.q)[:, 0].copy()
    # With zero momenta, each Fourier mode evolves as q(t) = q(0) cos(omega t).
    frequencies = 2 * system.nm.omegan * np.sin(np.arange(3) * np.pi / 4)
    expected = np.fft.irfft(
        np.fft.rfft(q) * np.cos(frequencies * system.motion.integrator.qdt), n=4
    )
    system.motion.integrator.free_qstep()
    np.testing.assert_allclose(system.beads.q[:, 0], expected, atol=2e-14, rtol=0)
    assert np.max(np.abs(expected - q)) > 1e-6


@pytest.mark.parametrize(
    "constraints", ["<fixcentroid>True</fixcentroid>", "<fixbeads>[0]</fixbeads>"]
)
def test_nve_constrained_harmonic_energy_and_reversibility(simulations, constraints):
    sim = simulations(constraints=constraints)
    system = sim.syslist[0]
    q0 = dstrip(system.beads.q).copy()
    p0 = dstrip(system.beads.p).copy()
    energy0 = (
        system.nm.kin + system.nm.vspring + system.forces.pot + system.ensemble.eens
    )
    for _ in range(100):
        system.motion.step()
    energy1 = (
        system.nm.kin + system.nm.vspring + system.forces.pot + system.ensemble.eens
    )
    assert abs(energy1 - energy0) < 5e-4 * energy0
    system.beads.p *= -1
    for _ in range(100):
        system.motion.step()
    np.testing.assert_allclose(system.beads.q, q0, atol=2e-12, rtol=0)
    np.testing.assert_allclose(system.beads.p, -p0, atol=2e-12, rtol=0)


def test_centroid_constraint_uses_pile_local_and_propagates_internal_modes(simulations):
    sim = simulations(
        constraints="<fixcentroid>True</fixcentroid>", mode="nvt", thermostat="pile_l"
    )
    system = sim.syslist[0]
    q0 = dstrip(system.nm.qnm).copy()
    for _ in range(30):
        system.motion.step()
    np.testing.assert_allclose(system.nm.qnm[0], q0[0], atol=1e-14)
    assert np.linalg.norm(system.nm.qnm[1:] - q0[1:]) > 0.01


@pytest.mark.parametrize("thermostat", ["pile_g", "svr", "gle", "nm_gle", "nm_gle_g"])
def test_centroid_constraint_retains_previous_thermostat_options(
    tmp_path, monkeypatch, thermostat
):
    monkeypatch.chdir(tmp_path)
    xml = constraint_xml(
        constraints="<fixcentroid>True</fixcentroid>", mode="nvt", thermostat=thermostat
    )
    if thermostat == "gle":
        xml = xml.replace("<tau>0.5</tau>", "<A shape='(1,1)'>[2]</A>")
    elif thermostat.startswith("nm_gle"):
        xml = xml.replace(
            "<tau>0.5</tau>", "<tau>0.5</tau><A shape='(4,1,1)'>[2,2,2,2]</A>"
        )
    sim = InteractiveSimulation(io.StringIO(xml))
    try:
        system = sim.syslist[0]
        q0 = dstrip(system.nm.qnm).copy()
        for _ in range(30):
            system.motion.step()
        np.testing.assert_allclose(system.nm.qnm[0], q0[0], atol=2e-14, rtol=0)
        np.testing.assert_allclose(system.nm.pnm[0], 0, atol=2e-14)
        assert np.linalg.norm(system.nm.qnm[1:] - q0[1:]) > 0.01
        assert np.isfinite(system.ensemble.eens)
    finally:
        sim.stop()
        softexit.reset()


@pytest.mark.parametrize(
    "constraints", ["<fixcentroid>True</fixcentroid>", "<fixbeads>[0,2,0]</fixbeads>"]
)
def test_motion_input_roundtrip_preserves_options(simulations, constraints):
    original = simulations(constraints=constraints).syslist[0].motion
    stored = InputMotion()
    stored.store(original)
    restored = stored.fetch()
    assert restored.enstype == "nve"
    assert restored.fixcentroid == original.fixcentroid
    np.testing.assert_array_equal(restored.fixbeads, original.fixbeads)
    assert restored.fixcom is False
    if not restored.fixcentroid:
        np.testing.assert_array_equal(restored.fixbeads, [0, 2])


@pytest.mark.parametrize(
    "constraints,thermostat",
    [
        ("<fixcentroid>True</fixcentroid>", "langevin"),
        ("<fixcentroid>True</fixcentroid>", "pile_g"),
        ("<fixbeads>[0,2]</fixbeads>", "langevin"),
        ("<fixbeads>[0,2]</fixbeads>", "pile_l"),
    ],
)
def test_restart_matches_continuous_nvt(simulations, tmp_path, constraints, thermostat):
    first = simulations(
        constraints=constraints, mode="nvt", splitting="baoab", thermostat=thermostat
    )
    first.run(steps=7, write_outputs=False)
    restart = (tmp_path / "RESTART").read_text()
    first.run(steps=8, write_outputs=False)
    expected_q = dstrip(first.syslist[0].beads.q).copy()
    expected_p = dstrip(first.syslist[0].beads.p).copy()
    expected_work = first.syslist[0].ensemble.eens
    first.stop()
    softexit.reset()
    second = InteractiveSimulation(io.StringIO(restart))
    try:
        second.run(steps=8, write_outputs=False)
        np.testing.assert_allclose(
            second.syslist[0].beads.q, expected_q, atol=2e-13, rtol=0
        )
        np.testing.assert_allclose(
            second.syslist[0].beads.p, expected_p, atol=2e-13, rtol=0
        )
        assert second.syslist[0].ensemble.eens == pytest.approx(
            expected_work, abs=2e-12
        )
    finally:
        second.stop()
        softexit.reset()


@pytest.mark.parametrize(
    "constraints", ["<fixbeads>[0]</fixbeads>", "<fixcentroid>True</fixcentroid>"]
)
def test_one_replica_is_fully_frozen(simulations, constraints):
    system = simulations(constraints=constraints, nbeads=1).syslist[0]
    q0 = dstrip(system.beads.q).copy()
    for _ in range(5):
        system.motion.step()
    np.testing.assert_array_equal(system.beads.q, q0)
    np.testing.assert_allclose(system.beads.p, 0, atol=1e-14)
    assert system.motion.get_fixdof() == 6


def test_one_replica_with_fixed_centroid_and_pile_g_is_fully_frozen(simulations):
    system = simulations(
        constraints="<fixcentroid>True</fixcentroid>",
        mode="nvt",
        thermostat="pile_g",
        nbeads=1,
    ).syslist[0]
    q0 = dstrip(system.beads.q).copy()
    for _ in range(5):
        system.motion.step()
    np.testing.assert_array_equal(system.beads.q, q0)
    np.testing.assert_array_equal(system.beads.p, 0)
    assert np.isfinite(system.ensemble.eens)


@pytest.mark.parametrize(
    "kwargs,match",
    [
        (
            {"constraints": "<fixcentroid>True</fixcentroid><fixbeads>[0]</fixbeads>"},
            "cannot be used together",
        ),
        ({"constraints": "<fixbeads>[-1]</fixbeads>"}, "out of bounds"),
        ({"constraints": "<fixbeads>[4]</fixbeads>"}, "out of bounds"),
        (
            {"constraints": "<fixbeads>[0]</fixbeads>", "mode": "npt"},
            "not implemented for NPT/NST",
        ),
        (
            {"constraints": "<fixcentroid>True</fixcentroid>", "mode": "sc"},
            "require NVE or NVT",
        ),
        (
            {"constraints": "<fixcentroid>True</fixcentroid>", "motion_mode": "dummy"},
            "motion mode='dynamics'",
        ),
        ({"mode": "nvt-cc"}, "nvt-cc"),
        (
            {
                "constraints": "<fixbeads>[0]</fixbeads>",
                "mode": "nvt",
                "thermostat": "pile_g",
            },
            "Cartesian Langevin or PILE-L",
        ),
        (
            {
                "constraints": "<fixbeads>[0]</fixbeads>",
                "normal_modes": "<normal_modes><open_paths>[0]</open_paths></normal_modes>",
            },
            "open paths or bosonic",
        ),
        (
            {
                "constraints": "<fixbeads>[0]</fixbeads>",
                "normal_modes": "<normal_modes propagator='cayley'/>",
            },
            "exact spring propagator",
        ),
        (
            {
                "constraints": "<fixbeads>[0]</fixbeads>",
                "normal_modes": "<normal_modes><frequencies style='manual'>[1,1,1]</frequencies></normal_modes>",
            },
            "RPMD dynamical masses",
        ),
        (
            {
                "constraints": "<fixbeads>[0]</fixbeads>",
                "normal_modes": "<normal_modes><frequencies style='eco'>[1]</frequencies></normal_modes>",
            },
            "RPMD dynamical masses",
        ),
        (
            {
                "constraints": "<fixcentroid>True</fixcentroid>",
                "normal_modes": "<normal_modes propagator='bab'><bosons id='index'>[0]</bosons></normal_modes>",
            },
            "bosonic exchange",
        ),
    ],
)
def test_invalid_constraint_combinations_fail_before_dynamics(
    simulations, kwargs, match
):
    with pytest.raises((ValueError, NotImplementedError), match=match):
        simulations(**kwargs)


def test_planetary_uses_standard_nvt_with_centroid_option():
    planetary = Planetary(timestep=0.1, thermostat=ThermoPILE_L())
    assert planetary.ccdyn.enstype == "nvt"
    assert planetary.ccdyn.fixcentroid
    assert not planetary.ccdyn.fixcom


@pytest.mark.parametrize(
    "constraints", ["<fixcentroid>True</fixcentroid>", "<fixbeads>[0,2]</fixbeads>"]
)
@pytest.mark.parametrize(
    "selection", [{}, {"bead": "0"}, {"bead": "1"}, {"nm": "0"}, {"nm": "1"}]
)
def test_temperature_restores_only_missing_kinetic_directions(
    simulations, constraints, selection
):
    system = simulations(
        constraints=constraints + "<fixatoms_dof>[0]</fixatoms_dof>"
    ).syslist[0]
    system.beads.p[:] = 0
    # In this zero-momentum state the reported temperature is exactly the
    # constrained fraction of the reference temperature (i-PI convention).
    if system.motion.fixcentroid:
        fraction = (
            (1.0 if selection["nm"] == "0" else 0.0) if "nm" in selection else 0.25
        )
    elif "bead" in selection:
        fraction = 1.0 if selection["bead"] == "0" else 0.0
    elif selection.get("nm") == "1":
        # Real cosine mode has support only at beads 0 and 2 for P=4.
        fraction = 1.0
    else:
        fraction = 0.5
    for atom, fixed_atoms, ndof in [("", 1, 6), ("0", 1, 3), ("O", 0, 3)]:
        expected = (
            system.ensemble.temp
            * (fixed_atoms + (ndof - fixed_atoms) * fraction)
            / ndof
        )
        assert system.properties.get_temp(atom=atom, **selection) == pytest.approx(
            expected
        )


@pytest.mark.parametrize(
    "centroid,thermostat",
    [(False, "langevin"), (False, "pile_l"), (True, "langevin"), (True, "pile_g")],
)
def test_harmonic_conditional_statistics(simulations, centroid, thermostat):
    """Sample mean and variance match the constrained harmonic density."""
    constraints = (
        "<fixcentroid>True</fixcentroid>" if centroid else "<fixbeads>[0]</fixbeads>"
    )
    system = simulations(
        constraints=constraints, mode="nvt", splitting="baoab", thermostat=thermostat
    ).syslist[0]
    system.motion.dt = 0.2
    system.motion.thermostat.tau = 3.0
    omega2 = system.nm.omegan**2
    thermal_energy = 4 * Constants.kb * system.ensemble.temp
    initial_q = dstrip(system.beads.q).copy()
    if centroid:
        # The first non-centroid real Fourier mode for P=4 has lambda=2.
        expected_mean = np.zeros(6)
        expected_var = thermal_energy / (0.2 + 2 * 2.0 * omega2)
    else:
        # Three free coordinates with two spring links to the fixed replica.
        stiffness = np.diag(np.full(3, 0.2 + 4 * omega2))
        stiffness += np.diag(np.full(2, -2 * omega2), 1)
        stiffness += np.diag(np.full(2, -2 * omega2), -1)
        expected_mean = (
            np.linalg.solve(stiffness, 2 * omega2 * np.array([1, 0, 1]))[:, None]
            * initial_q[0]
        )
        expected_var = thermal_energy * np.diag(np.linalg.inv(stiffness))
    samples = []
    temperatures = []
    for step in range(14000):
        system.motion.step()
        if step >= 2000 and step % 5 == 0:
            samples.append(
                dstrip(system.nm.qnm)[1].copy()
                if centroid
                else dstrip(system.beads.q)[1:].copy()
            )
            temperatures.append(system.properties.get_temp())
    samples = np.asarray(samples)
    mean_error = samples.mean(axis=0) - expected_mean
    assert np.max(np.abs(mean_error)) < 0.18 * np.sqrt(np.max(expected_var))
    variance = (
        samples.var(axis=0).mean() if centroid else samples.var(axis=0).mean(axis=1)
    )
    np.testing.assert_allclose(variance, expected_var, rtol=0.16)
    assert np.mean(temperatures) == pytest.approx(system.ensemble.temp, rel=0.07)


@pytest.mark.parametrize("centroid", [False, True])
@pytest.mark.parametrize("multiple_levels", [False, True])
def test_harmonic_trajectory_with_bias_and_mts_matches_exact_solution(
    simulations, centroid, multiple_levels
):
    constraints = (
        "<fixcentroid>True</fixcentroid>" if centroid else "<fixbeads>[0]</fixbeads>"
    )
    kwargs = {}
    if multiple_levels:
        kwargs.update(
            nmts="1,2",
            forces="""<force forcefield='harm' weight='0.5'><mts_weights>[1,0]</mts_weights></force>
          <force forcefield='harm' weight='0.5'><mts_weights>[0,1]</mts_weights></force>""",
        )
    system = simulations(
        constraints=constraints,
        nbeads=2,
        bias="<bias><force forcefield='bias'/></bias>",
        **kwargs,
    ).syslist[0]
    mass = 2.0
    if centroid:
        q0 = dstrip(system.nm.qnm)[1].copy()
        p0 = dstrip(system.nm.pnm)[1].copy()
        equilibrium = 0.0
        spring_stiffness = 4 * mass * system.nm.omegan**2
    else:
        q0 = dstrip(system.beads.q)[1].copy()
        p0 = dstrip(system.beads.p)[1].copy()
        spring_stiffness = 2 * mass * system.nm.omegan**2
        equilibrium = (
            spring_stiffness * dstrip(system.beads.q)[0] / (0.3 + spring_stiffness)
        )
    frequency = np.sqrt((0.3 + spring_stiffness) / mass)
    for _ in range(80):
        system.motion.step()
    time = 80 * system.motion.dt
    expected_q = (
        equilibrium
        + (q0 - equilibrium) * np.cos(frequency * time)
        + p0 * np.sin(frequency * time) / (mass * frequency)
    )
    expected_p = p0 * np.cos(frequency * time) - mass * frequency * (
        q0 - equilibrium
    ) * np.sin(frequency * time)
    actual_q = system.nm.qnm[1] if centroid else system.beads.q[1]
    actual_p = system.nm.pnm[1] if centroid else system.beads.p[1]
    # Second-order splitting error at dt=0.05 atomic units.
    np.testing.assert_allclose(actual_q, expected_q, atol=3e-4, rtol=0)
    np.testing.assert_allclose(actual_p, expected_p, atol=3e-4, rtol=0)


def test_constraints_on_multi_parent_are_rejected(simulations):
    with pytest.raises(ValueError, match="each dynamics motion"):
        simulations(constraints="<fixcentroid>True</fixcentroid>", motion_mode="multi")


def test_unconstrained_default_motion_still_fixes_com(simulations):
    system = simulations().syslist[0]
    assert system.motion.fixcom
    assert not system.motion.fixcentroid
    assert not len(system.motion.fixbeads)
    assert system.motion.get_fixdof() == 3
    assert system.motion.thermostat.ndof == 21


def test_planetary_sampling_tracks_updated_parent_centroid(tmp_path, monkeypatch):
    import re

    monkeypatch.chdir(tmp_path)
    # Isolate the pre-existing planetary call to the obsolete Forces.copy API.
    monkeypatch.setattr(Forces, "copy", Forces.clone, raising=False)
    xml = constraint_xml()
    xml = re.sub(
        r"<motion.*?</motion>",
        """<motion mode='planetary'><planetary>
      <timestep>0.05</timestep><nsamples>3</nsamples>
      <thermostat mode='pile_l'><tau>0.5</tau></thermostat>
      </planetary></motion>""",
        xml,
        flags=re.S,
    )
    sim = InteractiveSimulation(io.StringIO(xml))
    try:
        system = sim.syslist[0]
        # Exercise the real sampler; isolate its existing text/binary matrix
        # serializer, which is independent of the constraint implementation.
        matrices = []
        monkeypatch.setattr(system.motion, "save_matrix", matrices.append)
        for shift in (0.0, 0.5):
            system.beads.q += shift
            expected = dstrip(system.beads.qc).copy()
            system.motion.step(0)
            np.testing.assert_allclose(system.motion.dbeads.qc, expected, atol=1e-14)
            np.testing.assert_allclose(system.motion.dnm.pnm[0], 0.0, atol=1e-14)
        assert len(matrices) == 2
        assert all(np.isfinite(matrix.a).all() for matrix in matrices)
        stored = InputMotion()
        stored.store(system.motion)
        assert not stored.fixcentroid.fetch()
        assert stored.fetch().ccdyn.fixcentroid
    finally:
        sim.stop()
        softexit.reset()


@pytest.mark.parametrize("mode", ["npt", "nst", "scnpt"])
@pytest.mark.parametrize(
    "constraint", ["<fixcentroid>True</fixcentroid>", "<fixbeads>[0]</fixbeads>"]
)
def test_constant_pressure_constraints_raise_explicit_not_implemented(
    simulations, mode, constraint
):
    """Unsupported barostat coupling is rejected during dynamics binding."""
    with pytest.raises(NotImplementedError, match="not implemented for NPT/NST"):
        simulations(mode=mode, constraints=constraint)


@pytest.mark.parametrize(
    "first,second",
    [
        ("<fixbeads>[0]</fixbeads>", "<fixbeads>[0,2]</fixbeads>"),
        ("<fixbeads>[0,2]</fixbeads>", "<fixbeads>[0,3]</fixbeads>"),
        ("<fixbeads>[0]</fixbeads>", ""),
        ("", "<fixbeads>[0]</fixbeads>"),
    ],
)
def test_different_fixed_bead_sets_in_combined_motions_are_not_implemented(
    simulations, first, second
):
    import re

    xml = constraint_xml(constraints=first)
    other = constraint_xml(constraints=second)
    child = re.search(r"<motion.*?</motion>", xml, re.S).group()
    other = re.search(r"<motion.*?</motion>", other, re.S).group()
    xml = xml.replace(child, f"<motion mode='multi'>{child}{other}</motion>")
    with pytest.raises(NotImplementedError, match="Different fixed-bead sets"):
        InteractiveSimulation(io.StringIO(xml))


def test_identical_fixed_beads_in_combined_motions_match_separate_dynamics(simulations):
    """Children with identical constraints share one cached spring propagator."""
    import re

    xml = constraint_xml(constraints="<fixbeads>[0]</fixbeads>")
    first = re.search(r"<motion.*?</motion>", xml, re.S).group()
    xml = xml.replace(first, f"<motion mode='multi'>{first}{first}</motion>")
    sim = InteractiveSimulation(io.StringIO(xml))
    try:
        system = sim.syslist[0]
        references = [
            simulations(constraints="<fixbeads>[0]</fixbeads>").syslist[0]
            for _ in range(2)
        ]
        for _ in range(10):
            for child, reference in zip(system.motion.mlist, references):
                reference.beads.q[:] = system.beads.q
                reference.beads.p[:] = system.beads.p
                child.step()
                reference.motion.step()
                np.testing.assert_allclose(
                    system.beads.q, reference.beads.q, atol=2e-14
                )
                np.testing.assert_allclose(
                    system.beads.p, reference.beads.p, atol=2e-14
                )
    finally:
        sim.stop()


@pytest.mark.parametrize(
    "first,second,common",
    [
        (
            "<fixcentroid>True</fixcentroid>",
            "<fixcentroid>True</fixcentroid>",
            "<fixcentroid>True</fixcentroid>",
        ),
        (
            "<fixbeads>[0,2]</fixbeads>",
            "<fixbeads>[0,2]</fixbeads>",
            "<fixbeads>[0,2]</fixbeads>",
        ),
        ("<fixcentroid>True</fixcentroid>", "", ""),
    ],
)
def test_combined_motion_temperature_accounts_for_common_constraints(
    simulations, first, second, common
):
    import re

    xml = constraint_xml(constraints=first, mode="nvt")
    child = re.search(r"<motion.*?</motion>", xml, re.S).group()
    xml = xml.replace(
        child, f"<motion mode='multi'>{child}{child.replace(first, second, 1)}</motion>"
    )
    sim = InteractiveSimulation(io.StringIO(xml))
    try:
        system = sim.syslist[0]
        reference = simulations(
            constraints=common + "<fixcom>False</fixcom>", mode="nvt"
        ).syslist[0]
        assert system.motion.fixcentroid == reference.motion.fixcentroid
        np.testing.assert_array_equal(system.motion.fixbeads, reference.motion.fixbeads)
        for selection in ({}, {"bead": "0"}, {"bead": "1"}, {"nm": "0"}, {"nm": "1"}):
            system.beads.p[:] = reference.beads.p[:] = 0
            assert system.properties.get_temp(**selection) == pytest.approx(
                reference.properties.get_temp(**selection)
            )
        stored = InputMotion()
        stored.store(system.motion)
        restored = stored.fetch()
        assert restored.fixcentroid == system.motion.fixcentroid
        np.testing.assert_array_equal(restored.fixbeads, system.motion.fixbeads)
    finally:
        sim.stop()

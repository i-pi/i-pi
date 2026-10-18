"""Isotope and kinetic estimators in independently solvable limits."""

import numpy as np
import pytest
from ipi.engine.properties import Trajectories
from ipi_tests.unit_tests.engine.test_simulation_invariants import harmonic_input


@pytest.mark.parametrize("method", ["get_isotope_zetatd", "get_isotope_zetatd_4th"])
@pytest.mark.parametrize(
    "alpha", [0.5, 1.0, 2.0], ids=["half-mass", "same-mass", "double-mass"]
)
@pytest.mark.parametrize(
    "atom,indices",
    [("", [0, 1]), ("H", [0]), ("1", [1])],
    ids=["all-atoms", "hydrogen", "oxygen-index"],
)
def test_free_particle_isotope_weight_equals_spring_boltzmann_ratio(
    simulation, method, alpha, atom, indices
):
    xml = harmonic_input(nbeads=4).replace("{k1: 0.5}", "{k1: 0.0}")
    with simulation(xml) as sim:
        s = sim.syslist[0]
        q = s.beads.q.copy().reshape(4, 2, 3)
        # The discretized free-particle action is a quadratic form in adjacent displacements.
        differences = q - np.roll(q, 1, axis=0)
        temp = s.ensemble.temp
        action = (
            0.5
            * np.asarray([2.0, 3.0])
            * (4 * temp) ** 2
            * np.sum(differences**2, axis=(0, 2))
        )
        delta = (alpha - 1) * action[indices]
        weight = alpha**4.5 * np.mean(np.exp(-delta / (4 * temp)))
        expected = [delta.mean(), np.mean(delta**2), weight]
        if method.endswith("4th"):
            expected += [weight, weight]
        np.testing.assert_allclose(
            getattr(s.properties, method)(str(alpha), atom),
            expected,
            rtol=1e-12,
            atol=1e-15,
        )
        trajectory = Trajectories()
        trajectory.bind(s)
        per_atom_delta = (alpha - 1) * action
        per_atom_delta[[i for i in range(2) if i not in indices]] = 0
        per_atom_expected = np.column_stack(
            [per_atom_delta, per_atom_delta**2, np.exp(-per_atom_delta / (4 * temp))]
        )
        np.testing.assert_allclose(
            trajectory.get_isotope_zetatd(str(alpha), atom).reshape(2, 3),
            per_atom_expected,
            atol=1e-14,
        )


@pytest.mark.parametrize("method", ["get_isotope_zetasc", "get_isotope_zetasc_4th"])
@pytest.mark.parametrize("atom", ["", "H", "1"])
def test_unchanged_isotope_has_unit_weight(simulation, method, atom):
    with simulation(harmonic_input(nbeads=4)) as sim:
        s = sim.syslist[0]
        original = s.beads.q.copy()
        actual = getattr(s.properties, method)("1", atom)
        np.testing.assert_allclose(actual, [0, 0] + [1] * (len(actual) - 2), atol=1e-13)
        np.testing.assert_array_equal(s.beads.q, original)


@pytest.mark.parametrize("alpha", [0.5, 2.0])
def test_scaled_isotope_harmonic_weight_from_variance(simulation, alpha):
    with simulation(harmonic_input(nbeads=4)) as sim:
        s = sim.syslist[0]
        q = s.beads.q.copy()
        fluctuations = q - q.mean(axis=0)
        delta = 0.25 * (1 / alpha - 1) * np.sum(fluctuations[:, :3] ** 2)
        expected = [delta, delta**2, np.exp(-delta / (4 * s.ensemble.temp))]
        np.testing.assert_allclose(
            s.properties.get_isotope_zetasc(str(alpha), "H"), expected, atol=1e-13
        )
        np.testing.assert_array_equal(s.beads.q, q)


@pytest.mark.parametrize(
    "method",
    [
        "get_isotope_zetatd",
        "get_isotope_zetatd_4th",
        "get_isotope_zetasc",
        "get_isotope_zetasc_4th",
        "get_rg",
        "get_ktens",
    ],
)
@pytest.mark.parametrize("atom", ["7", "missing"])
def test_atom_selectors_reject_absent_atoms(simulation, method, atom):
    with simulation(harmonic_input(nbeads=4)) as sim:
        prop = getattr(sim.syslist[0].properties, method)
        if method == "get_ktens" and atom == "missing":
            np.testing.assert_array_equal(prop(atom=atom), np.zeros(6))
        else:
            with pytest.raises(IndexError):
                prop(atom=atom)


@pytest.mark.parametrize("atom,indices", [("", [0, 1]), ("H", [0]), ("1", [1])])
def test_gyration_and_kinetic_tensor_of_harmonic_path(simulation, atom, indices):
    with simulation(harmonic_input(nbeads=4)) as sim:
        s = sim.syslist[0]
        q = s.beads.q.copy().reshape(4, 2, 3)
        centered = q - q.mean(axis=0)
        radii = np.sqrt(np.sum(centered**2, axis=(0, 2)) / 4)
        assert s.properties.get_rg(atom) == pytest.approx(radii[indices].mean())
        tensor = sum(
            0.5 * s.ensemble.temp * np.eye(3)
            + 0.25 * centered[:, i].T @ centered[:, i] / 4
            for i in indices
        )
        np.testing.assert_allclose(
            s.properties.get_ktens(atom),
            tensor[[0, 1, 2, 0, 0, 1], [0, 1, 2, 1, 2, 2]],
            atol=1e-13,
        )


@pytest.mark.parametrize(
    "method", ["get_yama_estimators", "get_scyama_estimators", "get_scop_estimators"]
)
def test_free_particle_energy_and_temperature_derivative(simulation, method):
    with simulation(harmonic_input(nbeads=4).replace("{k1: 0.5}", "{k1: 0.0}")) as sim:
        s = sim.syslist[0]
        t = s.ensemble.temp
        value = getattr(s.properties, method)(fd_delta=0.001)
        sign = 1 if method == "get_scop_estimators" else -1
        np.testing.assert_allclose(value, [3 * t, sign * 3 * t * t], atol=1e-13)


@pytest.mark.parametrize("atom,mass", [("H", 2), ("O", 3)])
def test_free_particle_end_to_end_distribution(simulation, atom, mass):
    with simulation(harmonic_input(nbeads=4).replace("{k1: 0.5}", "{k1: 0.0}")) as sim:
        s = sim.syslist[0]
        q = s.beads.q.copy()
        actual = s.properties.get_linlin(".2", "-.3", ".1", atom)
        assert actual == pytest.approx(
            np.exp(-float(mass) * (0.04 + 0.09 + 0.01) * s.ensemble.temp / 2)
        )
        np.testing.assert_array_equal(s.beads.q, q)
        with pytest.raises(IndexError):
            s.properties.get_linlin(atom="absent")


def test_harmonic_takahashi_imada_correction_has_force_squared_scaling(simulation):
    with simulation(harmonic_input(nbeads=4)) as sim:
        s = sim.syslist[0]
        q = s.beads.q.copy()
        temp = s.ensemble.temp
        expected = 0.25 * np.sum(q * q / s.beads.m3) / (24 * (4 * temp) ** 2 * 4)
        assert s.properties.get_ti_term() == pytest.approx(expected)
        correction = s.properties.get_ti_correction()
        exponent = -expected / temp
        np.testing.assert_allclose(
            correction, [exponent, exponent**2, np.exp(exponent)]
        )
        assert s.properties.get_ti_term("H") + s.properties.get_ti_term(
            "O"
        ) == pytest.approx(expected)
        assert s.properties.get_ti_term("absent") == 0
        with pytest.raises(IndexError):
            s.properties.get_ti_term("7")


@pytest.fixture
def property_system():
    """Bind a real harmonic system without output files or background services."""
    from ipi.engine.simulation import Simulation

    return Simulation.load_from_xml(harmonic_input(nbeads=4), read_only=True).syslist[0]


@pytest.mark.parametrize("atom,indices", [("", [0, 1]), ("H", [0]), ("1", [1])])
@pytest.mark.parametrize("mode", ["", "0", "1"])
def test_kinetic_temperature_and_normal_mode_selection(
    property_system, atom, indices, mode
):
    s = property_system
    p = np.array(
        [
            [1, 2, 3, 4, 5, 6],
            [2, 3, 4, 5, 6, 7],
            [3, 4, 5, 6, 7, 8],
            [4, 5, 6, 7, 8, 9],
        ],
        dtype=float,
    )
    s.beads.p[:] = p
    masses = np.array([2.0, 3.0])
    # The centroid mode is the normalized sum; the first real mode for four beads is (p0-p2)/sqrt(2).
    if mode == "0":
        momentum = p.sum(axis=0).reshape(2, 3) / 2
        kinetic = np.sum(momentum[indices] ** 2 / (2 * masses[indices, None]))
    elif mode == "1":
        momentum = ((p[0] - p[2]) / np.sqrt(2)).reshape(2, 3)
        kinetic = np.sum(momentum[indices] ** 2 / (2 * masses[indices, None]))
    else:
        kinetic = (
            np.sum(p.reshape(4, 2, 3)[:, indices] ** 2 / (2 * masses[indices, None]))
            / 4
        )
    # Orthogonal normal-mode transformation preserves the total momentum covariance.
    vectors = p.reshape(4, 2, 3)
    expected_stress = np.einsum("bai,baj,a->ij", vectors, vectors, 1 / masses)
    np.testing.assert_allclose(s.properties.get_kstress_md(), expected_stress)
    np.testing.assert_array_equal(s.properties.get_atom_ids(atom), indices)
    assert s.properties.get_kinmd(atom, nm=mode) == pytest.approx(kinetic)
    assert s.properties.get_temp(atom, nm=mode) == pytest.approx(
        2 * kinetic / (3 * len(indices) * 4)
    )


@pytest.mark.parametrize("atom", ["0", "1"])
@pytest.mark.parametrize("bead", [-1, 2])
def test_atom_vectors_and_bead_potential_match_coordinates(property_system, atom, bead):
    s = property_system
    q = s.beads.q.copy()
    expected = q[:, 3 * int(atom) : 3 * int(atom) + 3]
    np.testing.assert_array_equal(
        s.properties.get_atom_vec(q, atom, str(bead)),
        expected.mean(axis=0) if bead < 0 else expected[bead],
    )
    energy = 0.25 * np.sum(q * q, axis=1)
    assert s.properties.get_pot(bead=str(bead)) == pytest.approx(
        energy.mean() if bead < 0 else energy[bead]
    )


@pytest.mark.parametrize("reference", ["zero", "file", "wrong-size"])
def test_force_position_virial_uses_requested_reference(
    property_system, tmp_path, reference
):
    s = property_system
    q = s.beads.q.copy()
    path = tmp_path / "reference.dat"
    ref = np.arange(6) / 10 if reference != "zero" else np.zeros(6)
    np.savetxt(path, ref[:3] if reference == "wrong-size" else ref)
    if reference == "wrong-size":
        with pytest.raises(ValueError, match="reference"):
            s.properties.get_fqvirial(str(path))
    else:
        expected = -0.25 * np.sum(q * (q - ref)) / 4
        assert s.properties.get_fqvirial(
            "" if reference == "zero" else str(path)
        ) == pytest.approx(expected)


@pytest.mark.parametrize("suzuki_chin", [False, True])
def test_centroid_virial_kinetic_stress_matches_harmonic_covariance(
    property_system, suzuki_chin
):
    s = property_system
    q = s.beads.q.copy().reshape(4, 2, 3)
    centered = q - q.mean(axis=0)
    momenta = s.beads.p.copy().reshape(4, 2, 3).mean(axis=0)
    weights = np.full(4, 0.5)
    if suzuki_chin:
        # For alpha=0 and k=1/2, F_SC adds +/-F/3 and -2*k^2*q/(9*m*omega_P^2) on odd beads.
        weights = np.array([1 / 3, 2 / 3, 1 / 3, 2 / 3])
    expected = np.zeros((3, 3))
    for atom, mass in enumerate([2.0, 3.0]):
        stiffness = weights.copy()
        if suzuki_chin:
            stiffness[1::2] += 0.5 / (9 * mass * s.nm.omegan2)
        expected += np.einsum("b,bi,bj->ij", stiffness, centered[:, atom], q[:, atom])
        expected += np.diag(4 * momenta[atom] ** 2 / mass)
    actual = s.properties.kstress_sctd() if suzuki_chin else s.properties.kstress_cv()
    np.testing.assert_allclose(actual, np.triu(expected), rtol=1e-10, atol=1e-12)
    np.testing.assert_array_equal(
        s.properties.tensor2vec(actual), actual[[0, 1, 2, 0, 0, 1], [0, 1, 2, 1, 2, 2]]
    )


@pytest.mark.parametrize("alpha", [0.5, 1.0, 2.0])
@pytest.mark.parametrize("atom,indices", [("", [0, 1]), ("H", [0]), ("1", [1])])
def test_per_atom_scaled_isotope_weights_match_harmonic_variance(
    property_system, alpha, atom, indices
):
    s = property_system
    q = s.beads.q.copy()
    centered = (q - q.mean(axis=0)).reshape(4, 2, 3)
    delta = 0.25 * (1 / alpha - 1) * np.sum(centered**2, axis=(0, 2)) / 4
    delta[[i for i in range(2) if i not in indices]] = 0
    expected = np.column_stack([delta, delta**2, np.exp(-delta / s.ensemble.temp)])
    np.testing.assert_allclose(
        s.trajs.get_isotope_zetasc(str(alpha), atom).reshape(2, 3), expected, atol=1e-13
    )
    np.testing.assert_array_equal(s.beads.q, q)


@pytest.mark.parametrize(
    "method",
    [
        "get_kincv",
        "get_sckinop",
        "get_sckintd",
        "get_kintd",
        "get_isotope_yama",
        "get_isotope_thermo",
        "get_linlin",
    ],
)
@pytest.mark.parametrize("atom", ["2", "missing"])
def test_additional_estimator_atom_selection(property_system, method, atom):
    prop = getattr(property_system.properties, method)
    if atom == "missing" and method in [
        "get_kincv",
        "get_sckinop",
        "get_sckintd",
        "get_kintd",
    ]:
        assert prop(atom=atom) == 0
    else:
        with pytest.raises(IndexError):
            prop(atom=atom)


@pytest.mark.parametrize(
    "method",
    ["get_exchange_distinct_prob", "get_exchange_longest_prob", "get_fermionic_sign"],
)
def test_exchange_properties_require_bosons(property_system, method):
    with pytest.raises(Exception, match="No bosons"):
        getattr(property_system.properties, method)()


@pytest.mark.parametrize("bosons", [[0], [0, 1]])
def test_exchange_estimators_at_coincident_paths(bosons):
    from ipi.engine.simulation import Simulation

    xml = (
        harmonic_input(nbeads=4)
        .replace("[2, 3]", "[2, 2]")
        .replace(
            "<system>",
            f"<system><normal_modes><bosons id='index'>{bosons}</bosons></normal_modes>",
        )
    )
    s = Simulation.load_from_xml(xml, read_only=True).syslist[0]
    s.beads.q[:] = 0
    # With all spring actions zero, the two permutations have equal weight.
    assert s.properties.get_exchange_distinct_prob() == pytest.approx(
        1 if len(bosons) == 1 else 0.5
    )
    assert s.properties.get_exchange_longest_prob() == pytest.approx(
        1
    )  # normalized by 1/N
    assert s.properties.get_fermionic_sign() == pytest.approx(
        1 if len(bosons) == 1 else 0, abs=1e-14
    )
    assert s.properties.get_kintd() == pytest.approx(12 * s.ensemble.temp)
    with pytest.raises(IndexError, match="not applicable to bosons"):
        s.properties.get_kincv()
    if len(bosons) == 2:
        with pytest.raises(IndexError, match="subset"):
            s.properties.get_kintd("0")


def test_open_path_primitive_estimator_includes_ring_closure():
    # Characterization only: this reproduces the documented cyclic estimator.
    # It does not establish that averaging it over open paths gives physical KE.
    from ipi.engine.simulation import Simulation

    xml = harmonic_input(nbeads=4).replace(
        "<system>", "<system><normal_modes><open_paths>[0]</open_paths></normal_modes>"
    )
    s = Simulation.load_from_xml(xml, read_only=True).syslist[0]
    q = s.beads.q[:, :3].copy()
    springs = np.sum((q - np.roll(q, 1, axis=0)) ** 2)
    expected = 6 * s.ensemble.temp - 0.5 * 2 * s.nm.omegan2 * springs / 4
    assert s.properties.get_kintd("0") == pytest.approx(expected)


@pytest.mark.parametrize(
    "encoding",
    [
        "json",
        pytest.param(
            "legacy",
            marks=pytest.mark.xfail(
                strict=True,
                raises=AssertionError,
                reason="Legacy dipole parsing returns after the first bead, leaving later beads as NaN",
            ),
        ),
    ],
)
@pytest.mark.parametrize("bead", ["", "2"])
def test_driven_field_dipole_and_born_charge_properties(
    property_system, bead, encoding
):
    from ipi.engine.motion.driven_dynamics import DrivenDynamics, ElectricField, BEC

    s = property_system
    dipoles = s.beads.q[:, :3] - s.beads.q[:, 3:]
    charges = np.arange(72, dtype=float).reshape(4, 6, 3)
    motion = DrivenDynamics(
        timestep=0.01,
        mode="eda-nve",
        efield=ElectricField(amp=np.array([1.0, 2.0, 3.0]), freq=2.0),
        bec=BEC(bec=charges),
    )
    motion.Electric_Dipole.bind(s.ensemble)
    import json

    harmonic = s.forces.ff["harmonic"].driver

    def charged_harmonic(cell, q):
        energy, force, virial, _ = harmonic(cell, q)
        dipole = q[0] - q[1]
        extra = (
            json.dumps({"dipole": dipole.tolist()})
            if encoding == "json"
            else "Total dipole moment " + " ".join(str(value) for value in dipole)
        )
        return energy, force, virial, extra

    s.forces.ff["harmonic"].driver = charged_harmonic
    s.motion = s.properties.motion = motion
    s.ensemble.time = np.pi / 6
    field = np.array([1.0, 2.0, 3.0]) * 0.5
    np.testing.assert_allclose(s.properties.get_Efield(), field, atol=1e-14)
    np.testing.assert_array_equal(
        s.properties.get_dipole(bead),
        dipoles.mean(axis=0) if bead == "" else dipoles[2],
    )
    for axis, index in zip("xyz", range(3)):
        np.testing.assert_array_equal(s.trajs._get_bec(axis), charges[:, :, index])


@pytest.mark.parametrize("method", ["get_Efield", "get_dipole"])
def test_driven_properties_request_exit_for_ordinary_dynamics(
    property_system, monkeypatch, method
):
    import ipi.engine.properties as module
    from types import SimpleNamespace

    messages = []

    def stop(**kw):
        messages.append(kw["message"])
        raise SystemExit()

    monkeypatch.setattr(module, "softexit", SimpleNamespace(trigger=stop))
    with pytest.raises(SystemExit):
        getattr(property_system.properties, method)()
    assert "only when using" in messages[0]
    assert np.isnan(property_system.trajs._get_bec("x")).all()


@pytest.mark.parametrize("method", ["get_isotope_zetatd", "get_isotope_zetasc"])
def test_trajectory_isotope_rejects_invalid_atom(property_system, method):
    with pytest.raises(IndexError):
        getattr(property_system.trajs, method)(atom="2")


@pytest.mark.parametrize("ni,nj", [("2", "0"), ("0", "2")])
def test_off_diagonal_kinetic_tensor_rejects_invalid_atoms(property_system, ni, nj):
    with pytest.raises(IndexError):
        property_system.properties.get_kij(ni, nj)


def test_output_property_units_and_arguments_are_preserved(property_system):
    s = property_system
    value, dimension, unit = s.properties["potential{electronvolt}(bead=2)"]
    assert value == pytest.approx(0.25 * np.sum(s.beads.q[2] ** 2))
    assert (dimension, unit) == ("energy", "electronvolt")
    value, dimension, unit = s.trajs["positions{angstrom}"]
    np.testing.assert_array_equal(value, s.beads.q)
    assert (dimension, unit) == ("length", "angstrom")
    from ipi.engine.properties import getall

    with pytest.raises(ValueError, match="argument list"):
        getall("potential(bead=2")


@pytest.mark.parametrize(
    "method,arguments,error",
    [
        ("get_atom_ids", {"atom": "2"}, IndexError),
        ("get_kinmd", {"nm": "invalid"}, ValueError),
    ],
)
def test_kinetic_selection_rejects_invalid_arguments(
    property_system, method, arguments, error
):
    with pytest.raises(error):
        getattr(property_system.properties, method)(**arguments)


def test_kinetic_selection_with_no_matching_atoms(property_system):
    assert property_system.properties.get_kinmd(atom="missing", return_count=True) == (
        0,
        0,
    )


def test_bead_kinetic_energy_rejects_scaled_normal_mode_masses():
    from ipi.engine.simulation import Simulation

    xml = harmonic_input(nbeads=4).replace(
        "<system>",
        '<system><normal_modes><frequencies style="pa-cmd">[1]</frequencies></normal_modes>',
    )
    s = Simulation.load_from_xml(xml, read_only=True).syslist[0]
    with pytest.raises(ValueError, match="mass matrix"):
        s.properties.get_kinmd(bead="0")

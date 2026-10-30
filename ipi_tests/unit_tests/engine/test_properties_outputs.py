"""Property and trajectory outputs: geometry, per-atom values, units and parsing."""

import json
from types import SimpleNamespace

import numpy as np
import pytest

import ipi.engine.properties as properties_module
from ipi.engine.motion.driven_dynamics import BEC, DrivenDynamics, ElectricField
from ipi.engine.properties import getall, help_latex, help_rst
from ipi.engine.simulation import Simulation
from ipi_tests.unit_tests.engine.inputs import harmonic_input
from ipi_tests.unit_tests.engine.mathematical_references import polynomial


@pytest.fixture
def harmonic_system():
    """Read-only two-atom, four-bead system in U = 0.25 |q|^2."""
    return Simulation.load_from_xml(harmonic_input(nbeads=4), read_only=True).syslist[0]


@pytest.mark.parametrize("count", [2, 6])
def test_cell_and_density_outputs(reference_system, count):
    s, r = reference_system(count)
    assert s.properties["density"][0] == pytest.approx(5 / 120)
    assert s.properties["volume"][0] == pytest.approx(120)
    np.testing.assert_allclose(s.properties["cell_h"][0], [6, 5, 4, 1, 0.5, 0.7])
    lengths = np.linalg.norm(r.h, axis=0)
    angles = [
        np.degrees(np.arccos(np.dot(r.h[:, i], r.h[:, j]) / (lengths[i] * lengths[j])))
        for i, j in [(1, 2), (0, 2), (0, 1)]
    ]
    np.testing.assert_allclose(s.properties["cell_abcABC"][0], np.r_[lengths, angles])


@pytest.mark.parametrize("atom", [0, 1])
def test_atom_and_centroid_outputs(reference_system, atom):
    s, r = reference_system()
    _, force, _ = polynomial(r.q, r.a, r.b)
    columns = slice(3 * atom, 3 * atom + 3)
    velocity = r.p / np.repeat(r.masses, 3)
    np.testing.assert_allclose(
        s.properties[f"atom_x_path({atom})"][0], r.q[:, columns].ravel()
    )
    np.testing.assert_allclose(
        s.properties[f"atom_f_path({atom})"][0], force[:, columns].ravel()
    )
    for key, values in [("atom_v", velocity), ("atom_p", r.p), ("atom_f", force)]:
        for bead in [-1, 2]:
            expected = (
                values[:, columns].mean(axis=0) if bead < 0 else values[bead, columns]
            )
            np.testing.assert_allclose(
                s.properties[f"{key}({atom};{bead})"][0], expected
            )
    for key, values in [
        ("velocities", velocity),
        ("v_centroid", velocity.mean(axis=0)),
        ("x_centroid", r.q.mean(axis=0)),
        ("p_centroid", r.p.mean(axis=0)),
        ("f_centroid", force.mean(axis=0)),
    ]:
        np.testing.assert_allclose(s.trajs[key][0], values)
    for parity, selected in [
        ("even", slice(None, None, 2)),
        ("odd", slice(1, None, 2)),
    ]:
        np.testing.assert_allclose(
            s.trajs["x_centroid_" + parity][0], r.q[selected].mean(axis=0)
        )
        np.testing.assert_allclose(
            s.trajs["v_centroid_" + parity][0], velocity[selected].mean(axis=0)
        )


@pytest.mark.parametrize("species,indices", [("", [0, 1]), ("H", [0]), ("O", [1])])
def test_centroid_center_of_mass_velocity(harmonic_system, species, indices):
    s = harmonic_system
    momenta = np.arange(24).reshape(4, 2, 3) / 10
    s.beads.p_[:] = momenta.reshape(4, 6)
    expected = (
        momenta.mean(axis=0)[indices].sum(axis=0) / np.array([2.0, 3.0])[indices].sum()
    )
    np.testing.assert_allclose(s.properties.get_vcom(species), expected)


@pytest.mark.xfail(
    strict=True,
    raises=TypeError,
    reason="Bead-specific vcom indexes the Atoms object instead of its momenta",
)
def test_single_bead_center_of_mass_velocity(harmonic_system):
    s = harmonic_system
    expected = s.beads.p[2].reshape(2, 3).sum(axis=0) / 5
    np.testing.assert_allclose(s.properties.get_vcom(bead=2), expected)


@pytest.mark.parametrize("time", [0.1, 0.7, 1.2])
def test_driven_field_envelope_and_forces(reference_system, time):
    s, r = reference_system()
    charge = np.arange(r.count * 6 * 3).reshape(r.count, 6, 3) / 17
    motion = DrivenDynamics(
        timestep=0.01,
        mode="eda-nve",
        efield=ElectricField(
            amp=np.array([0.2, -0.3, 0.4]), freq=0.7, phase=0.2, peak=0.6, sigma=0.3
        ),
        bec=BEC(bec=charge),
    )
    motion.bind(
        s.ensemble, s.beads, s.nm, s.cell, s.forces, s.simul.prng, s.simul.output_maker
    )
    s.motion = s.properties.motion = motion
    s.ensemble.time = time
    envelope = np.exp(-0.5 * ((time - 0.6) / 0.3) ** 2)
    field = envelope * np.cos(0.7 * time + 0.2) * np.array([0.2, -0.3, 0.4])
    assert s.properties["Eenvelope"][0] == pytest.approx(envelope)
    np.testing.assert_allclose(s.trajs["Eforces"][0], charge @ field, atol=1e-13)


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
def test_driven_field_dipole_and_born_charges(harmonic_system, bead, encoding):
    s = harmonic_system
    dipoles = s.beads.q[:, :3] - s.beads.q[:, 3:]
    charges = np.arange(72, dtype=float).reshape(4, 6, 3)
    motion = DrivenDynamics(
        timestep=0.01,
        mode="eda-nve",
        efield=ElectricField(amp=np.array([1.0, 2.0, 3.0]), freq=2.0),
        bec=BEC(bec=charges),
    )
    motion.Electric_Dipole.bind(s.ensemble)
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
    np.testing.assert_allclose(
        s.properties.get_Efield(), np.array([1.0, 2.0, 3.0]) * 0.5, atol=1e-14
    )
    np.testing.assert_array_equal(
        s.properties.get_dipole(bead),
        dipoles.mean(axis=0) if bead == "" else dipoles[2],
    )
    for index, axis in enumerate("xyz"):
        np.testing.assert_array_equal(s.trajs._get_bec(axis), charges[:, :, index])


@pytest.mark.parametrize("method", ["get_Efield", "get_dipole"])
def test_driven_properties_require_driven_dynamics(
    harmonic_system, monkeypatch, method
):
    messages = []

    def stop(**kw):
        messages.append(kw["message"])
        raise SystemExit()

    monkeypatch.setattr(properties_module, "softexit", SimpleNamespace(trigger=stop))
    with pytest.raises(SystemExit):
        getattr(harmonic_system.properties, method)()
    assert "only when using" in messages[0]
    assert np.isnan(harmonic_system.trajs._get_bec("x")).all()


def test_output_units_and_arguments(harmonic_system):
    s = harmonic_system
    value, dimension, unit = s.properties["potential{electronvolt}(bead=2)"]
    assert value == pytest.approx(0.25 * np.sum(s.beads.q[2] ** 2))
    assert (dimension, unit) == ("energy", "electronvolt")
    value, dimension, unit = s.trajs["positions{angstrom}"]
    np.testing.assert_array_equal(value, s.beads.q)
    assert (dimension, unit) == ("length", "angstrom")
    with pytest.raises(ValueError, match="argument list"):
        getall("potential(bead=2")


@pytest.mark.parametrize("text", ["potential{electronvolt", "potential}"])
@pytest.mark.xfail(
    strict=True,
    raises=pytest.fail.Exception,
    reason="Property parser accepts unmatched unit braces",
)
def test_parser_rejects_unmatched_unit_braces(text):
    with pytest.raises(ValueError):
        getall(text)


@pytest.mark.parametrize("standalone", [False, True])
def test_property_help_text(standalone):
    metadata = {
        "a_value": {
            "help": "short",
            "longhelp": "Long description",
            "dimension": "energy",
            "size": 3,
        },
        "z_value": {"help": "Short description", "dimension": "undefined"},
    }
    latex = help_latex(metadata, standalone)
    rst = help_rst(metadata, standalone)
    for rendered in [latex, rst]:
        assert "Long description" in rendered
        assert "Short description" in rendered
        assert "dimension: energy" in rendered
        assert "size: 3" in rendered
        assert "undefined" not in rendered
    assert latex.index(r"a\_value") < latex.index(r"z\_value")
    assert rst.index("a_value") < rst.index("z_value")
    assert (r"\begin{document}" in latex) == standalone
    assert (r"\end{document}" in latex) == standalone

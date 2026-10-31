"""Driven dynamics: field shape, field-driven motion and client response data."""

import xml.etree.ElementTree as ET
from types import SimpleNamespace

import numpy as np
import pytest

from ipi.engine.beads import Beads
from ipi.engine.motion.driven_dynamics import BEC, ElectricDipole, ElectricField
from ipi.utils.depend import depend_value
from ipi.utils.units import unit_to_internal
from ipi_tests.unit_tests.engine.inputs import harmonic_input


@pytest.mark.parametrize("pulsed", [False, True])
def test_electric_field_scalar_vector_and_disabled_values(pulsed):
    amp = np.array([0.2, -0.3, 0.4])
    t = np.array([0.0, 0.5, 1.0, 1.5])
    field = ElectricField(
        amp=amp, freq=0.7, phase=0.2, peak=1.0, sigma=0.4 if pulsed else np.inf
    )
    expected = np.cos(0.7 * t + 0.2)
    if pulsed:
        expected *= np.exp(-0.5 * ((t - 1) / 0.4) ** 2)
    np.testing.assert_allclose(field.Efield(t), expected[:, None] * amp)
    for i, time in enumerate(t):
        np.testing.assert_allclose(field.Efield(time), expected[i] * amp)
    field.enabled = False
    np.testing.assert_array_equal(field.Efield(t), np.zeros((4, 3)))
    np.testing.assert_array_equal(field.Efield(0.5), np.zeros(3))


def test_charges_in_constant_field_accelerate_uniformly(simulation):
    root = ET.fromstring(
        harmonic_input(timestep=0.01).replace("{k1: 0.5}", "{k1: 0.0}")
    )
    motion = root.find("system/motion")
    motion.set("mode", "driven_dynamics")
    dynamics = motion.find("dynamics")
    dynamics.tag = "driven_dynamics"
    dynamics.set("mode", "eda-nve")
    field = ET.SubElement(dynamics, "efield")
    ET.SubElement(field, "amp", shape="(3)").text = "[.2,-.1,.3]"
    charges = np.vstack([np.eye(3), -np.eye(3)])
    ET.SubElement(dynamics, "bec", mode="manual", shape="(6,3)").text = str(
        charges.ravel().tolist()
    )
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        s = sim.syslist[0]
        q0 = s.beads.q.copy()
        p0 = s.beads.p.copy()
        m = s.beads.m3.copy()
        force = np.array([[0.2, -0.1, 0.3, -0.2, 0.1, -0.3]])
        sim.run(steps=10, write_outputs=False)
        np.testing.assert_allclose(s.beads.p, p0 + 0.1 * force, atol=1e-12)
        np.testing.assert_allclose(
            s.beads.q, q0 + 0.1 * p0 / m + 0.5 * 0.1**2 * force / m, atol=1e-12
        )
        expected = s.beads.p[0, :3] / 2 - s.beads.p[0, 3:] / 3
        np.testing.assert_allclose(s.motion.integrator.dipole_dtime(0), expected)
        np.testing.assert_allclose(
            s.motion.integrator.dipole_dtime(), expected[None, :]
        )
        with pytest.raises(ValueError, match="outermost"):
            s.motion.integrator.td_pstep(0, 1)


def response_objects(extras):
    """BEC and dipole objects that read the given force-client extras for two beads."""
    ensemble = SimpleNamespace(beads=Beads(2, 2), forces=SimpleNamespace(extras=extras))
    bec = BEC(cbec=True)
    bec.bind(ensemble, "eda-nve", depend_value(name="threshold", value=1e-8))
    dipole = ElectricDipole()
    dipole.bind(ensemble)
    return bec, dipole


@pytest.mark.parametrize(
    "data,match",
    [
        ({}, "not returned"),
        ({"BEC": [np.zeros((3, 3))] * 2}, "number of BEC"),
        ({"BEC": [np.zeros((6, 2))] * 2}, "wrong shape"),
        ({"BEC": [np.ones((6, 3))] * 2}, "charge conservation"),
    ],
)
def test_invalid_client_charges_are_rejected(data, match):
    bec, _ = response_objects(data)
    with pytest.raises(ValueError, match=match):
        _ = bec.bec


@pytest.mark.parametrize("index", [-1, 2])
def test_response_bead_selection_rejects_invalid_indices(index):
    bec, dipole = response_objects({})
    with pytest.raises(ValueError, match="bead"):
        bec._get_driver_BEC(index)
    with pytest.raises(ValueError, match="bead"):
        dipole._get_dipole(index)


def test_structured_dipoles_preserve_bead_order():
    values = [[0.1, -0.2, 0.3], [0.4, 0.5, -0.6]]
    _, dipole = response_objects({"dipole": values})
    np.testing.assert_allclose(dipole.dipole, values)


@pytest.mark.parametrize(
    "extras", [{}, {"dipole": [[1, 2], [3, 4]]}, {"raw": ["invalid", "invalid"]}]
)
def test_missing_or_malformed_dipoles_cannot_look_like_zero(extras):
    _, dipole = response_objects(extras)
    assert np.isnan(dipole.dipole).all()


def test_legacy_dipole_text_converts_angstrom_to_atomic_units():
    _, dipole = response_objects(
        {"raw": ["Total dipole moment [eAng]: 1.0 -2.0 3.0"] * 2}
    )
    expected = np.tile([1, -2, 3], (2, 1)) * unit_to_internal("length", "angstrom", 1)
    np.testing.assert_allclose(dipole.dipole, expected)

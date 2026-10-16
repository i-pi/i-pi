"""Electric-field propagation checked against constant-acceleration motion."""

import xml.etree.ElementTree as ET
import numpy as np
import pytest
from ipi.engine.motion.driven_dynamics import ElectricField
from ipi_tests.unit_tests.engine.test_simulation_invariants import harmonic_input


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


def test_opposite_charges_in_constant_field_follow_exact_acceleration(simulation):
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

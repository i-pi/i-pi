"""Pressure dynamics must preserve the complete trajectory through a restart."""

import xml.etree.ElementTree as ET
import numpy as np
import pytest
from ipi_tests.unit_tests.engine.test_simulation_invariants import harmonic_input


def pressure_input(barostat, nbeads=1, constant_volume=False):
    root = ET.fromstring(
        harmonic_input(nbeads=nbeads, thermostat="langevin", timestep=0.001)
    )
    ensemble = root.find("system/ensemble")
    ET.SubElement(ensemble, "pressure").text = ".000001"
    dynamics = root.find("system/motion/dynamics")
    dynamics.set("mode", "nst" if barostat == "anisotropic" else "npt")
    if barostat == "anisotropic":
        ET.SubElement(ensemble, "stress", shape="(3,3)").text = str(
            (np.eye(3) * 0.000001).ravel().tolist()
        )
    baro = ET.SubElement(dynamics, "barostat", mode=barostat)
    ET.SubElement(baro, "tau").text = "100"
    thermo = ET.SubElement(baro, "thermostat", mode="langevin")
    ET.SubElement(thermo, "tau").text = "100"
    if barostat == "anisotropic":
        ET.SubElement(baro, "h0", shape="(3,3)").text = str(
            (np.eye(3) * 20).ravel().tolist()
        )
    if constant_volume:
        ET.SubElement(baro, "vol_constraint").text = "True"
    return ET.tostring(root, encoding="unicode")


def state(sim):
    s = sim.syslist[0]
    return [
        s.beads.q.copy(),
        s.beads.p.copy(),
        s.cell.h.copy(),
        np.asarray(s.motion.barostat.p).copy(),
        s.motion.barostat.ebaro,
    ]


@pytest.mark.parametrize("barostat", ["isotropic", "flexible", "anisotropic"])
@pytest.mark.parametrize("nbeads", [1, 4])
def test_pressure_restart_preserves_cell_momenta_and_energy(
    restart_states, barostat, nbeads
):
    xml = pressure_input(barostat, nbeads)
    expected, actual = restart_states(xml, state)
    for a, b in zip(actual, expected):
        np.testing.assert_allclose(a, b, rtol=1e-10, atol=1e-12)
    assert np.linalg.det(actual[2]) > 0


def test_flexible_barostat_volume_constraint(simulation):
    with simulation(pressure_input("flexible", constant_volume=True)) as sim:
        initial = sim.syslist[0].cell.V
        sim.run(steps=8, write_outputs=False)
        assert sim.syslist[0].cell.V == pytest.approx(initial, rel=1e-12)

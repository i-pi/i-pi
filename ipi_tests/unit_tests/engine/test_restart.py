"""A run restarted from its checkpoint must continue exactly like an uninterrupted run."""

import xml.etree.ElementTree as ET

import numpy as np
import pytest

from ipi_tests.unit_tests.engine.inputs import harmonic_input, pressure_input


def thermostat_input(mode, nbeads, explicit_covariance=False):
    root = ET.fromstring(harmonic_input(nbeads=nbeads, thermostat=mode))
    bath = root.find("system/motion/dynamics/thermostat")
    if mode and "gle" in mode:
        drift = np.array([[2.0, -0.2], [0.2, 1.0]])
        covariance = np.eye(2) * 0.004
        if mode != "gle":
            drift = np.tile(drift, (nbeads, 1, 1))
            covariance = np.tile(covariance, (nbeads, 1, 1))
        ET.SubElement(bath, "A", shape=str(drift.shape)).text = str(
            drift.ravel().tolist()
        )
        if explicit_covariance:
            ET.SubElement(bath, "C", shape=str(covariance.shape)).text = str(
                covariance.ravel().tolist()
            )
    if mode == "multi":
        for child_mode in ["langevin", "svr"]:
            child = ET.SubElement(bath, "thermostat", mode=child_mode)
            ET.SubElement(child, "tau").text = "0.5"
    return ET.tostring(root, encoding="unicode")


def assert_same(actual, expected):
    for a, b in zip(actual, expected):
        np.testing.assert_allclose(a, b, rtol=1e-10, atol=1e-12)


@pytest.mark.parametrize(
    "mode,nbeads,explicit_covariance",
    [
        (m, n, False)
        for m in [None, "langevin", "svr", "pile_l", "pile_g"]
        for n in [1, 4]
    ]
    + [(m, 4, c) for m in ["gle", "nm_gle", "nm_gle_g"] for c in [False, True]]
    + [(m, 4, False) for m in ["multi", "cl", "ffl"]],
)
def test_thermostat_restart(restart_states, mode, nbeads, explicit_covariance):
    """Positions, momenta, heat and auxiliary momenta use the same random draws."""

    def snapshot(sim):
        system = sim.syslist[0]
        bath = system.motion.thermostat
        values = [system.beads.q.copy(), system.beads.p.copy(), bath.ethermo]
        if mode and "gle" in mode:
            values.append(bath.s.copy())
        return values + [system.ensemble.time]

    xml = thermostat_input(mode, nbeads, explicit_covariance)
    expected, actual = restart_states(xml, snapshot, steps=(3, 3))
    assert_same(actual, expected)


@pytest.mark.parametrize("barostat", ["isotropic", "flexible", "anisotropic"])
@pytest.mark.parametrize("nbeads", [1, 4])
def test_barostat_restart(restart_states, barostat, nbeads):
    def snapshot(sim):
        s = sim.syslist[0]
        return [
            s.beads.q.copy(),
            s.beads.p.copy(),
            s.cell.h.copy(),
            np.asarray(s.motion.barostat.p).copy(),
            s.motion.barostat.ebaro,
        ]

    expected, actual = restart_states(pressure_input(barostat, nbeads), snapshot)
    assert_same(actual, expected)
    assert np.linalg.det(actual[2]) > 0


@pytest.mark.parametrize("mode", ["sc", "scnpt", "nvt-cc"])
@pytest.mark.parametrize("splitting", ["obabo", "baoab"])
def test_higher_order_and_constrained_centroid_restart(
    simulation, restart_states, mode, splitting
):
    root = ET.fromstring(
        pressure_input("isotropic", 4)
        if mode == "scnpt"
        else harmonic_input(nbeads=4, thermostat="pile_l", timestep=0.001)
    )
    dynamics = root.find("system/motion/dynamics")
    dynamics.set("mode", mode)
    dynamics.set("splitting", splitting)
    if mode == "scnpt":
        dynamics.find("barostat").set("mode", "sc-isotropic")
    xml = ET.tostring(root, encoding="unicode")

    def snapshot(sim):
        s = sim.syslist[0]
        return [s.beads.q.copy(), s.beads.p.copy(), s.cell.h.copy(), s.ensemble.econs]

    expected, actual = restart_states(xml, snapshot)
    assert_same(actual, expected)
    if mode == "nvt-cc":
        with simulation(xml) as sim:
            initial_centroid = sim.syslist[0].beads.qc.copy()
        np.testing.assert_allclose(
            expected[0].mean(axis=0), initial_centroid, atol=1e-13
        )


@pytest.mark.parametrize("kind,attribute", [("t", "temp"), ("p", "pext")])
@pytest.mark.parametrize("logscale", [False, True])
def test_ramp_restart(simulation, restart_states, kind, attribute, logscale):
    root = ET.fromstring(harmonic_input(timestep=0.01))
    system = root.find("system")
    dynamics = system.find("motion")
    system.remove(dynamics)
    multi = ET.SubElement(system, "motion", mode="multi")
    ramp = ET.SubElement(multi, "motion", mode=kind + "_ramp")
    options = ET.SubElement(ramp, kind + "_ramp")
    for key, value in {
        kind + "_start": 0.1,
        kind + "_end": 1.6,
        "total_steps": 4,
        "logscale": str(logscale),
    }.items():
        ET.SubElement(options, key).text = str(value)
    multi.append(dynamics)
    xml = ET.tostring(root, encoding="unicode")
    # Four steps from 0.1 to 1.6, then constant.
    schedule = [0.2, 0.4, 0.8, 1.6] if logscale else [0.475, 0.85, 1.225, 1.6]
    with simulation(xml) as sim:
        for expected in schedule + [1.6, 1.6]:
            sim.run(steps=1, write_outputs=False)
            assert getattr(sim.syslist[0].ensemble, attribute) == pytest.approx(
                expected
            )

    def snapshot(sim):
        s = sim.syslist[0]
        return [
            s.beads.q.copy(),
            s.beads.p.copy(),
            getattr(s.ensemble, attribute),
            s.ensemble.time,
        ]

    expected, actual = restart_states(xml, snapshot)
    assert_same(actual, expected)
    assert actual[2] == pytest.approx(1.6)
    assert actual[3] == pytest.approx(0.06)

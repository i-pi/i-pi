"""Self-consistent phonons for a separable harmonic potential, for which the exact
dynamical matrix is a fixed point and the minimum is found in one iteration."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

import json
import xml.etree.ElementTree as ET

import numpy as np
import pytest
from scipy.optimize import brentq

from ipi.utils.array_backend import to_numpy
from ipi_tests.unit_tests.engine.inputs import motion_input

# U = 0.25 |q|^2 for two atoms of masses 2 and 3 (see inputs.harmonic_input)
MASSES = np.repeat([2.0, 3.0], 3)
FREQUENCIES = np.sqrt(0.5 / MASSES)
MAX_STEPS = 8
NPARALLEL = 2
# reset, force evaluations, and displacement of the mean position
STEPS_PER_ITERATION = 1 + MAX_STEPS // NPARALLEL + 1


def scp_input(displace_mode, mode, scale=1.0):
    """Harmonic input with SCP starting from the exact dynamical matrix, with the
    initial positions of the harmonic input multiplied by scale."""

    options = ET.Element("scp", mode=mode)
    for key, value in {
        "displace_mode": displace_mode,
        "max_steps": MAX_STEPS,
        "max_iter": 2,
        "nparallel": NPARALLEL,
        "random_type": "pseudo",
        "asr": "none",
    }.items():
        ET.SubElement(options, key).text = str(value)
    ET.SubElement(options, "dynmat", shape="(6,6)").text = str(
        np.diag(0.5 / MASSES).ravel().tolist()
    )
    root = motion_input("scp", options)
    q = root.find("system/beads/q")
    q.text = str((scale * np.array(json.loads(q.text))).tolist())
    return ET.tostring(root, encoding="unicode")


@pytest.mark.parametrize("scale", [1.0, 1.3, 2.0])
@pytest.mark.parametrize(
    "displace_mode,mode", [("ik", "qn"), ("ik", "cl"), ("nmik", "qn"), ("nmik", "cl")]
)
def test_harmonic_dynamical_matrix_is_a_fixed_point(
    simulation, displace_mode, mode, scale
):
    """The mean position moves to the minimum after the first iteration, and
    the frequencies remain the harmonic ones, also when the samples of the
    first iteration have a negligible weight at the minimum"""

    with simulation(scp_input(displace_mode, mode, scale)) as sim:
        scp = sim.syslist[0].motion
        assert np.abs(to_numpy(scp.beads.q)).max() > 0.1
        for iteration in [1, 2]:
            sim.run(steps=STEPS_PER_ITERATION, write_outputs=False)
            assert scp.isc == iteration
            np.testing.assert_allclose(to_numpy(scp.beads.q), 0.0, atol=1e-12)
            np.testing.assert_allclose(
                np.sort(to_numpy(scp.w)), np.sort(FREQUENCIES), rtol=1e-10
            )
        # the error of the force from the samples drawn at the minimum
        _, f_err, _ = scp.phononator.weighted_force()
        assert np.isfinite(to_numpy(f_err)).all()


def test_steepest_descent_moves_the_mean_position(simulation):
    """With the default check of the significance of the forces, "sd" moves the
    mean position by tau D <f>, with D the displacement correlation"""

    with simulation(scp_input("sd", "qn")) as sim:
        scp = sim.syslist[0].motion
        q = to_numpy(scp.beads.q).ravel().copy()
        sim.run(steps=STEPS_PER_ITERATION, write_outputs=False)
        assert scp.isc == 1
        # the harmonic average force is -K q, and D and K are unchanged
        displacement = scp.tau * to_numpy(scp.D) @ to_numpy(scp.K) @ q
        np.testing.assert_allclose(
            to_numpy(scp.beads.q).ravel(), q - displacement, rtol=1e-10
        )


class DoubleWell:
    """Driver for U = sum(-a x^2 + b x^4) over the Cartesian coordinates."""

    def __init__(self, a, b):
        self.a, self.b = a, b

    def __call__(self, cell, positions):
        x = positions.reshape(-1)
        energy = np.sum(-self.a * x**2 + self.b * x**4)
        force = 2 * self.a * x - 4 * self.b * x**3
        return energy, force.reshape(positions.shape), np.zeros((3, 3)), ""


@pytest.mark.parametrize("mode", ["cl", "qn"])
def test_double_well_converges_to_the_optimal_gaussian(simulation, mode):
    """For U = -a x^2 + b x^4, the Gaussian of mean x0 and width s^2 that
    minimizes the free energy has curvature K = <U''> = 4a - 24 b s^2 and
    x0^2 = a / 2b - 3 s^2, with s^2 = kT / K (classical) or
    coth(w / 2kT) / 2 m w, w^2 = K / m (quantum)"""

    a, b, temperature = 0.5, 0.05, 0.3
    max_steps, iterations = 1024, 10
    options = ET.Element("scp", mode=mode)
    for key, value in {
        "displace_mode": "nmik",
        "max_steps": max_steps,
        "max_iter": iterations,
        "nparallel": 64,
        "random_type": "pseudo",
        "asr": "none",
    }.items():
        ET.SubElement(options, key).text = str(value)
    # starts from the minima, with the exact harmonic dynamical matrix there
    ET.SubElement(options, "dynmat", shape="(6,6)").text = str(
        np.diag(4 * a / MASSES).ravel().tolist()
    )
    root = motion_input("scp", options)
    root.find("system/beads/q").text = str([float(np.sqrt(a / (2 * b)))] * 6)
    temp = root.find("system/ensemble/temperature")
    temp.attrib.clear()
    temp.text = str(temperature)

    def width(K, m):
        if mode == "cl":
            return temperature / K
        w = np.sqrt(K / m)
        return 1 / np.tanh(w / (2 * temperature)) / (2 * m * w)

    # the root that tends to the harmonic curvature 4a at low temperature
    K = np.array(
        [
            brentq(lambda K: K - 4 * a + 24 * b * width(K, m), 2 * a, 4 * a)
            for m in MASSES
        ]
    )
    s2 = np.array([width(k, m) for k, m in zip(K, MASSES)])
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        sim.syslist[0].forces.ff["harmonic"].driver = DoubleWell(a, b)
        scp = sim.syslist[0].motion
        sim.run(steps=iterations * (max_steps // 64 + 2), write_outputs=False)
        assert scp.isc == iterations
        np.testing.assert_allclose(
            to_numpy(scp.beads.q).ravel(), np.sqrt(a / (2 * b) - 3 * s2), rtol=1e-2
        )
        # the frequencies are sorted, the coordinates are not
        np.testing.assert_allclose(
            np.sort(to_numpy(scp.w)), np.sort(np.sqrt(K / MASSES)), rtol=1e-2
        )

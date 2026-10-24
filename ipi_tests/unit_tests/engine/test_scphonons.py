"""Self-consistent phonons for a separable harmonic potential, for which the exact
dynamical matrix is a fixed point and the minimum is found in one iteration."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

import xml.etree.ElementTree as ET

import numpy as np
import pytest

from ipi.utils.array_backend import to_numpy
from ipi_tests.unit_tests.engine.inputs import motion_input

# U = 0.25 |q|^2 for two atoms of masses 2 and 3 (see inputs.harmonic_input)
MASSES = np.repeat([2.0, 3.0], 3)
FREQUENCIES = np.sqrt(0.5 / MASSES)
MAX_STEPS = 8
NPARALLEL = 2
# reset, force evaluations, and displacement of the mean position
STEPS_PER_ITERATION = 1 + MAX_STEPS // NPARALLEL + 1


def scp_input(displace_mode, mode):
    """Harmonic input with SCP starting from the exact dynamical matrix."""

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
    return ET.tostring(motion_input("scp", options), encoding="unicode")


# "ik" estimates the Hessian after moving the mean position, by reweighting the
# samples, which with classical statistics are too narrow to have any weight there
@pytest.mark.parametrize(
    "displace_mode,mode", [("ik", "qn"), ("nmik", "qn"), ("nmik", "cl")]
)
def test_harmonic_dynamical_matrix_is_a_fixed_point(simulation, displace_mode, mode):
    """The mean position moves to the minimum after the first iteration, and
    the frequencies remain the harmonic ones"""

    with simulation(scp_input(displace_mode, mode)) as sim:
        scp = sim.syslist[0].motion
        assert np.abs(to_numpy(scp.beads.q)).max() > 0.1
        for iteration in [1, 2]:
            sim.run(steps=STEPS_PER_ITERATION, write_outputs=False)
            assert scp.isc == iteration
            np.testing.assert_allclose(to_numpy(scp.beads.q), 0.0, atol=1e-12)
            np.testing.assert_allclose(
                np.sort(to_numpy(scp.w)), np.sort(FREQUENCIES), rtol=1e-10
            )

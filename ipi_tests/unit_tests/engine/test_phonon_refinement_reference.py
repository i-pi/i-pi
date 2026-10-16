"""Finite differences must recover the force constants, even from a poor guess."""

import xml.etree.ElementTree as ET

import numpy as np
import pytest

from ipi_tests.unit_tests.engine.test_mode_and_path_reference import motion_input


@pytest.mark.parametrize(
    "mode",
    [
        "fd",
        *[
            pytest.param(
                mode,
                marks=pytest.mark.xfail(
                    strict=True,
                    reason="NMFD transform discards the refined matrix and restores the input matrix",
                ),
            )
            for mode in ("nmfd", "enmfd")
        ],
    ],
)
def test_phonon_refinement_recovers_harmonic_force_constants(simulation, mode):
    options = ET.Element("vibrations", mode=mode)
    masses = np.repeat([2.0, 3.0], 3)
    ET.SubElement(options, "dynmat", shape="(6,6)").text = str(
        np.diag(0.8 / masses).ravel().tolist()
    )
    ET.SubElement(options, "pos_shift").text = ".001"
    ET.SubElement(options, "energy_shift").text = ".0001"
    root = motion_input("vibrations", options)
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        motion = sim.syslist[0].motion
        motion.exit_on_convergence = False
        sim.run(steps=7, write_outputs=False)
        np.testing.assert_allclose(
            motion.refdynmatrix, np.diag(0.5 / masses), atol=1e-10
        )

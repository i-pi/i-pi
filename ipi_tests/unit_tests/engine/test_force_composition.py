"""Weighted and centroid-contracted forces preserve energies and cache invalidation."""

import xml.etree.ElementTree as ET

import numpy as np
import pytest
from ipi_tests.unit_tests.engine.test_simulation_invariants import harmonic_input


@pytest.mark.parametrize("weight", [0.0, 0.3, -0.2])
def test_weighted_centroid_force_and_saved_state(simulation, weight):
    root = ET.fromstring(harmonic_input(nbeads=4))
    forces = root.find("system/forces")
    ET.SubElement(
        forces, "force", forcefield="harmonic", nbeads="1", weight=str(weight)
    )
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        s = sim.syslist[0]
        q = s.beads.q.copy()
        centroid = q.mean(axis=0)
        expected_force = -0.5 * (q + weight * centroid)
        expected_pots = 0.25 * (np.sum(q * q, axis=1) + weight * np.sum(centroid**2))
        np.testing.assert_allclose(s.forces.f, expected_force)
        np.testing.assert_allclose(s.forces.pots, expected_pots)
        snapshot = s.forces.dump_state()
        for interpolate in [False, True]:
            count = 4 if interpolate else 1
            for weighted in [False, True]:
                w = weight if weighted else 1.0
                np.testing.assert_allclose(
                    s.forces.forces_component(1, weighted, interpolate),
                    np.tile(-0.5 * w * centroid, (count, 1)),
                )
                np.testing.assert_allclose(
                    s.forces.pots_component(1, weighted, interpolate),
                    np.full(count, 0.25 * w * np.sum(centroid**2)),
                )
                np.testing.assert_array_equal(
                    s.forces.virs_component(1, weighted, interpolate),
                    np.zeros((count, 3, 3)),
                )
        s.beads.q[:] = q * 2
        np.testing.assert_allclose(s.forces.f, 2 * expected_force)
        s.beads.q[:] = q
        s.forces.load_state(snapshot)
        np.testing.assert_allclose(s.forces.f, expected_force)
        s.beads.q[:] = q * 3
        np.testing.assert_allclose(s.forces.f, 3 * expected_force)
        np.testing.assert_allclose(s.forces.pots, 9 * expected_pots)
        with pytest.raises(ValueError, match="components"):
            s.forces.load_state(snapshot[:1])

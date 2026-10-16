"""Geometry optimization must leave prescribed coordinates unchanged."""

import xml.etree.ElementTree as ET
import numpy as np
import pytest
from ipi_tests.unit_tests.engine.test_simulation_invariants import harmonic_input


@pytest.mark.parametrize(
    "method", ["sd", "cg", "bfgs", "bfgstrm", "lbfgs", "damped_bfgs"]
)
@pytest.mark.parametrize("selection", ["fixatoms", "fixatoms_dof"])
def test_fixed_coordinates_survive_optimization(simulation, method, selection, request):
    if method == "damped_bfgs" and selection == "fixatoms_dof":
        request.node.add_marker(
            pytest.mark.xfail(
                strict=True,
                raises=ValueError,
                reason="Damped_BFGS reshapes five free coordinates into triples",
            )
        )
    root = ET.fromstring(harmonic_input())
    motion = root.find("system/motion")
    motion.clear()
    motion.set("mode", "minimize")
    ET.SubElement(motion, "fixcom").text = "False"
    ET.SubElement(motion, selection, shape="(1)").text = "[0]"
    optimizer = ET.SubElement(motion, "optimizer", mode=method)
    ET.SubElement(optimizer, "exit_on_convergence").text = "False"
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        system = sim.syslist[0]
        initial = system.beads.q.copy()
        fixed = [0, 1, 2] if selection == "fixatoms" else [0]
        free = sorted(set(range(6)) - set(fixed))
        sim.run(steps=150, write_outputs=False)
        np.testing.assert_array_equal(system.beads.q[:, fixed], initial[:, fixed])
        np.testing.assert_allclose(system.beads.q[:, free], 0, atol=2e-4)
        assert system.forces.pot == pytest.approx(
            0.25 * np.sum(initial[:, fixed] ** 2), abs=1e-7
        )

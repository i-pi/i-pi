"""Property documentation and center-of-mass velocities have explicit references."""

import numpy as np
import pytest
from ipi.engine.properties import help_latex, help_rst
from ipi_tests.unit_tests.engine.test_simulation_invariants import harmonic_input


@pytest.mark.parametrize("standalone", [False, True])
def test_property_help_preserves_descriptions_dimensions_and_sizes(standalone):
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


@pytest.mark.parametrize("species,indices", [("", [0, 1]), ("H", [0]), ("O", [1])])
def test_centroid_center_of_mass_velocity(simulation, species, indices):
    with simulation(harmonic_input(nbeads=4)) as sim:
        s = sim.syslist[0]
        momenta = np.arange(24).reshape(4, 2, 3) / 10
        s.beads.p[:] = momenta.reshape(4, 6)
        expected = (
            momenta.mean(axis=0)[indices].sum(axis=0)
            / np.array([2.0, 3.0])[indices].sum()
        )
        np.testing.assert_allclose(s.properties.get_vcom(species), expected)


@pytest.mark.xfail(
    strict=True,
    raises=TypeError,
    reason="Bead-specific vcom indexes the Atoms object instead of its momenta",
)
def test_single_bead_center_of_mass_velocity(simulation):
    with simulation(harmonic_input(nbeads=4)) as sim:
        s = sim.syslist[0]
        expected = s.beads.p[2].reshape(2, 3).sum(axis=0) / 5
        np.testing.assert_allclose(s.properties.get_vcom(bead=2), expected)

"""Mass-scaled free normal modes must follow the specified harmonic dynamics."""

import xml.etree.ElementTree as ET

import numpy as np
import pytest
from scipy.linalg import expm
from ipi_tests.unit_tests.engine.test_simulation_invariants import harmonic_input


@pytest.mark.parametrize(
    "style,freqs,propagator",
    [
        pytest.param(
            style,
            freqs,
            propagator,
            marks=(
                pytest.mark.xfail(
                    strict=True,
                    raises=AssertionError,
                    reason="Cayley propagation ignores normal-mode mass scaling",
                )
                if propagator == "cayley" and style != "rpmd"
                else []
            ),
        )
        for style, freqs in [
            ("rpmd", []),
            ("manual", [0.02, 0.03, 0.04]),
            ("pa-cmd", [0.02]),
            ("wmax-cmd", [0.1, 0.02]),
        ]
        for propagator in ["exact", "cayley"]
    ],
)
def test_open_and_closed_mode_propagators_match_harmonic_matrix(
    simulation, style, freqs, propagator
):
    root = ET.fromstring(harmonic_input(nbeads=4, timestep=0.1))
    normal = ET.SubElement(
        root.find("system"), "normal_modes", transform="matrix", propagator=propagator
    )
    ET.SubElement(normal, "frequencies", style=style).text = str(freqs)
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        nm = sim.syslist[0].nm
        for omega, factors, prop in [
            (nm.omegak, nm.nm_factor, nm.prop_pq),
            (nm.o_omegak, nm.o_nm_factor, nm.o_prop_pq),
        ]:
            for k in range(1, 4):
                generator = np.array([[0.0, -omega[k] ** 2], [1 / factors[k], 0.0]])
                if propagator == "exact":
                    expected = expm(nm.dt * generator)
                else:
                    half = nm.dt * generator / 2
                    expected = np.linalg.solve(np.eye(2) - half, np.eye(2) + half)
                np.testing.assert_allclose(prop[k], expected, atol=1e-12)
                assert np.linalg.det(prop[k]) == pytest.approx(1.0)
        cloned = nm.copy()
        assert cloned.mode == nm.mode
        np.testing.assert_array_equal(cloned.nm_freqs, nm.nm_freqs)

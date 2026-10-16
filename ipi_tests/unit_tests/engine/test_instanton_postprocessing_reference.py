"""Postprocessing formulas checked on explicitly specified quadratic configurations."""

import re
import xml.etree.ElementTree as ET
import numpy as np
import pytest
from ipi.utils.tools.instanton_postproc import instanton_compute
from ipi_tests.unit_tests.engine.test_simulation_invariants import harmonic_input


@pytest.mark.parametrize("case", ["TS", "splitting", "rate"])
def test_postprocessing_uses_saved_energy_and_spring_action(simulation, capsys, case):
    count = 1 if case == "TS" else 4
    root = ET.fromstring(harmonic_input(nbeads=count))
    temperature = root.find("system/ensemble/temperature")
    temperature.attrib.clear()
    temperature.text = ".2"
    motion = root.find("system/motion")
    motion.clear()
    motion.set("mode", "instanton")
    ET.SubElement(motion, "fixcom").text = "False"
    options = ET.SubElement(
        motion, "instanton", mode="splitting" if case == "splitting" else "rate"
    )
    q = np.tile([0.2, -0.3, 0.1, -0.4, 0.2, 0.3], (count, 1))
    q[:, 0] += 0.01 * np.arange(count)
    h = np.eye(6) * 0.5
    if case == "TS":
        h[0, 0] = -0.5
    ET.SubElement(options, "hessian", shape=f"(6,{6*count})").text = str(
        np.tile(h, (1, count)).ravel().tolist()
    )
    potential = 0.25 * np.sum(q * q, axis=1)
    ET.SubElement(options, "old_pot", shape=f"({count})").text = str(potential.tolist())
    ET.SubElement(options, "energy_shift").text = "0"
    ET.SubElement(options, "hessian_asr").text = "none"
    with open("configuration.xml", "w") as f:
        f.write(ET.tostring(root, encoding="unicode"))
    frequencies = np.sqrt(0.5 / np.repeat([2.0, 3.0], 3))
    np.savetxt("minimum.dat", frequencies[None, :])
    with pytest.raises(SystemExit) as stopped:
        instanton_compute(
            "configuration.xml",
            "TS" if case == "TS" else "instanton",
            0.2,
            asr="none",
            input_freq="minimum.dat",
        )
    assert stopped.value.code == 0
    output = capsys.readouterr().out
    assert "nan" not in output.lower()
    if case == "TS":
        value = float(re.search(r"logQvib: ([^\n]+)", output).group(1))
        expected = -np.log(2 * np.sinh(np.delete(frequencies, 0) / 0.4)).sum()
        assert value == pytest.approx(expected)
    elif case == "splitting":
        a, b, total = map(
            float,
            re.search(r"S1/hbar (\S+) ,S2/hbar (\S+) ,S/hbar (\S+)", output).groups(),
        )
        expected_potential = potential.sum() / (count * 0.2)
        expected_spring = (
            0.5
            * (count * 0.2) ** 2
            * np.sum(np.diff(q, axis=0) ** 2 * np.repeat([2.0, 3.0], 3))
            / (count * 0.2)
        )
        assert a == pytest.approx(expected_potential)
        assert b == pytest.approx(expected_spring)
        assert total == pytest.approx(a + b)
    else:
        # Mirroring the half path doubles the bead count and preserves physical curvatures.
        eig = np.loadtxt("freq.dat")
        assert eig.size == 2 * count * 6
        assert "Instanton rate. Nbeads 8" in output

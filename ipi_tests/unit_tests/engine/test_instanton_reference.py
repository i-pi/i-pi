"""Instanton friction terms and postprocessing on quadratic actions."""

import re
import xml.etree.ElementTree as ET
from types import SimpleNamespace

import numpy as np
import pytest

from ipi.engine.motion.instanton import FrictionMapper
from ipi.utils.tools.instanton_postproc import instanton_compute
from ipi.utils.units import unit_to_internal
from ipi_tests.unit_tests.engine.inputs import harmonic_input


@pytest.mark.parametrize("count", [4, 7])
def test_constant_friction_action_gradient_and_hessian(simulation, count):
    with simulation(harmonic_input(nbeads=count)) as sim:
        s = sim.syslist[0]
        eta = np.diag(np.arange(1.0, 7.0))
        mapper = FrictionMapper(False, eta)
        mapper.bind(
            SimpleNamespace(
                beads=s.beads,
                cell=s.cell,
                forces=s.forces,
                nm=s.nm,
                rp_factor=1.0,
                fix=SimpleNamespace(nactive=2),
                coef=np.ones((count + 1, 1)),
                options={"max_ms": 0, "max_e": 0},
            )
        )
        # Constant spectral density over every normal-mode frequency.
        conversion = unit_to_internal("frequency", "inversecm", 1.0)
        data = np.column_stack(
            [np.linspace(0, 2 * max(mapper.omegak) / conversion, 8), np.ones(8)]
        )
        mapper.set_fric_spec_dens(data, 0)
        q = np.linspace(-0.2, 0.4, count)[:, None] * np.arange(1.0, 7.0)[None, :]
        mapper.dbeads.q[:] = q
        mapper.initialize(q, mapper.dforces)
        energy, gradient = mapper.compute_friction_terms()
        n = np.arange(count)
        transform = np.sqrt(2 / count) * np.cos(
            np.pi * n[:, None] * (n[None, :] + 0.5) / count
        )
        transform[0] /= np.sqrt(2)
        omega = 2 * s.nm.omegan * np.sin(np.pi * n / (2 * count))
        bead_matrix = transform.T @ np.diag(omega) @ transform
        expected_gradient = bead_matrix @ q @ eta
        assert energy == pytest.approx(np.sum(q * expected_gradient) / 2, rel=1e-7)
        np.testing.assert_allclose(gradient, expected_gradient, rtol=1e-7, atol=1e-12)
        hessian = mapper.get_fric_rp_hessian(
            np.zeros((count, 6, 6, 6)), mapper.eta, False
        )
        np.testing.assert_allclose(
            hessian, np.kron(bead_matrix, eta), rtol=1e-7, atol=1e-12
        )
        assert np.linalg.eigvalsh(hessian).min() > -1e-12
        # A rigid shift of the whole path must not change the friction action.
        mapper.dbeads.q[:] = q + 0.7
        shifted_energy, shifted_gradient = mapper.compute_friction_terms()
        assert shifted_energy == pytest.approx(energy, abs=1e-12)
        np.testing.assert_allclose(shifted_gradient, gradient, atol=1e-12)


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

"""Instanton friction terms and postprocessing on quadratic actions."""

import re
import xml.etree.ElementTree as ET
from types import SimpleNamespace

import numpy as np
import pytest

from ipi.engine.motion.instanton import FrictionMapper
from ipi.utils.array_backend import to_numpy, xp
from ipi.utils.tools.instanton_postproc import instanton_compute
from ipi.utils.units import unit_to_internal
from ipi_tests.unit_tests.engine.inputs import harmonic_input, motion_input


def friction_mapper(s, eta):
    """A FrictionMapper with a constant friction tensor, bound to system s."""
    mapper = FrictionMapper(False, eta)
    mapper.bind(
        SimpleNamespace(
            beads=s.beads,
            cell=s.cell,
            forces=s.forces,
            nm=s.nm,
            rp_factor=1.0,
            fix=SimpleNamespace(nactive=2),
            coef=xp.ones((s.beads.nbeads + 1, 1)),
            options={"max_ms": 0, "max_e": 0},
        )
    )
    return mapper


@pytest.mark.parametrize("count", [4, 7])
def test_constant_friction_action_gradient_and_hessian(simulation, count):
    with simulation(harmonic_input(nbeads=count)) as sim:
        s = sim.syslist[0]
        eta = np.diag(np.arange(1.0, 7.0))
        mapper = friction_mapper(s, eta)
        # Constant spectral density over every normal-mode frequency.
        conversion = unit_to_internal("frequency", "inversecm", 1.0)
        omegak = to_numpy(mapper.omegak)
        data = np.column_stack(
            [np.linspace(0, 2 * max(omegak) / conversion, 8), np.ones(8)]
        )
        mapper.set_fric_spec_dens(data, 0)
        q = np.linspace(-0.2, 0.4, count)[:, None] * np.arange(1.0, 7.0)[None, :]
        mapper.dbeads.q_[:] = xp.asarray(q)
        mapper.initialize(q, mapper.dforces)
        energy, gradient = mapper.compute_friction_terms()
        n = np.arange(count)
        transform = np.sqrt(2 / count) * np.cos(
            np.pi * n[:, None] * (n[None, :] + 0.5) / count
        )
        transform[0] /= np.sqrt(2)
        omega = 2 * float(s.nm.omegan) * np.sin(np.pi * n / (2 * count))
        bead_matrix = transform.T @ np.diag(omega) @ transform
        expected_gradient = bead_matrix @ q @ eta
        assert to_numpy(energy) == pytest.approx(
            np.sum(q * expected_gradient) / 2, rel=1e-7
        )
        np.testing.assert_allclose(
            to_numpy(gradient), expected_gradient, rtol=1e-7, atol=1e-12
        )
        hessian = mapper.get_fric_rp_hessian(
            np.zeros((count, 6, 6, 6)), mapper.eta, False
        )
        np.testing.assert_allclose(
            hessian, np.kron(bead_matrix, eta), rtol=1e-7, atol=1e-12
        )
        assert np.linalg.eigvalsh(hessian).min() > -1e-12
        # A rigid shift of the whole path must not change the friction action.
        mapper.dbeads.q_[:] = xp.asarray(q + 0.7)
        shifted_energy, shifted_gradient = mapper.compute_friction_terms()
        assert to_numpy(shifted_energy) == pytest.approx(to_numpy(energy), abs=1e-12)
        np.testing.assert_allclose(
            to_numpy(shifted_gradient), to_numpy(gradient), atol=1e-12
        )


def test_friction_without_laplace_transform_is_frequency_independent(simulation):
    """Without a tabulated Laplace transform of the friction kernel, the friction
    is frequency independent (Litman et al., JCP 156, 194106 (2022), Eq. 20)"""
    with simulation(harmonic_input(nbeads=4)) as sim:
        mapper = friction_mapper(sim.syslist[0], np.eye(6))
        mapper.set_fric_spec_dens(np.ones(0), 0.0)
        np.testing.assert_allclose(mapper.fric_LTwk[:, 0], to_numpy(mapper.omegak))


def test_laplace_transform_is_normalised_at_the_reference_energy(simulation):
    """The tabulated Laplace transform of the friction kernel (frequencies in
    cm^-1) is divided by its value at the energy at which the friction tensor
    was computed (in atomic units)"""
    with simulation(harmonic_input(nbeads=4)) as sim:
        mapper = friction_mapper(sim.syslist[0], np.eye(6))
        conversion = unit_to_internal("frequency", "inversecm", 1.0)
        omegak = to_numpy(mapper.omegak)
        top = 2 * max(omegak)
        freq = np.linspace(0, top, 9)
        data = np.column_stack([freq / conversion, 1 + freq / top])
        energy = 0.3 * top
        mapper.set_fric_spec_dens(data, energy)
        expected = omegak * (1 + omegak / top) / (1 + energy / top)
        np.testing.assert_allclose(mapper.fric_LTwk[:, 0], expected, rtol=1e-10)


def test_friction_energy_is_a_term_of_the_ring_polymer(simulation):
    """The energies of the beads, which the optimizers sum, include the friction
    energy of the ring polymer once"""
    with simulation(harmonic_input(nbeads=4)) as sim:
        mapper = friction_mapper(sim.syslist[0], np.diag(np.arange(1.0, 7.0)))
        mapper.set_fric_spec_dens(np.ones(0), 0.0)
        q = np.linspace(-0.2, 0.4, 4)[:, None] * np.arange(1.0, 7.0)[None, :]
        mapper.dbeads.q_[:] = xp.asarray(q)
        mapper.initialize(q, mapper.dforces)
        friction, _ = mapper.compute_friction_terms()
        energy, _ = mapper.evaluate()
        assert friction > 0
        assert np.sum(to_numpy(energy)) == pytest.approx(
            np.sum(to_numpy(mapper.dforces.pots)) + to_numpy(friction), rel=1e-12
        )


FRICTION_WELL = """
__DRIVER_NAME__ = "friction_well"
__DRIVER_CLASS__ = "FrictionWell"

import json
import numpy as np
from ipi.pes.dummy import Dummy_driver


class FrictionWell(Dummy_driver):
    # double well along x of the first atom, harmonic elsewhere, with a
    # position-dependent friction along x of the first atom only

    def compute_structure(self, cell, pos):
        x = np.asarray(pos, dtype=float).reshape(-1)
        pot = -0.5 * x[0] ** 2 + x[0] ** 4 / 16 + 0.5 * np.sum(x[1:] ** 2)
        g = x.copy()
        g[0] = -x[0] + x[0] ** 3 / 4
        eta = np.zeros((len(x), len(x)))
        eta[0, 0] = 0.3 * (1 + 0.5 * x[0] ** 2)
        extras = json.dumps({"friction": eta.tolist()})
        return pot, -g.reshape(np.shape(pos)), np.zeros((3, 3)), extras
"""


def friction_instanton_input(natoms, sd, update, driver):
    """Rate instanton with friction (position dependent if sd) in the
    FRICTION_WELL potential, with all atoms but the first one fixed."""
    count = 4
    options = ET.Element("instanton", mode="rate")
    ET.SubElement(options, "opt").text = "nichols"
    ET.SubElement(options, "hessian_asr").text = "none"
    ET.SubElement(options, "energy_shift").text = "0"
    ET.SubElement(options, "hessian_init").text = "True"
    ET.SubElement(options, "hessian_update").text = update
    ET.SubElement(options, "hessian_final").text = "True"
    ET.SubElement(options, "friction").text = "True"
    ET.SubElement(options, "frictionSD").text = str(sd)
    eta = np.zeros((3 * natoms, 3 * natoms))
    eta[0, 0] = 0.3
    ET.SubElement(options, "eta", shape=f"{eta.shape}").text = str(eta.ravel().tolist())
    h = np.eye(3 * natoms)
    h[0, 0] = -1.0
    ET.SubElement(options, "hessian", shape=f"{h.shape}").text = str(h.ravel().tolist())
    ET.SubElement(options, "delta").text = "0.3"
    tolerances = ET.SubElement(options, "tolerances")
    for name in ["energy", "force", "position"]:
        ET.SubElement(tolerances, name).text = "1e-7"
    root = motion_input("instanton", options, nbeads=count)
    if natoms > 1:
        ET.SubElement(root.find("system/motion"), "fixatoms").text = str(
            list(range(1, natoms))
        )
    temperature = root.find("system/ensemble/temperature")
    temperature.attrib.clear()
    temperature.text = "0.1"
    root.remove(root.find("ffdirect"))
    ff = ET.SubElement(root, "ffdirect", name="harmonic")
    ET.SubElement(ff, "pes").text = "custom"
    ET.SubElement(ff, "pes_path").text = str(driver)
    beads = root.find("system/beads")
    beads.set("natoms", str(natoms))
    for tag, shape, value in [
        ("q", f"({count},{3*natoms})", [0.0] * (3 * natoms * count)),
        ("p", f"({count},{3*natoms})", [0.0] * (3 * natoms * count)),
        ("m", f"({natoms})", [1.0] * natoms),
        ("names", f"({natoms})", ["H"] * natoms),
    ]:
        beads.find(tag).set("shape", shape)
        beads.find(tag).text = str(value)
    return ET.tostring(root, encoding="unicode")


@pytest.mark.parametrize("sd,update", [(False, "powell"), (True, "powell")])
def test_friction_instanton_with_a_fixed_atom(simulation, tmp_path, sd, update):
    """A fixed atom without friction leaves the friction instanton of the other
    atom unchanged"""
    driver = tmp_path / "friction_well.py"
    driver.write_text(FRICTION_WELL)
    paths = []
    for natoms in [1, 2]:
        with simulation(friction_instanton_input(natoms, sd, update, driver)) as sim:
            s = sim.syslist[0]
            s.motion.exit_on_convergence = False
            for step in range(50):
                s.motion.step(step)
                if s.motion.finished:
                    break
            assert s.motion.finished
            paths.append(to_numpy(s.beads.q)[:, :3].copy())
    assert np.ptp(paths[0][:, 0]) > 0.5
    np.testing.assert_allclose(paths[1], paths[0], atol=1e-7)


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

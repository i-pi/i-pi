"""Checks that the numpy, torch and JAX backends give the same dynamics, and
that an MD step does not synchronize the host with a CUDA device."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

import xml.etree.ElementTree as ET

import numpy as np
import pytest

import ipi.utils.array_backend as ab
from ipi.utils.prng import Random
from ipi_tests.unit_tests.engine.inputs import harmonic_input

try:
    import torch

    cuda_available = torch.cuda.is_available()
except ImportError:
    cuda_available = False

LANGEVIN = "<thermostat mode='langevin'><tau units='femtosecond'>100</tau></thermostat>"
SVR = "<thermostat mode='svr'><tau units='femtosecond'>100</tau></thermostat>"
PILE_L = "<thermostat mode='pile_l'><tau units='femtosecond'>100</tau></thermostat>"
PILE_G = "<thermostat mode='pile_g'><tau units='femtosecond'>100</tau></thermostat>"
GLE = """<thermostat mode='gle'><A shape='(3,3)'>
  [ 1e-3, 2e-3, 1e-4, -2e-3, 3e-3, 1e-4, -1e-4, -1e-4, 2e-3 ]</A></thermostat>"""
FFL = "<thermostat mode='ffl'><tau units='femtosecond'>100</tau></thermostat>"
CL = """<thermostat mode='cl'><tau units='femtosecond'>100</tau>
  <intau units='femtosecond'>200</intau></thermostat>"""
MULTI = f"<thermostat mode='multi'>{LANGEVIN}{SVR}</thermostat>"


def barostat_xml(mode, thermostat=LANGEVIN, extra=""):
    return f"""<barostat mode='{mode}'><tau units='femtosecond'>200</tau>
  {thermostat}{extra}</barostat>"""


ISOTROPIC = barostat_xml("isotropic")
FLEXIBLE = barostat_xml("flexible")
# a second forcefield, for multiple time steps, contraction and bias
SECOND = {
    "ff": "<ffdirect name='ff2'><pes>harmonic</pes>"
    "<parameters>{k1: 5e-4}</parameters></ffdirect>"
}
MTS = {
    **SECOND,
    "forces": "<force forcefield='ff'><mts_weights>[0, 1]</mts_weights></force>"
    "<force forcefield='ff2'><mts_weights>[1, 0]</mts_weights></force>",
    "dynamics": "<nmts>[1, 2]</nmts>",
}

# (number of beads, ensemble, thermostat, barostat[, extra input sections])
CASES = {
    "nve": (1, "nve", "", ""),
    "nvt_langevin": (1, "nvt", LANGEVIN, ""),
    "nvt_svr": (1, "nvt", SVR, ""),
    "nvt_gle": (1, "nvt", GLE, ""),
    "nvt_ffl": (1, "nvt", FFL, ""),
    "nvt_cl": (1, "nvt", CL, ""),
    "nvt_multi": (1, "nvt", MULTI, ""),
    "nvt_mts_bias": (
        1,
        "nvt",
        LANGEVIN,
        "",
        {**MTS, "ensemble": "<bias><force forcefield='ff2'/></bias>"},
    ),
    "npt_isotropic": (1, "npt", LANGEVIN, ISOTROPIC),
    "npt_flexible": (1, "npt", LANGEVIN, FLEXIBLE),
    "npt_svr_barostat": (1, "npt", LANGEVIN, barostat_xml("flexible", SVR)),
    "npt_gle_barostat": (1, "npt", LANGEVIN, barostat_xml("flexible", GLE)),
    "nst_anisotropic": (
        1,
        "nst",
        LANGEVIN,
        barostat_xml(
            "anisotropic", extra="<h0 shape='(3,3)'>[10, 0, 0, 0, 10, 0, 0, 0, 10]</h0>"
        ),
        {
            "ensemble": "<stress shape='(3,3)'>[1e-6, 0, 0, 0, 2e-6, 0, 0, 0, 1e-6]</stress>"
        },
    ),
    "pimd_nvt": (4, "nvt", PILE_L, ""),
    "pimd_pile_g": (4, "nvt", PILE_G, ""),
    "pimd_npt": (4, "npt", PILE_L, ISOTROPIC),
    "pimd_nvt_cc": (4, "nvt-cc", PILE_L, ""),
    "pimd_fixatoms": (
        4,
        "nvt",
        PILE_G,
        "",
        {"motion": "<fixcom>False</fixcom><fixatoms>[0, 3]</fixatoms>"},
    ),
    "pimd_sc": (4, "sc", PILE_L, ""),
    "pimd_sc_forward": (
        4,
        "sc",
        PILE_L,
        "",
        {"forces": "<force forcefield='ff' fd_epsilon='1e-3'/>"},
    ),
    "pimd_scnpt": (4, "scnpt", PILE_L, barostat_xml("sc-isotropic")),
    "pimd_contracted": (
        4,
        "nvt",
        PILE_L,
        "",
        {
            **SECOND,
            "forces": "<force forcefield='ff' nbeads='2'/><force forcefield='ff2'/>",
        },
    ),
}

HARMONIC = (
    "<ffdirect name='ff'><pes>harmonic</pes>"
    "<parameters>{k1: 1e-3}</parameters></ffdirect>",
    "<force forcefield='ff'/>",
)
# a force component with zero weight is never evaluated
NO_FORCES = (
    "<fflj name='ff' pbc='False'><parameters>{eps: 0.0, sigma: 1.0}</parameters></fflj>",
    "<force forcefield='ff' weight='0.0'/>",
)


def write_input(case, forcefield, natoms=8):
    """Writes `init.xyz` in the current folder and returns the XML input."""

    nbeads, mode, thermostat, barostat, *extra = CASES[case]
    extra = extra[0] if extra else {}
    box = (natoms / 0.05) ** (1.0 / 3.0)
    pos = np.random.default_rng(12345).random((natoms, 3)) * box
    with open("init.xyz", "w") as f:
        f.write(f"{natoms}\n# CELL(abcABC): {box} {box} {box} 90 90 90 ")
        f.write("cell{angstrom} positions{angstrom}\n")
        for x, y, z in pos:
            f.write(f"H {x:.8f} {y:.8f} {z:.8f}\n")
    pressure = "<pressure units='bar'>1</pressure>" if barostat else ""
    ffxml, forcexml = forcefield
    forcexml = extra.get("forces", forcexml)
    return f"""
<simulation verbosity='quiet' threading='false' safe_stride='100000000'>
  {ffxml}{extra.get("ff", "")}
  <total_steps>100000000</total_steps>
  <output prefix='simulation'/>
  <prng><seed>12345</seed></prng>
  <system>
    <initialize nbeads='{nbeads}'>
      <file mode='xyz'>init.xyz</file>
      <velocities mode='thermal' units='kelvin'>300</velocities>
    </initialize>
    <forces>{forcexml}</forces>
    <ensemble>
      <temperature units='kelvin'>300</temperature>{pressure}{extra.get("ensemble", "")}
    </ensemble>
    <motion mode='dynamics'>{extra.get("motion", "")}
      <dynamics mode='{mode}'>
        <timestep units='femtosecond'>0.5</timestep>{thermostat}{barostat}
        {extra.get("dynamics", "")}
      </dynamics>
    </motion>
  </system>
</simulation>
"""


def mts_input(forcefield):
    """PIMD with a forcefield and the dummy potential applied at two levels of
    multiple time stepping."""

    xml = write_input("pimd_nvt", forcefield)
    return (
        xml.replace(
            forcefield[0],
            forcefield[0] + "<ffdirect name='ff2'><pes>dummy</pes></ffdirect>",
        )
        .replace(
            "<force forcefield='ff'/>",
            "<force forcefield='ff'><mts_weights>[1,0]</mts_weights></force>"
            "<force forcefield='ff2'><mts_weights>[0,1]</mts_weights></force>",
        )
        .replace(
            "<timestep units='femtosecond'>0.5</timestep>",
            "<timestep units='femtosecond'>0.5</timestep><nmts>[1,3]</nmts>",
        )
    )


def nst_input(forcefield):
    """Classical dynamics at constant stress, with the anisotropic barostat."""

    box = (8 / 0.05) ** (1.0 / 3.0)  # the box of write_input, in angstrom
    return (
        write_input("npt_isotropic", forcefield)
        .replace(
            "<pressure units='bar'>1</pressure>",
            "<stress units='bar'>[1, 0, 0, 0, 1, 0, 0, 0, 1]</stress>",
        )
        .replace("mode='npt'", "mode='nst'")
        .replace(
            "<barostat mode='isotropic'>",
            "<barostat mode='anisotropic'><h0 units='angstrom'>"
            f"[{box}, 0, 0, 0, {box}, 0, 0, 0, {box}]</h0>",
        )
    )


INPUTS = {"mts": mts_input, "nst": nst_input}


@pytest.fixture
def numpy_draws(monkeypatch):
    """Makes the torch and JAX backends draw the same random numbers as
    numpy, and restores the numpy backend afterwards."""

    def gvec(self, shape):
        return ab.xp.asarray(self.gvec_serial(shape))

    def gammavec(self, k, shape):
        return ab.xp.asarray(Random.gammavec(self, k, shape))

    def gfill(self, out):
        out[...] = ab.xp.asarray(self.gvec_serial(tuple(out.shape)))

    monkeypatch.setattr(Random, "gvec_torch", gvec)
    monkeypatch.setattr(Random, "gfill_torch", gfill)
    monkeypatch.setattr(Random, "gvec_jax", gvec)
    monkeypatch.setattr(Random, "gammavec_torch", gammavec)
    monkeypatch.setattr(Random, "gammavec_jax", gammavec)
    yield
    ab.set_backend("numpy")


def trajectory(xml, backend, device, nsteps):
    """Runs nsteps and returns positions, momenta, cell and conserved energy."""

    from ipi.scripting import InteractiveSimulation

    sim = InteractiveSimulation(xml, backend=backend, device=device)
    frames = []
    for step in range(nsteps):
        sim.run_step(step)
        system = sim.syslist[0]
        frames.append(
            [
                np.array(ab.to_numpy(system.beads.q)),
                np.array(ab.to_numpy(system.beads.p)),
                np.array(ab.to_numpy(system.cell.h)),
                float(sim.properties("conserved")),
            ]
        )
    return frames


BACKENDS = [("torch", "cpu"), ("jax", "cpu")] + (
    [("torch", "cuda")] if cuda_available else []
)


@pytest.mark.parametrize("backend,device", BACKENDS)
@pytest.mark.parametrize("case", list(CASES))
def test_backend_equivalence(
    case, backend, device, numpy_draws, request, tmp_path, monkeypatch
):
    """torch and JAX follow the numpy trajectory when they draw the same
    random numbers"""

    pytest.importorskip(backend)
    if backend == "torch":
        request.getfixturevalue("strict_torch")
    if backend == "jax" and case in ("nvt_gle", "npt_gle_barostat"):
        pytest.xfail("the GLE thermostat modifies its auxiliary momenta in place")
    monkeypatch.chdir(tmp_path)
    xml = write_input(case, HARMONIC)
    reference = trajectory(xml, "numpy", "cpu", 10)
    result = trajectory(xml, backend, device, 10)
    for ref, res in zip(reference, result):
        for a, b in zip(ref, res):
            assert np.allclose(a, b, rtol=1e-10, atol=1e-12)


# estimators computed from the positions and the forces of the beads
ESTIMATORS = [
    ("properties", "kinetic_cv"),
    ("properties", "kinetic_cv(1)"),
    ("properties", "kinetic_tens"),
    ("properties", "kinetic_tens(2)"),
    ("properties", "kinetic_ij(1;2)"),
    ("properties", "kinetic_ij(3;3)"),
    ("trajs", "kinetic_cv"),
    ("trajs", "kinetic_od"),
    ("trajs", "r_gyration"),
]


@pytest.mark.parametrize("backend,device", BACKENDS)
@pytest.mark.parametrize("case", ["nvt_langevin", "pimd_nvt"])
def test_estimators_backend_equivalence(
    case, backend, device, numpy_draws, request, tmp_path, monkeypatch
):
    """torch and JAX compute the same kinetic energy estimators and gyration
    radii as numpy"""

    from ipi.scripting import InteractiveSimulation

    pytest.importorskip(backend)
    if backend == "torch":
        request.getfixturevalue("strict_torch")
    monkeypatch.chdir(tmp_path)
    xml = write_input(case, HARMONIC)
    values = []
    for name, dev in [("numpy", "cpu"), (backend, device)]:
        sim = InteractiveSimulation(xml, backend=name, device=dev)
        sim.run_step(0)
        system = sim.syslist[0]
        values.append(
            [
                np.array(ab.to_numpy(getattr(system, kind)[key][0]))
                for kind, key in ESTIMATORS
            ]
        )
    for estimator, a, b in zip(ESTIMATORS, *values):
        assert np.allclose(a, b, rtol=1e-10, atol=1e-12), estimator


@pytest.mark.parametrize("multi", [False, True])
def test_jax_refuses_unsupported_motion(multi, tmp_path, monkeypatch):
    """A motion that does not work with JAX, also within a multi motion, stops
    with an error when the simulation is set up"""

    from ipi.scripting import InteractiveSimulation

    pytest.importorskip("jax")
    monkeypatch.chdir(tmp_path)
    root = ET.fromstring(harmonic_input())
    system = root.find("system")
    dynamics = system.find("motion")
    system.remove(dynamics)
    motion = ET.SubElement(system, "motion", mode="minimize")
    ET.SubElement(motion, "optimizer", mode="sd")
    if multi:
        system.remove(motion)
        ET.SubElement(system, "motion", mode="multi").extend([dynamics, motion])
    try:
        with pytest.raises(ValueError, match="GeopMotion cannot run with the JAX"):
            InteractiveSimulation(
                ET.tostring(root, encoding="unicode"), backend="jax", device="cpu"
            )
    finally:
        ab.set_backend("numpy")


@pytest.mark.skipif(not cuda_available, reason="CUDA not available")
@pytest.mark.parametrize(
    "case",
    [
        "nve",
        "nvt_langevin",
        "nvt_svr",
        "pimd_nvt",
        "pimd_pile_g",
        "pimd_npt",
        "npt_flexible",
    ],
)
def test_cuda_step_is_sync_free(case, tmp_path, monkeypatch, device_random):
    """An MD step without force evaluation does not wait for the device"""

    from ipi.scripting import InteractiveSimulation

    monkeypatch.chdir(tmp_path)
    sim = InteractiveSimulation(
        write_input(case, NO_FORCES), backend="torch", device="cuda"
    )
    try:
        # the first steps compute and cache the integrator constants
        for step in range(2):
            sim.run_step(step)
        torch.cuda.set_sync_debug_mode("error")
        for step in range(2, 5):
            sim.run_step(step)
    finally:
        torch.cuda.set_sync_debug_mode("default")
        ab.set_backend("numpy")

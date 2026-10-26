"""Checks that the numpy, torch and JAX backends give the same dynamics, and
that an MD step does not synchronize the host with a CUDA device."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

import numpy as np
import pytest

import ipi.utils.array_backend as ab
from ipi.utils.prng import Random

try:
    import torch

    cuda_available = torch.cuda.is_available()
except ImportError:
    cuda_available = False

LANGEVIN = "<thermostat mode='langevin'><tau units='femtosecond'>100</tau></thermostat>"
SVR = "<thermostat mode='svr'><tau units='femtosecond'>100</tau></thermostat>"
PILE_L = "<thermostat mode='pile_l'><tau units='femtosecond'>100</tau></thermostat>"
GLE = """<thermostat mode='gle'><A shape='(3,3)'>
  [ 1e-3, 2e-3, 1e-4, -2e-3, 3e-3, 1e-4, -1e-4, -1e-4, 2e-3 ]</A></thermostat>"""
ISOTROPIC = f"<barostat mode='isotropic'><tau units='femtosecond'>200</tau>{LANGEVIN}</barostat>"
FLEXIBLE = (
    f"<barostat mode='flexible'><tau units='femtosecond'>200</tau>{LANGEVIN}</barostat>"
)

# (number of beads, ensemble, thermostat, barostat)
CASES = {
    "nve": (1, "nve", "", ""),
    "nvt_langevin": (1, "nvt", LANGEVIN, ""),
    "nvt_svr": (1, "nvt", SVR, ""),
    "nvt_gle": (1, "nvt", GLE, ""),
    "npt_isotropic": (1, "npt", LANGEVIN, ISOTROPIC),
    "npt_flexible": (1, "npt", LANGEVIN, FLEXIBLE),
    "pimd_nvt": (4, "nvt", PILE_L, ""),
    "pimd_npt": (4, "npt", PILE_L, ISOTROPIC),
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

    nbeads, mode, thermostat, barostat = CASES[case]
    box = (natoms / 0.05) ** (1.0 / 3.0)
    pos = np.random.default_rng(12345).random((natoms, 3)) * box
    with open("init.xyz", "w") as f:
        f.write(f"{natoms}\n# CELL(abcABC): {box} {box} {box} 90 90 90 ")
        f.write("cell{angstrom} positions{angstrom}\n")
        for x, y, z in pos:
            f.write(f"H {x:.8f} {y:.8f} {z:.8f}\n")
    pressure = "<pressure units='bar'>1</pressure>" if barostat else ""
    ffxml, forcexml = forcefield
    return f"""
<simulation verbosity='quiet' threading='false' safe_stride='100000000'>
  {ffxml}
  <total_steps>100000000</total_steps>
  <output prefix='simulation'/>
  <prng><seed>12345</seed></prng>
  <system>
    <initialize nbeads='{nbeads}'>
      <file mode='xyz'>init.xyz</file>
      <velocities mode='thermal' units='kelvin'>300</velocities>
    </initialize>
    <forces>{forcexml}</forces>
    <ensemble><temperature units='kelvin'>300</temperature>{pressure}</ensemble>
    <motion mode='dynamics'>
      <dynamics mode='{mode}'>
        <timestep units='femtosecond'>0.5</timestep>{thermostat}{barostat}
      </dynamics>
    </motion>
  </system>
</simulation>
"""


@pytest.fixture
def numpy_draws(monkeypatch):
    """Makes the torch and JAX backends draw the same Gaussian numbers as
    numpy, and restores the numpy backend afterwards."""

    def gvec(self, shape):
        return ab.xp.asarray(self.gvec_serial(shape))

    def gfill(self, out):
        out[...] = ab.xp.asarray(self.gvec_serial(tuple(out.shape)))

    monkeypatch.setattr(Random, "gvec_torch", gvec)
    monkeypatch.setattr(Random, "gfill_torch", gfill)
    monkeypatch.setattr(Random, "gvec_jax", gvec)
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
    if backend == "jax" and case == "nvt_gle":
        pytest.xfail("the GLE thermostat modifies its auxiliary momenta in place")
    monkeypatch.chdir(tmp_path)
    xml = write_input(case, HARMONIC)
    reference = trajectory(xml, "numpy", "cpu", 10)
    result = trajectory(xml, backend, device, 10)
    for ref, res in zip(reference, result):
        for a, b in zip(ref, res):
            assert np.allclose(a, b, rtol=1e-10, atol=1e-12)


@pytest.mark.skipif(not cuda_available, reason="CUDA not available")
@pytest.mark.parametrize("case", ["nve", "nvt_langevin", "pimd_nvt", "pimd_npt"])
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

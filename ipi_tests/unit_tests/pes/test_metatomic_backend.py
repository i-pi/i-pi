"""Checks that metatomic models evaluated through a direct forcefield work on
the arrays of the torch backend without converting them to numpy, and give the
same results as with numpy arrays."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

import shutil
from pathlib import Path

import numpy as np
import pytest

import ipi.engine.forcefields as forcefields
import ipi.utils.array_backend as ab
from ipi.utils.units import unit_to_internal

torch = pytest.importorskip("torch")
pytest.importorskip("metatomic.torch")
pytest.importorskip("vesin.metatomic")
ase_io = pytest.importorskip("ase.io")

# Lennard-Jones model for nickel, in float64 and without extensions
LJ = Path(__file__).resolve().parents[3] / "examples" / "clients" / "metatomic" / "lj"
devices = ["cpu"] + (["cuda"] if torch.cuda.is_available() else [])


@pytest.fixture
def numpy_backend():
    """Restores the numpy backend after the test."""
    yield
    ab.set_backend("numpy")


def nickel():
    """Returns the cell (lattice vectors as columns) and the positions of the
    nickel structure, in atomic units."""

    atoms = ase_io.read(LJ / "nickel.xyz")
    cell = unit_to_internal("length", "angstrom", atoms.cell.array.T.copy())
    pos = unit_to_internal("length", "angstrom", atoms.positions.copy())
    return cell, pos


@pytest.mark.parametrize("device", devices)
def test_driver_accepts_tensors(device):
    """Tensors give the same results as numpy arrays, as tensors on their device"""

    from ipi.pes.metatomic import MetatomicDriver

    driver = MetatomicDriver(
        template=str(LJ / "nickel.xyz"), model=str(LJ / "nickel-lj.pt"), device=device
    )
    cell, pos = nickel()
    cells = [cell, cell]
    positions = [pos, pos * 1.01]

    reference = [driver(cell, pos)] + driver(cells, positions)
    tensor = lambda x: torch.as_tensor(x, device=device)
    results = [driver(tensor(cell), tensor(pos))] + driver(
        [tensor(c) for c in cells], [tensor(p) for p in positions]
    )
    for ref, res in zip(reference, results):
        for a, b in zip(ref[:3], res[:3]):
            assert isinstance(b, torch.Tensor)
            assert b.device.type == device and b.dtype == torch.float64
            assert np.allclose(np.asarray(a), b.cpu().numpy(), rtol=1e-12, atol=1e-14)
        assert ref[3] == res[3]


def write_input(nbeads, device, rotations=False):
    """Writes the nickel structure in the current folder, and returns the
    input of a PIMD simulation with the LJ model in a direct forcefield."""

    shutil.copy(LJ / "nickel.xyz", "init.xyz")
    model = LJ / "nickel-lj.pt"
    ff = (
        "<ffdirect name='lj' pbc='false'><pes>metatomic</pes>"
        f"<parameters>{{template:init.xyz,model:{model},device:{device}}}</parameters>"
        f"<batch_size>{2 if rotations else nbeads}</batch_size></ffdirect>"
    )
    if rotations:
        ff = (
            f"<ffrotations name='ff' pbc='false'>{ff}<random>True</random>"
            "<inversion>True</inversion><prng><seed>12</seed></prng></ffrotations>"
        )
    else:
        ff = ff.replace("name='lj'", "name='ff'")
    return f"""
<simulation verbosity='quiet' threading='false'>
  {ff}
  <total_steps>10</total_steps>
  <output prefix='lj'/>
  <prng><seed>12345</seed></prng>
  <system>
    <initialize nbeads='{nbeads}'>
      <file mode='ase' units='angstrom'>init.xyz</file>
      <velocities mode='thermal' units='kelvin'>300</velocities>
    </initialize>
    <forces><force forcefield='ff'/></forces>
    <ensemble><temperature units='kelvin'>300</temperature></ensemble>
    <motion mode='dynamics'>
      <dynamics mode='nve'><timestep units='femtosecond'>1.0</timestep></dynamics>
    </motion>
  </system>
</simulation>
"""


def forces(xml, backend, device):
    """Returns forces, potentials and virial of the initial configuration."""

    from ipi.scripting import InteractiveSimulation

    sim = InteractiveSimulation(xml, backend=backend, device=device)
    try:
        f = sim.syslist[0].forces
        return [np.array(ab.to_numpy(x)) for x in (f.f, f.pots, f.vir)]
    finally:
        sim.stop()


@pytest.mark.parametrize("device", devices)
@pytest.mark.parametrize("nbeads", [1, 4])
def test_forces_stay_on_device(
    nbeads, device, numpy_backend, strict_torch, tmp_path, monkeypatch
):
    """With the torch backend, positions and forces are not converted to numpy"""

    monkeypatch.chdir(tmp_path)
    xml = write_input(nbeads, device)
    reference = forces(xml, "numpy", "cpu")

    def no_conversion(x):
        raise AssertionError("the request was converted to numpy")

    monkeypatch.setattr(forcefields, "to_numpy", no_conversion)
    results = forces(xml, "torch", device)
    for a, b in zip(reference, results):
        assert np.allclose(a, b, rtol=1e-12, atol=1e-14)


@pytest.mark.parametrize("device", devices)
def test_rotations_around_backend_forcefield(
    device, numpy_backend, strict_torch, tmp_path, monkeypatch
):
    """The rotational average works on the results of a forcefield that works
    on the arrays of the backend, and leaves the LJ forces unchanged"""

    monkeypatch.chdir(tmp_path)
    reference = forces(write_input(2, device), "torch", device)
    results = forces(write_input(2, device, rotations=True), "torch", device)
    for a, b in zip(reference, results):
        assert np.allclose(a, b, rtol=1e-10, atol=1e-12)

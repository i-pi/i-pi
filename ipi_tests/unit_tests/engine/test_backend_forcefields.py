"""Forcefields whose drivers work on the arrays of the array backend
(`backend_arrays`) give the same results with torch tensors as with numpy
arrays, alone, in a committee and in a rotational average."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

import numpy as np
import pytest

import ipi.utils.array_backend as ab
from ipi.pes.harmonic import Harmonic_driver
from ipi_tests.unit_tests.engine.inputs import write_input

torch = pytest.importorskip("torch")
devices = ["cpu"] + (["cuda"] if torch.cuda.is_available() else [])


def harmonic(name, k):
    return (
        f"<ffdirect name='{name}' pbc='false'><pes>harmonic</pes>"
        f"<parameters>{{k1: {k}}}</parameters></ffdirect>"
    )


# the isotropic harmonic potential does not change under rotations
FORCEFIELDS = {
    "plain": harmonic("ff", 1e-3),
    "committee": "<ffcommittee name='ff'>"
    + harmonic("a", 1e-3)
    + harmonic("b", 2e-3)
    + "</ffcommittee>",
    "rotations": "<ffrotations name='ff' pbc='false'>"
    + harmonic("h", 1e-3)
    + "<random>True</random><inversion>True</inversion>"
    + "<prng><seed>12</seed></prng></ffrotations>",
}


@pytest.fixture
def received(monkeypatch):
    """Makes the harmonic driver work on the arrays of the backend, and returns
    the types of the positions it receives."""

    types = []
    compute_structure = Harmonic_driver.compute_structure

    def record(self, cell, pos):
        types.append(type(pos))
        return compute_structure(self, cell, pos)

    monkeypatch.setattr(Harmonic_driver, "backend_arrays", True, raising=False)
    monkeypatch.setattr(Harmonic_driver, "compute_structure", record)
    return types


def forces(simulation, xml, backend, device):
    """Returns forces, potentials and virial of the initial configuration."""

    with simulation(xml, backend=backend, device=device) as sim:
        f = sim.syslist[0].forces
        return [np.array(ab.to_numpy(x)) for x in (f.f, f.pots, f.vir)]


@pytest.mark.parametrize("device", devices)
@pytest.mark.parametrize("kind", list(FORCEFIELDS))
def test_backend_forcefields(kind, device, simulation, received, strict_torch):
    """The drivers receive tensors, and give the forces computed with numpy"""

    xml = write_input("pimd_nvt", (FORCEFIELDS[kind], "<force forcefield='ff'/>"))
    reference = forces(simulation, xml, "numpy", "cpu")
    received.clear()
    results = forces(simulation, xml, "torch", device)
    assert received and all(t is torch.Tensor for t in received)
    for a, b in zip(reference, results):
        assert np.allclose(a, b, rtol=1e-10, atol=1e-12)

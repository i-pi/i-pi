"""Tests for molecular dynamics with holonomic constraints."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

import numpy as np
import pytest

import ipi.utils.array_backend as ab
from ipi.utils.prng import Random
from ipi.utils.units import Elements, unit_to_internal

try:
    import torch

    devices = ["cpu"] + (["cuda"] if torch.cuda.is_available() else [])
except ImportError:
    devices = []

# four water molecules: bond and angle constraints with values taken from the
# initial geometry, bond and angle constraints with explicit values,
# bond constraints only, and Eckart constraints
XML = """
<simulation verbosity='quiet' threading='false' safe_stride='100000000'>
  <ffdirect name='ff'><pes>harmonic</pes><parameters>{k1: 1e-2}</parameters></ffdirect>
  <total_steps>1000</total_steps>
  <output prefix='simulation'/>
  <prng><seed>31415</seed></prng>
  <system>
    <initialize nbeads='1'>
      <file mode='xyz'>init.xyz</file>
      <velocities mode='thermal' units='kelvin'>300</velocities>
    </initialize>
    <forces><force forcefield='ff'/></forces>
    <ensemble><temperature units='kelvin'>300</temperature></ensemble>
    <motion mode='constrained_dynamics'>
      <fixcom>false</fixcom>
      <constrained_dynamics mode='nvt' splitting='baoab'>
        <nsteps_geo>3</nsteps_geo>
        <nsteps_o>2</nsteps_o>
        <csolver><maxit>1000</maxit><norm_order>2</norm_order><tolerance>1e-10</tolerance></csolver>
        <timestep units='femtosecond'>1.0</timestep>
        <thermostat mode='langevin'><tau units='femtosecond'>10</tau></thermostat>
        <constraint mode='multi'>
          <constraint mode='angle'><atoms>[0,1,2]</atoms></constraint>
          <constraint mode='distance'><atoms>[0,1,0,2]</atoms></constraint>
        </constraint>
        <constraint mode='multi'>
          <constraint mode='angle'><atoms>[3,4,5]</atoms><values units='degree'>[104.0]</values></constraint>
          <constraint mode='distance'><atoms>[3,4,3,5]</atoms><values units='angstrom'>[0.97,0.95]</values></constraint>
        </constraint>
        <constraint mode='distance'><atoms>[6,7,6,8]</atoms></constraint>
        <constraint mode='eckart'><atoms>[9,10,11]</atoms></constraint>
      </constrained_dynamics>
    </motion>
  </system>
</simulation>
"""

WATER = np.array([[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]])


def write_init():
    """Writes `init.xyz` in the current folder and returns the positions in bohr."""

    pos = np.concatenate(
        [WATER + [3.0 * m + 2.0, 2.0 + 0.1 * m, 2.0] for m in range(4)]
    )
    with open("init.xyz", "w") as f:
        f.write(
            "12\n# CELL(abcABC): 12 12 12 90 90 90 cell{angstrom} positions{angstrom}\n"
        )
        for name, x in zip("OHH" * 4, pos):
            f.write(f"{name} {x[0]:.6f} {x[1]:.6f} {x[2]:.6f}\n")
    return unit_to_internal("length", "angstrom", 1.0) * pos


def run(backend, device, nsteps=10):
    """Runs nsteps of constrained dynamics, returning the positions (as an
    array of shape (natoms, 3)) and the conserved quantity at each step."""

    from ipi.scripting import InteractiveSimulation

    sim = InteractiveSimulation(XML, backend=backend, device=device)
    traj = []
    for step in range(nsteps):
        sim.run_step(step)
        q = np.array(ab.to_numpy(sim.syslist[0].beads.q)).reshape(-1, 3)
        traj.append((q, float(sim.properties("conserved"))))
    return traj


def angle(q, i, j, k):
    a, b = q[j] - q[i], q[k] - q[i]
    return np.arccos(a @ b / np.linalg.norm(a) / np.linalg.norm(b))


@pytest.fixture
def numpy_backend(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    yield
    ab.set_backend("numpy")


def test_constraints_are_enforced(numpy_backend):
    """Bond lengths and angles keep their target values"""

    q0 = write_init()
    angstrom = unit_to_internal("length", "angstrom", 1.0)
    for q, _ in run("numpy", "cpu"):
        # targets taken from the initial geometry
        for i, j in [(0, 1), (0, 2), (6, 7), (6, 8)]:
            assert np.isclose(
                np.linalg.norm(q[j] - q[i]), np.linalg.norm(q0[j] - q0[i]), rtol=1e-8
            )
        assert np.isclose(angle(q, 0, 1, 2), angle(q0, 0, 1, 2), rtol=1e-8)
        # explicit targets
        assert np.isclose(np.linalg.norm(q[4] - q[3]), 0.97 * angstrom, rtol=1e-8)
        assert np.isclose(np.linalg.norm(q[5] - q[3]), 0.95 * angstrom, rtol=1e-8)
        assert np.isclose(angle(q, 3, 4, 5), np.radians(104.0), rtol=1e-8)
        # Eckart: the centre of mass of the last molecule does not move
        m = np.array([Elements.mass("O"), Elements.mass("H"), Elements.mass("H")])
        assert np.allclose(m @ q[9:12] / m.sum(), m @ q0[9:12] / m.sum(), atol=1e-8)


@pytest.mark.parametrize("device", devices)
def test_constrained_backend_equivalence(
    device, numpy_backend, strict_torch, monkeypatch
):
    """torch follows the numpy trajectory when it draws the same random numbers"""

    def gvec(self, shape):
        return ab.xp.asarray(self.gvec_serial(shape))

    monkeypatch.setattr(Random, "gvec_torch", gvec)
    write_init()
    reference = run("numpy", "cpu")
    for (qr, cr), (q, c) in zip(reference, run("torch", device)):
        assert np.allclose(q, qr, rtol=1e-10, atol=1e-12)
        assert np.isclose(c, cr, rtol=1e-10)

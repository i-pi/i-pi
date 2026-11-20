"""Systems stepped in tasks, which switch while they wait for forces, give
the same results as systems stepped in turn, also when they make different
numbers of force requests, and when their outputs request forces."""

import glob

import numpy as np
import pytest

from ipi.engine.forcefields import FFDirect
from ipi.engine.simulation import Simulation
from ipi.utils.softexit import softexit

NSYSTEMS = 3
NSTEPS = 4


def systems_xml(interleave, motion, batch_size=1, nbeads=2, outputs=""):
    """A simulation with NSYSTEMS two-atom systems, at different positions,
    sharing a harmonic forcefield."""

    systems = []
    for i in range(NSYSTEMS):
        q = np.tile([0.2, -0.3, 0.1, -0.4, 0.2, 0.3], (nbeads, 1)) * (1 + 0.1 * i)
        q[:, 0] += np.arange(nbeads) * 0.01
        p = np.tile([0.1, 0.2, -0.1, 0.3, -0.2, 0.1], (nbeads, 1))
        systems.append(f"""<system prefix="s{i}">
          <forces><force forcefield="harmonic"/></forces>
          <ensemble><temperature units="kelvin">300</temperature></ensemble>
          {motion}
          <beads natoms="2" nbeads="{nbeads}">
            <q shape="({nbeads},6)">{q.ravel().tolist()}</q>
            <p shape="({nbeads},6)">{p.ravel().tolist()}</p>
            <m shape="(2)">[2, 3]</m><names shape="(2)">[H, O]</names>
          </beads><cell shape="(3,3)">[20,0,0,0,20,0,0,0,20]</cell>
        </system>""")
    return f"""<simulation verbosity="quiet" interleave="{interleave}" floatformat="%24.16e">
      <output prefix="sim">{outputs}</output>
      <total_steps>{NSTEPS}</total_steps><prng><seed>317</seed></prng>
      <ffdirect name="harmonic"><pes>harmonic</pes><parameters>{{k1: 0.5}}</parameters>
        <batch_size>{batch_size}</batch_size></ffdirect>
      {"".join(systems)}
    </simulation>"""


NVE = """<motion mode="dynamics"><fixcom>False</fixcom>
  <dynamics mode="nve"><timestep>0.02</timestep></dynamics></motion>"""

SWAPS = """<motion mode="multi">
  <motion mode="dynamics"><fixcom>False</fixcom>
    <dynamics mode="nvt"><timestep>0.02</timestep>
      <thermostat mode="svr"><tau>0.5</tau></thermostat></dynamics></motion>
  <motion mode="atomswap"><atomswap><nxc>2</nxc><names>[H, O]</names></atomswap></motion>
</motion>"""

ESTIMATORS = """<properties stride="1" filename="out">
  [ step, conserved, potential, kinetic_cv, isotope_scfep(1.5), isotope_tdfep(1.5) ]
</properties>
<trajectory stride="1" filename="pos">positions</trajectory>"""

TRAJECTORY = """<properties stride="1" filename="out">[ step, potential ]</properties>
<trajectory stride="1" filename="pos">positions</trajectory>"""


def run(tmp_path, name, xml, monkeypatch):
    """Runs a simulation in its own folder, and returns its outputs."""

    folder = tmp_path / name
    folder.mkdir()
    monkeypatch.chdir(folder)
    softexit.reset()
    simulation = Simulation.load_from_xml(xml)
    simulation.run()
    softexit.cleanup()
    softexit.reset()
    outputs = {}
    for filename in sorted(glob.glob("*sim.*")):
        with open(filename) as output:
            outputs[filename] = output.read()
    return outputs


def assert_same_outputs(first, second):
    assert first.keys() == second.keys()
    for name in first:
        assert first[name] == second[name], name


def count_driver_calls(monkeypatch):
    """Counts the calls to the drivers of direct forcefields, and the number of
    structures in each."""

    calls = []
    call_driver = FFDirect._call_driver

    def counting(self, cell, pos):
        calls.append(len(pos) if isinstance(cell, list) else 1)
        return call_driver(self, cell, pos)

    monkeypatch.setattr(FFDirect, "_call_driver", counting)
    return calls


@pytest.mark.parametrize("batch_size", [1, NSYSTEMS * 2])
def test_systems_in_tasks_match_systems_in_turn(tmp_path, monkeypatch, batch_size):
    """With outputs that compute forces (the isotope estimators displace each
    atom), and no random numbers, the tasks give the same results as the
    systems stepped one after the other."""

    reference = run(
        tmp_path, "turn", systems_xml("false", NVE, outputs=ESTIMATORS), monkeypatch
    )
    calls = count_driver_calls(monkeypatch)
    tasks = run(
        tmp_path,
        "tasks",
        systems_xml("true", NVE, batch_size=batch_size, outputs=ESTIMATORS),
        monkeypatch,
    )
    # the properties, and the positions of the two beads, of each system
    assert len(reference) == 3 * NSYSTEMS
    assert_same_outputs(reference, tasks)
    if batch_size > 1:
        # the beads of all the systems are evaluated together at each step
        assert calls.count(NSYSTEMS * 2) >= NSTEPS


@pytest.mark.parametrize("batch_size", [1, NSYSTEMS])
def test_systems_with_different_numbers_of_requests(tmp_path, monkeypatch, batch_size):
    """Each system tries a random number of atom swaps per step, each with
    its own force evaluation: the runs complete, and are reproducible."""

    xml = systems_xml(
        "true", SWAPS, batch_size=batch_size, nbeads=1, outputs=TRAJECTORY
    )
    calls = count_driver_calls(monkeypatch)
    first = run(tmp_path, "first", xml, monkeypatch)
    nstructures = sum(calls)
    second = run(tmp_path, "second", xml, monkeypatch)
    assert len(first) == 2 * NSYSTEMS
    assert_same_outputs(first, second)
    # more evaluations than the dynamics alone makes
    assert nstructures > NSYSTEMS * (NSTEPS + 1)

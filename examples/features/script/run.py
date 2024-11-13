#!/usr/bin/python

from ipi.engine.simulation import Simulation
from ipi.utils.messages import verbosity
import subprocess

import numpy as np
from ipi.utils.units import Elements
from ipi.inputs.beads import InputBeads
from ipi.inputs.cell import InputCell
from ipi.engine.beads import Beads
from ipi.engine.cell import Cell
from ipi.utils.io.inputs.io_xml import (
    xml_node,
    xml_parse_string,
    xml_parse_file,
    xml_write,
    xml_edit,
)

import ase, ase.io


def setup_system(xml_template, atoms):
    pass


data = ase.io.read("./water-64.xyz", 0)
beads = Beads(nbeads=1, natoms=len(data))
beads.q = data.positions.flatten() / 0.529177
beads.names = np.array(data.symbols)
beads.m = np.array([Elements.mass(n) for n in beads.names])

input_beads = InputBeads()
input_beads.store(beads)
xml_beads = xml_parse_string(input_beads.write("beads")).fields[0][1]

print("CELL ", data.cell)
cell = Cell(h=np.array(data.cell) / 0.529177)
input_cell = InputCell()
input_cell.store(cell)
xml_cell = xml_parse_string(input_cell.write("cell")).fields[0][1]

xml_template = xml_parse_file("input.xml")

print(xml_template)

xml_edit(xml_template, ["system", "initialize"], xml_beads, mode="remove")
xml_edit(xml_template, ["system"], xml_beads, mode="add")
xml_edit(xml_template, ["system"], xml_cell, mode="add")
xml_edited = xml_write(xml_template)
print(xml_edited)
sim = Simulation.load_from_xml(xml_edited)
verbosity.level = "quiet"

print(sim.syslist[0].beads.q[0, :3])
print(sim.syslist[0].beads.p[0, :3])
driver = subprocess.Popen(["i-pi-driver", "-u", "-a", "h2o-md.1", "-m", "qtip4pf"])

print("# istep V K")
for istep in range(10):
    sim.run_step(istep)
    print(
        istep,
        sim.syslist[0].properties["potential"][0],
        sim.syslist[0].properties["kinetic_md"][0],
    )
sim.stop()
driver.wait()

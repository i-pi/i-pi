#!/usr/bin/python

from ipi.engine.simulation import Simulation
from ipi.utils.messages import verbosity
import subprocess

import numpy as np
from ipi.inputs.beads import InputBeads
import ase, ase.io
def setup_system(xml_template, atoms):
    pass

data = ase.io.read("./water_64.xyz", 0)
ibeads = InputBeads()
q = data.positions.flatten()
p = np.zeros_like(q)
m = np.zeros(len(data))
ibeads.q.store(q)
ibeads.names.store(np.array(data.symbols))
ibeads.nbeads.store(1)
ibeads.natoms.store(len(data))
ibeads.p.store(p)
ibeads.m.store(m)
print(ibeads.write("beads"))
print(ibeads.fetch())

sim = Simulation.load_from_xml(open("input.xml"))
verbosity.level = "quiet"

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

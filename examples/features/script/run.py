#!/usr/bin/python

from ipi.engine.simulation import Simulation
from ipi.utils.messages import verbosity
import subprocess

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

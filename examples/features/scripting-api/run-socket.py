#!/usr/bin/python
from ipi.scripting import (
    simulation_xml,
    forcefield_xml,
    motion_nvt_xml,
    InteractiveSimulation,
)
import subprocess
import ase, ase.io

# Forces can also be computed by an external driver that connects to i-PI
# through a socket. Here we use a UNIX domain socket and the Python driver
# that ships with i-PI, but any i-PI client (LAMMPS, CP2K, ...) can be used
data = ase.io.read("water-64.xyz")
input_xml = simulation_xml(
    structures=data,
    forcefield=forcefield_xml(name="harmonic", mode="unix", address="scripting"),
    motion=motion_nvt_xml(timestep=0.5 * ase.units.fs, fixcom=True),
    temperature=400,
    prefix="socket",
)

print("Running with XML input:\n\n", input_xml)

# Creating the InteractiveSimulation opens the socket...
sim = InteractiveSimulation(input_xml)

# ... so the driver must be launched after this point, and before calling `run`,
# which waits for the driver to connect before computing the forces
driver = subprocess.Popen(
    ["i-pi-py_driver", "-u", "-a", "scripting", "-m", "harmonic", "-o", "1e-3"]
)

sim.run(10)
print("# After 10 steps")
print(
    f'Time:        {sim.properties("time") / ase.units.fs}  '
    f'Potential:   {sim.properties("potential")}  '
    f'Temperature: {sim.properties("temperature")}'
)

# the driver stays connected between calls to `run`
sim.run(10)
print("# After 20 steps")
print(
    f'Time:        {sim.properties("time") / ase.units.fs}  '
    f'Potential:   {sim.properties("potential")}  '
    f'Temperature: {sim.properties("temperature")}'
)

# `stop` closes the socket, sending an exit message to the driver.
# Without it, the driver would keep waiting for i-PI after the script ends
sim.stop()
driver.wait()

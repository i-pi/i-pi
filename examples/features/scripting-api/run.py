#!/usr/bin/python
from ipi.utils.scripting import (
    simulation_xml,
    forcefield_xml,
    motion_nvt_xml,
    InteractiveSimulation,
)
from ipi.utils.depend import dstrip
import ase, ase.io

# There are utilities to quickly set up XML inputs for commonly-used simulations
data = ase.io.read("water-64.xyz")
input_xml = simulation_xml(
    structures=data,
    forcefield=forcefield_xml(
        name="harmonic", mode="direct", pes="harmonic", parameters="{k1:1e-3}"
    ),
    motion=motion_nvt_xml(timestep=0.5 * ase.units.fs),
    temperature=400,
    prefix="script",
)

print("Running with XML input:\n\n", input_xml)

# The main object for scripting is `InteractiveSimulation`, that is initialized from
# and XML input and acts as a wrapper around an i-PI simulation object
sim = InteractiveSimulation(input_xml)

# `properties` accesses the (well) properties of the simulation object
print("# Initial properties")
print(
    f'Time:        {sim.properties("time") / ase.units.fs}  '
    f'Potential:   {sim.properties("potential")}  '
    f'Temperature: {sim.properties("temperature")}'
)
# `run` advances the interactive simulation by one (or the prescribed number) of steps
sim.run(10)
print("# After 10 steps")
print(
    f'Time:        {sim.properties("time") / ase.units.fs}  '
    f'Potential:   {sim.properties("potential")}  '
    f'Temperature: {sim.properties("temperature")}'
)
sim.run(10, write_outputs=False)  # we can suppress the outputs
print("# Without outputs")
print(
    f'Time:        {sim.properties("time") / ase.units.fs}  '
    f'Potential:   {sim.properties("potential")}  '
    f'Temperature: {sim.properties("temperature")}'
)
sim.run(10)
print("# With outputs")
print(
    f'Time:        {sim.properties("time") / ase.units.fs}  '
    f'Potential:   {sim.properties("potential")}  '
    f'Temperature: {sim.properties("temperature")}'
)

# `get_structures` dumps the state of the system as ASE Atoms objects, possibly listing
# all systems and beads
ase.io.write("final_positions.xyz", sim.get_structures())

# we can also set the simulation state
structure = sim.get_structures()
structure.positions[:] = data.positions
structure.arrays["momenta"][:] = 0
sim.set_structures(structure)
print("# Resetting position & momenta")
print(
    f'Time:        {sim.properties("time") / ase.units.fs}  '
    f'Potential:   {sim.properties("potential")}  '
    f'Temperature: {sim.properties("temperature")}'
)


# one can do more aggressive "interventions" by monkey-patching the
# InputSimulation object (although this requires an understanding
# of the i-PI internals)
def vv_obabo(self, step=None):
    self.thermostat.step()
    self.beads.p[:] += dstrip(self.forces.f) * self.dt * 0.5
    self.beads.q[:] += dstrip(self.beads.p) / dstrip(self.beads.m3) * self.dt
    self.beads.p[:] += dstrip(self.forces.f) * self.dt * 0.5
    self.thermostat.step()
    self.ensemble.time += self.dt


sim.set_motion_step(vv_obabo)
sim.run(10, write_outputs=True)
print("# Custom step")
print(
    f'Time:        {sim.properties("time") / ase.units.fs}  '
    f'Potential:   {sim.properties("potential")}  '
    f'Temperature: {sim.properties("temperature")}'
)

#!/usr/bin/python
from ipi.utils.scripting import (
    simulation_xml,
    forcefield_xml,
    motion_nvt_xml,
    InteractiveSimulation,
)
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
print(
    sim.properties("time") / ase.units.fs,
    sim.properties("potential"),
    sim.properties("temperature"),
)
# `run` advances the interactive simulation by one (or the prescribed number) of steps
sim.run(10)
print(
    sim.properties("time") / ase.units.fs,
    sim.properties("potential"),
    sim.properties("temperature"),
)
sim.run(10, write_outputs=False)  # we can suppress the outputs
print(
    sim.properties("time") / ase.units.fs,
    sim.properties("potential"),
    sim.properties("temperature"),
)
sim.run(10)
print(
    sim.properties("time") / ase.units.fs,
    sim.properties("potential"),
    sim.properties("temperature"),
)

# `get_structures` dumps the state of the system as ASE Atoms objects, possibly listing
# all systems and beads
ase.io.write("final_positions.xyz", sim.get_structures())

# we can also set the simulation state
structure = sim.get_structures()
structure.positions[:] = data.positions
structure.arrays["ipi_velocities"][:] = 0
sim.set_structures(structure)
print(sim.properties("potential"), sim.properties("kinetic_md"))

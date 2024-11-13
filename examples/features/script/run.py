#!/usr/bin/python
from ipi.utils.scripting import setup_forcefield, setup_nvt, InteractiveSimulation
import ase.io

# There are utilities to quickly set up XML inputs for commonly-used simulations
data = ase.io.read("water-64.xyz")
input_xml = setup_nvt(
    structures=data,
    forcefield=setup_forcefield(
        name="harmonic", mode="direct", pes="harmonic", parameters="{k1:1e-3}"
    ),
    temperature=400,
    timestep=0.5,
)

# The main object for scripting is `InteractiveSimulation`, that is initialized from
# and XML input and acts as a wrapper around an i-PI simulation object
sim = InteractiveSimulation(input_xml)

# `properties` accesses the (well) properties of the simulation object
print(sim.properties("potential"), sim.properties("kinetic_md"))
# `run` advances the interactive simulation by one (or the prescribed number) of steps
sim.run(10)
print(sim.properties("potential"), sim.properties("kinetic_md"))
sim.run(10)
print(sim.properties("potential"), sim.properties("kinetic_md"))

# `get_structures` dumps the state of the system as ASE Atoms objects, possibly listing
# all systems and beads
ase.io.write("final_positions.xyz", sim.get_structures())

# we can also set the simulation state
structure = sim.get_structures()
structure.positions[:] = data.positions
structure.arrays["ipi_velocities"][:] = 0
sim.set_structures(structure)
print(sim.properties("potential"), sim.properties("kinetic_md"))

#!/usr/bin/python
from ipi.utils.scripting import InteractiveSimulation
from ipi.utils.depend import dstrip
import numpy as np
import ase, ase.io

# Scripting simulations can also be initialized from an external XML template
with open("template.xml", "r", encoding="utf-8") as file:
    input_xml = file.read()

print("Running with XML input:\n\n", input_xml)


# It is possible to define custom properties by passing a dictionary
# of property definitions to the InteractiveSimulation constructor.
# Here we define a custom property "fancy_virial" that computes the
# scalar product of posisions and momenta as a dummy example.
def my_virial(self):
    # self is a reference to the `properties` object of a system, and
    # provides access to all physical quantities defining the system
    virial = np.dot(
        dstrip(self.beads.q).flatten(), dstrip(self.beads.p).flatten()
    )  # dummy operation to use positions
    return virial * 0.5


# and this returns a 3-vector
def my_virial_xyz(self):
    q = dstrip(self.beads.q).reshape(-1, 3)
    p = dstrip(self.beads.p).reshape(-1, 3)
    virial = (p * q).sum(axis=0)

    return virial


# in order to have a class method used as a property
# (useful e.g. for caching) we need to be a bit more elaborate
class MyPropertyCalculator:
    def __init__(self, scaling_factor):
        # can initialize some stuff in here
        self.scaling_factor = scaling_factor

    def class_virial(self, properties):
        return self.scaling_factor * my_virial(properties)


# we create an instance of the class...
my_calculator = MyPropertyCalculator(scaling_factor=2.0)
# and use a lambda to wrap the method so that it complies with the
# expected function signature
my_virial_class = lambda properties: my_calculator.class_virial(properties)

# properties can be just a function, or a full dictionary defining size, help string, etc.
props = {
    "fancy_virial": my_virial,
    "fancy_virial_xyz": {
        "func": my_virial_xyz,
        "dimension": "action",
        "size": 3,
        "help": "A dummy function returning a 3-vector",
    },
    "class_virial": my_virial_class,
}
sim = InteractiveSimulation(input_xml, custom_properties=props)

# `properties` accesses the (well) properties of the simulation object
print("# Initial properties")
print(
    f'Time:        {sim.properties("time") / ase.units.fs}  '
    f'Potential:   {sim.properties("potential")}  '
    f'Temperature: {sim.properties("temperature")}  '
    f'Virial:      {sim.properties("fancy_virial")}  '
    f'Virial_xyz:  {sim.properties("fancy_virial_xyz")}  '
    f'Class_virial: {sim.properties("class_virial")}  '
)
# `run` advances the interactive simulation by one (or the prescribed number) of steps
sim.run(10)
print("# After 10 steps")
print(
    f'Time:        {sim.properties("time") / ase.units.fs}  '
    f'Potential:   {sim.properties("potential")}  '
    f'Temperature: {sim.properties("temperature")}  '
    f'Virial:      {sim.properties("fancy_virial")}  '
    f'Virial_xyz:  {sim.properties("fancy_virial_xyz")}  '
    f'Class_virial: {sim.properties("class_virial")}  '
)

Scripting interface
===================

This is a simple demonstration of using i-PI within a Python script.
The `utils.scripting` module contain some helper functions to construct
an XML input file, and a `InteractiveSimulation` object that provides
an interface to run, and access the properties, of an i-PI simulation.
Any valid XML input file can also be used.

This example uses a simple `<ffdirect>` PES, but it is also possible 
to run with an external driver, that has to be launched after the 
simulation is initialized, but before it is run. 

The class outputs files if required in the XML file, and unless it is
disabled in the `run` method,  but also provides a way to dump the 
configurations and compute the properties interactively. 
ALL input and output properties here are taken to use ASE units, since
it is expected that this API will be most useful to users who are
familiar with an ASE workflow.

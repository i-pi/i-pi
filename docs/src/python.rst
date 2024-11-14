Tools and python utilities
==========================

Post-processing tools
~~~~~~~~~~~~~~~~~~~~~
Some observables (such as autocorrelation functions, isotope fractionation ratios, particle
momentum distributions) require post-processing the raw outputs of a simulation. 
The folder ``tools`` contains several post-processing scripts, that can 
be used to this end - each (tersely) self-documented with a help string and in
the code. 
Some of these tools are also accessible as python functions, so that they can 
be invoked from custom post-processing workflows. They can be found in 
the ``ipi/utils/tools`` folder, and accessed by importing the module, e.g.

.. code-block:: python
        
    from ipi.utils.tools import isra_deconvolute

is a routine to correct the correlation spectrum computed from a thermostatted simulation,
see this `example <https://atomistic-cookbook.org/examples/thermostats/thermostats.html>`_.

Parsing
~~~~~~~

i-PI trajectories can be output in 
`extended xyz format <https://wiki.fysik.dtu.dk/ase/ase/io/formatoptions.html#extxyz>`_, 
and can be read by `ASE <https://wiki.fysik.dtu.dk/ase/index.html>`_, while the property
outputs can be simply read with ``np.loadtxt``. However, more reliable parsing can be
achieved using the ``ipi.read_output`` and ``ipi.read_trajectory`` functions. 

Scripting i-PI
~~~~~~~~~~~~~~

If one wants to run an i-PI simulation within a Python script, it is also possible
to use a (somewhat primitive) scripting API, defined in the ``ipi.utils.scripting``
module. The core component is the ``InteractiveSimulation`` class, that can be
initialized from an *XML* input, advanced for a given number of steps using the
``run`` method. The calculation requires also the use of a driver, that can 
communicate through sockets (in which case it must be launched after 
initialization, and before running) or using an :ref:`ffdirect` block. 

Properties can be accessed using the ``properties`` method, and a snapshot
of the configuration by calling ``get_structures``. 
Several utility functions, also available within the module, facilitate
building on the fly an *XML* string to initialize the simulation.
An example of usage of this interface goes as follows:

.. code-block:: python 

    data = ase.io.read("init.xyz")
    input_xml = simulation_xml(
        structures=data,
        forcefield=forcefield_xml(name="dummy", mode="direct"),
        motion=motion_nvt_xml(timestep=2.0 * ase.units.fs),
        temperature=300,
        prefix="example",
    )

    sim = InteractiveSimulation(input_xml)
    sim.run(100)
    potential = sim.properties("potential")
    ase.io.write("final_structure.xyz", sim.get_structures())





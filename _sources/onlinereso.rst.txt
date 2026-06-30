.. _librarywebsites:

Tutorials, recipes, and on-line resources
=========================================

Massive Open Online Course (MOOC)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The theory behind classical and path-integral methods for structure 
and dynamics along with practical exercises can be found in an online 
course freely available to everyone, and accessible from this
`link <https://courseware.epfl.ch/courses/course-v1:EPFL+X+2022/about>`_.

i-PI resources
~~~~~~~~~~~~~~

For more information about i-PI and to download the source code go to
http://ipi-code.org/.

In http://gle4md.org/ one can also obtain colored-noise parameters to
run Path Integral with Generalized Langevin Equation thermostat
(PI+GLE/PIGLET) calculations.

Examples and demos 
~~~~~~~~~~~~~~~~~~~
The `examples` and `demos` folders in the i-PI `source code <https://github.com/i-pi/i-pi>`_
contain inputs for many different types of
calculations based on i-PI. Examples are typically minimal use-cases of specific
features, examples of client codes, or setups for high performance computing, 
while demos are more structured, tutorial-like examples that show how
to realize more complex setups, and also provide a brief discussion of the 
underlying algorithms.

To run these examples, you should typically start i-PI, redirecting the output to
a log file, and then run a couple of instances of the driver code. The progress
of the wrapper is followed by monitoring the log file with the `tail` Linux command.

Optionally, you can make a copy of the directory with the example somewhere
else if you want to keep the i-PI directory clean. For example::

    bash
    cd demos/para-h2-tutorial/tutorial-1/
    i-pi tutorial-1.xml > log &
    i-pi-driver -a localhost -p 31415 -m sg -o 15 &
    i-pi-driver -a localhost -p 31415 -m sg -o 15 &
    tail -f log

Tutorials
~~~~~~~~~
A simple tutorial on how to run i-PI can be found in the documentation:
`a simple tutorial <tutorials.html>`_.

You can try a set of guided examples that demonstrate some more 
advanced features of i-PI. These Tutorials, from the 2021 CECAM 
Flagship school, are available here:
`i-PI 2021 tutorials <https://github.com/i-pi/tutorials-schools/>`_.

Tutorials from the 2023 CECAM Flagship school are available here:
`i-PI 2023 tutorials <https://github.com/i-pi/piqm2023-tutorial>`_.

Tutorial from one of our developers on using machine learning interatomic potentials with i-PI are available here:
`MLIPs-with-i-PI <https://github.com/venkatkapil24/MLIPs-with-iPI>`_.

Note that these examples use some features 
(such as qTIP4P/F water and the Zundel cation potentials in the driver code, 
and optionally also CP2K). Instructions on how to install and run these codes 
with i-PI are contained in the tutorial.

Atomistic recipes
~~~~~~~~~~~~~~~~~

Several examples of usage of i-PI to perform advanced molecular 
simulations can be found in the 
`i-PI chapter <https://atomistic-cookbook.org/software/i-pi.html>`_ of the 
`atomistic cookbook <https://atomistic-cookbook.org/>`_.
These include an 
`introduction to path integral simulations <https://atomistic-cookbook.org/latest/examples/path-integrals/path-integrals.html>`_, 
a discussion of `how to compute quantum heat capacities <https://atomistic-cookbook.org/latest/examples/heat-capacity/heat-capacity.html>`_,
and an example of `path integral metadynamics <https://atomistic-cookbook.org/latest/examples/pi-metad/pi-metad.html>`_.


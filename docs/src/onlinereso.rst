.. _librarywebsites:

On-line resources
=================

Massive Open Online Course (MOOC)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The theory behind many types of simulation that i-PI can do and examples 
of simulations with i-PI can also be found in an online course freely 
available to everyone, and accessible from this 
`link <https://courseware.epfl.ch/courses/course-v1:EPFL+X+2022/about>`_.

i-PI resources
~~~~~~~~~~~~~~

For more information about i-PI and to download the source code go to
http://ipi-code.org/.

In http://gle4md.org/ one can also obtain colored-noise parameters to
run Path Integral with Generalized Langevin Equation thermostat
(PI+GLE/PIGLET) calculations.

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

Note that these examples use some features 
(such as qTIP4P/F water and the Zundel cation potentials in the driver code, 
and optionally also CP2K). Instructions on how to install and run these codes 
with i-PI are contained in the tutorial.

Atomistic recipes
~~~~~~~~~~~~~~~~~

Several examples of usage of i-PI to perform advanced molecular 
simulations can be found among the recipes of the 
`atomistic cookbook <https://atomistic-cookbook.org/>`_.
These include an 
`introduction to path integral simulations <https://atomistic-cookbook.org/latest/examples/path-integrals/path-integrals.html>`_, 
a discussion of `how to compute quantum heat capacities <https://atomistic-cookbook.org/latest/examples/heat-capacity/heat-capacity.html>`_,
and an example of `path integral metadynamics<https://atomistic-cookbook.org/latest/examples/pi-metad/pi-metad.html>`_.

Client code resources
~~~~~~~~~~~~~~~~~~~~~

Several codes provide out-of-the-box an i-PI interface, including 
ASE, 
CASTEP, 
CP2K,
DFTB+,
elphmod,
ffsGDML,
FHI-aims, 
LAMMPS, 
librascal, 
Quantum ESPRESSO, 
Siesta,
Yaff.

If you are interested in interfacing your code to i-PI please get in
touch, we are always glad to help!

There are several Fortran and C libraries that most client codes will
probably need to run, such as FFTW, BLAS and LAPACK. These can be found
at `www.fftw.org <http://www.fftw.org>`__,
`www.netlib.org/blas <http://www.netlib.org/blas>`__ and
`www.netlib.org/lapack <http://www.netlib.org/lapack>`__ respectively.

These codes do not come as part of the i-PI package, and must be
downloaded separately. 

.. See chapter :ref:`clientinstall` for more details of how to do this. 


==============================
i-PI: a Universal Force Engine
==============================

A Python interface for ab initio path integral molecular dynamics simulations.
i-PI is composed of a Python server (i-pi itself, that does not need to be
compiled but only requires a relatively recent version of Python and Numpy)
that propagates the (path integral) dynamics of the nuclei, and of an external
code that acts as a client and computes the electronic energy and forces.

This is typically a patched version of an electronic structure code, but a
simple self-contained Fortran driver that implements Lennard-Jones and
Silvera-Goldman potentials is included for test purposes.


Quick Setup and Test
====================

To use i-PI with already existing drivers, install and update using Pip::

   $ pip install -U i-PI

Test with Pytest::

   $ pip install pytest
   $ pytest --pyargs ipi.tests


Full installation
=================

To develop i-PI or test it with the self-contained driver, follow these
instructions. It is assumed that i-PI will
be run from a Linux environment, with a recent version of Python, Numpy and
gfortran, and that the terminal is initially in the i-pi package directory (the
directory containing this file), which you can obtain by cloning the repository::

   $ git clone https://github.com/i-pi/i-pi.git

Source the environment settings file :code:`env.sh` as :code:`$ source env.sh` or :code:`$ .
env.sh`.  It is useful to put this in your :code:`.bashrc` or other settings file if
you always want to have i-PI available.


Compile the driver code
-----------------------

::

  $ cd drivers
  $ make
  $ cd ..


Run one of the examples
-----------------------

This will first start the wrapper in the background, redirecting the output to
a log file, and then run a couple of instances of the driver code. The progress
of the wrapper is followed by monitoring the log file with the `tail` Linux
command.

Optionally, you can make a copy of the directory with the example somewhere
else if you want to keep the i-PI directory clean.

::

  $ cd examples/tutorial/tutorial-1/
  $ i-pi tutorial-1.xml > log &
  $ i-pi-driver -h localhost -p 31415 -m sg -o 15 &
  $ i-pi-driver -h localhost -p 31415 -m sg -o 15 &
  $ tail -f log

The monitoring can be interrupted with CTRL+C when the run has finished (5000 steps).


Run the automatic test suite
----------------------------

The automatic test suite can be run by calling the i-pi-test script. 
You need to have the `pytest` package installed

::

  $ i-pi-test

See more details in the README file inside the ipi_tests folder.

Contributing
================

i-PI is an open source project, and everyone is welcome to contribute
with bug fixes, documentation, examples, new features, regression tests, etc.

Your contribution should be based on the master branch. We kindly ask you to first fork the project,
make your changes, make sure you comply with all items in our checklist below, and finally create a pull request (PR).

Checklist to create a pull request:

- The PR follows our format compliance (based on `black` and `flake8` as explained above)
- All the new classes and functions include the corresponding docstrings

(If the PR adds a new functionally, please fulfill the next two requirements as well)

- Add a working example to the `examples` folder to showcase the new functionality
- Add a regression test to the `i-pi/ipi_tests/regression_tests` folder (see the corresponding README file for further details)
- Make sure that all the automatic checks pass without any error

We are looking forward to your contribution!

Format Compliance
-----------------

i-PI code should be compliant to a minimal subset of PEP-8 recommendations.
Currently, we require the use of `black` as formatter and linter.
We also ask for the usage of `flake8` for syntactic checks, which is also
part of linting.
In most systems, both packages can be easily installed using `pip`.
BEFORE proceeding to a pull request, the minimal requirement is that you run

::

  $ make lint
  $ make pretty 

This will ensure the formatting and linting requirement are applied in the whole 
directory tree. Please resolve any warnings or errors that may appear. Your
commit will not pass the CI tests otherwise.

For a more flexible setup, we also provide the script `i-pi-style`, for
which instructions can be obtained by typing 

::

  $ i-pi-style -h 

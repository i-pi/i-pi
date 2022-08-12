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
--------------------

To use i-PI with already existing drivers, install and update using Pip:

Last version::

   $ python -m pip install git+https://github.com/i-pi/i-pi.git

Last Release::

   $ pip install -U i-PI

Test with Pytest::

   $ pip install pytest
   $ pytest --pyargs ipi.tests


Full installation
-----------------

To develop i-PI or test it with the self-contained driver, follow these
instructions. It is assumed that i-PI will
be run from a Linux environment, with a recent version of Python, Numpy and
gfortran, and that the terminal is initially in the i-pi package directory (the
directory containing this file), which you can obtain by cloning the repository::

   $ git clone https://github.com/i-pi/i-pi.git

Source the environment settings file :code:`env.sh` as :code:`$ source env.sh` or :code:`$ .
env.sh`.  It is useful to put this in your :code:`.bashrc` or other settings file if
you always want to have i-PI available.


### Compile the driver code

::

  $ cd drivers
  $ make
  $ cd ..


### Run one of the examples

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


### Run the automatic test suite

The automatic test suite can be run by calling the i-pi-test script.
You need to have the `pytest` package installed

::

  $ i-pi-test

See more details in the README file inside the ipi_tests folder.

Contributing
------------

If you have new features you want to implement into i-PI, your contributions are much welcome.
See `CONTRIBUTING.md` for a brief set of style guidelines and best practices. Before embarking
into a substantial project, it might be good to get in touch with the developers, e.g. by opening
a wishlist issue.

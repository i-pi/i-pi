i-PI: a Universal Force Engine
==============================

A Python interface for ab initio path integral molecular dynamics simulations.
i-PI is composed of a Python server (i-pi itself, that does not need to be
compiled but only requires a relatively recent version of Python and Numpy)
that apply an algorithm that updates the positions of the nuclei, and of an external
code that acts as a client and computes the electronic energy and forces.

This is typically a patched version of an electronic structure code, but a
simple self-contained Fortran driver that implements Lennard-Jones and
Silvera-Goldman potentials is included for test purposes.

i-PI was originally developed to simulate the quantum mechanical nature of light
nuclei by performing path integral molecular dynamics simulations,
and it implements most of the state-of-the-art methods to accelerate this kind of 
calculations. It has since grown to also provide all sorts of simulation 
strategies, from replica exchange to geometry optimization. 

Quick Setup and Test
--------------------

To use i-PI with already existing drivers, install and update using Pip:

Last version::
```bash
python -m pip install git+https://github.com/i-pi/i-pi.git
```

Last Release::
```bash
pip install -U i-PI
```

Test with Pytest::
```bash
pip install pytest
pytest --pyargs ipi.tests
```

Full installation
-----------------

To develop i-PI or test it with the self-contained driver, follow these
instructions. It is assumed that i-PI will
be run from a Linux environment, with a recent version of Python, Numpy and
gfortran, and that the terminal is initially in the i-pi package directory (the
directory containing this file), which you can obtain by cloning the repository
```bash
git clone https://github.com/i-pi/i-pi.git
```

Source the environment settings file `env.sh` as `source env.sh` or `.
env.sh`.  It is useful to put this in your `.bashrc` or other settings file if
you always want to have i-PI available.


### Compile the driver code

```bash
cd drivers/f90
make
cd ../..
```

### Examples and demos

The `examples` and `demos` folders contain inputs for many different types of
calculations based on i-PI. Examples are typically minimal use-cases of specific
features, while demos are more structured, tutorial-like examples that show how
to realize more complex setups, and also provide a brief discussion of the 
underlying algorithms.

To run these examples, you should typically start i-PI, redirecting the output to
a log file, and then run a couple of instances of the driver code. The progress
of the wrapper is followed by monitoring the log file with the `tail` Linux command.

Optionally, you can make a copy of the directory with the example somewhere
else if you want to keep the i-PI directory clean. For example

```bash
cd demos/para-h2-tutorial/tutorial-1/
i-pi tutorial-1.xml > log &
i-pi-driver -h localhost -p 31415 -m sg -o 15 &
i-pi-driver -h localhost -p 31415 -m sg -o 15 &
tail -f log
```

The monitoring can be interrupted with CTRL+C when the run has finished (5000 steps).


### Run the automatic test suite

The automatic test suite can be run by calling the i-pi-test script.
You need to have the `pytest` package installed

```
i-pi-test
```

See more details in the README file inside the ipi_tests folder.

Contributing
------------

If you have new features you want to implement into i-PI, your contributions are much welcome.
See `CONTRIBUTING.md` for a brief set of style guidelines and best practices. Before embarking
into a substantial project, it might be good to get in touch with the developers, e.g. by opening
a wishlist issue.

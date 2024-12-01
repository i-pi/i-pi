# i-PI: a Universal Force Engine

A Python interface for ab initio path integral molecular dynamics simulations (and more).
i-PI is a Python server (that does not need to be compiled and only requires a relatively
recent version of Python and Numpy) that applies an algorithm to update the positions of
the nuclei. One of many compatible external codes acts as client, and computes the
electronic energy and forces.

This is typically a patched version of an electronic structure code, but a
simple self-contained Fortran driver that implements several simple interatomic
potentials is included for test purposes.

i-PI was originally developed to simulate the quantum mechanical nature of light
nuclei by performing path integral molecular dynamics simulations,
and it implements most of the state-of-the-art methods to accelerate this kind of
calculations. It has since grown to also provide all sorts of simulation
strategies, from replica exchange to geometry optimization.

If you use i-PI in your research, please cite the accompanying publication:
for version 3, the relevant paper is
[Litman et al., _J. Chem. Phys._ 161, 062504 (2024)](https://doi.org/10.1063/5.0215869)

```
@article{litman2024ipi,
title={i-PI 3.0: a flexible and efficient framework for advanced atomistic simulations},
author={Yair Litman and Venkat Kapil and Yotam M. Y. Feldman and Davide Tisi and Tomislav Begušić and Karen Fidanyan and Guillaume Fraux and Jacob Higer and Matthias Kellner and Tao E. Li and Eszter S. Pós and Elia Stocco and George Trenins and Barak Hirshberg and Mariana Rossi and Michele Ceriotti},
journal = {J. Chem. Phys.},
pages   = {062505},
volume  = {161},
year    = {2024}
}
```

## Quick Setup

To use i-PI with an existing driver, install and update using `pip`:

Last version:

```bash
python -m pip install git+https://github.com/i-pi/i-pi.git
```

Last Release:

```bash
pip install -U ipi
```

## Documentation 

You can find the online documentation at [https://docs.ipi-code.org](https://docs.ipi-code.org/). Alternatively, you can build it locally by following instructions in the `docs/README.md` file.

## Source installation

To develop i-PI or test it with the self-contained driver, follow these
instructions. It is assumed that i-PI will
be run from a Linux environment, with a recent version of Python, Numpy and
gfortran, and that the terminal is initially in the i-pi package directory (the
directory containing this file), which you can obtain by cloning the repository

```bash
git clone https://github.com/i-pi/i-pi.git
```

Source the environment settings file `env.sh` as `source env.sh` or `.
env.sh`. It is useful to put this in your `.bashrc` or other settings file if
you always want to have i-PI available.

## Compile the driver code

The built-in driver requires a FORTRAN compiler, and can be built as

```bash
cd drivers/f90
make
cd ../..
```

There is also a Python driver available in `drivers/py`, which however has limited
functionalities.

## Examples and demos

The `examples` and `demos` folders contain inputs for many different types of
calculations based on i-PI. Examples are typically minimal use-cases of specific
features, while demos are more structured, tutorial-like examples that show how
to realize more complex setups, and also provide a brief discussion of the
underlying algorithms.

To run these examples, you should typically start i-PI, redirecting the output to
a log file, and then run a couple of instances of the driver code. The progress
of the wrapper is followed by monitoring the log file with the `tail` Linux command.

Optionally, you can make a copy of the directory with the example somewhere
else if you want to keep the i-PI directory clean. For example, after sourcing the `env.sh` file, 

```bash
cd demos/para-h2-tutorial/tutorial-1/
i-pi tutorial-1.xml > log &
i-pi-driver -a localhost -p 31415 -m sg -o 15 &
i-pi-driver -a localhost -p 31415 -m sg -o 15 &
tail -f log
```

The monitoring can be interrupted with CTRL+C when the run has finished (5000 steps).

## Tutorials and online resources

The i-PI [documentation](https://docs.ipi-code.org/onlinereso.html) has a list of
available tutorials, recipes and other useful online resources.

## Run the automatic test suite

The automatic test suite can be run by calling the i-pi-tests script.
You need to have the `pytest` package installed

```
i-pi-tests
```

You may also need to install some dependencies, listed in `requirements.txt`.

See more details in the README file inside the `ipi_tests` folder.

## Contributing

If you have new features you want to implement into i-PI, your contributions are much welcome.
See `CONTRIBUTING.md` for a brief set of style guidelines and best practices. Before embarking
into a substantial project, it might be good to get in touch with the developers, e.g. by opening
a wishlist issue.

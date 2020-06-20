#!/usr/bin/env python2

import os
import sys

import ase.io
from ase.calculators.gaussian import Gaussian

from ipi.interfaces.clients import ClientASE


# get command-line arguments
try:
    nthreads = int(sys.argv[1])
    dir_calc = sys.argv[2]
except IndexError:
    print(
        "command line arguments: <number of Gaussian threads> <calculator directory name>"
    )
    sys.exit(1)


# print some info
print("Gaussian ASE i-PI socket runner")
print("-------------------------------")
print()
print("   number of threads: {:d}".format(nthreads))
print("calculator directory: {:s}".format(dir_calc))
print()

# fast scratch directory for Gaussian
os.environ["GAUSS_SCRDIR"] = "/dev/shm"

# create ASE atoms and calculator
atoms = ase.io.read("init.xyz")
atoms.set_calculator(
    Gaussian(
        method="MP2",
        # Note that 6-31G is for faster testing only, not a good basis set
        # for gas-phase W2. Use the bigger one instead.
        basis="6-31G",
        # basis = '6-311++G**',
        nproc=nthreads,
        label=dir_calc + "/calculator",
    )
)

# create the socket client...
client = ClientASE(atoms, address="ase")
# client = ClientASE(atoms, mode='inet', address='localhost', port=12345)

# ... and run it
client.run(t_max=35, fn_exit="EXIT_ASE")  # test max run time feature
client.run(fn_exit="EXIT_ASE")

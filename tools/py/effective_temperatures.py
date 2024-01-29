#!/usr/bin/python

__author__ = "Igor Poltavsky"

""" effective_temperatures.py
The script reads the forces from from standard i-PI output files and computes the effective temperatures
for all degrees of freedom. The output is saved to 'prefix.effective_temperatures.dat' file which is located
in the folder which contains the input files. The results are printed out in the format: "atom",
"effective temperatures for x, y, and z directions", and "average effective temperature".

The script assumes that the input files are in 'xyz' format and prefix.for_*.xyz (forces)
naming scheme. This requires the following lines in input.xml file:
<trajectory filename='force' stride='1' format='xyz' cell_units='angstrom'>
forces </trajectory>

Syntax:
   python effective_temperatures.py "prefix" "simulation temperature (in Kelvin)"
   "number of time frames to skip in the beginning of each file (default 0)"

The computed effective temperatures are in Kelvin.
"""

import numpy as np
import sys
import glob

from ipi.utils.units import unit_to_internal, unit_to_user, Constants
from ipi.utils.io import read_file


def effectiveTemperatures(prefix, temp, ss=0):
    """
    Computes effective temperatures for a given (PI)MD dynamics.
    """

    temperature = unit_to_internal(
        "temperature", "kelvin", float(temp)
    )  # simulation temperature
    skipSteps = int(ss)  # steps to skip for thermalization

    f2_av = None  # average square forces array
    names = None

    fns_for = sorted(glob.glob(prefix + ".for*"))
    fn_out = prefix + ".effective_temperatures.dat"

    nbeads = len(fns_for)

    # print some information
    print("temperature = {:f} K".format(float(temp)))
    print()
    print("number of beads = {:d}".format(nbeads))
    print()
    print("forces file names:")
    for fn_for in fns_for:
        print("{:s}".format(fn_for))
    print()
    print("output file name:")
    print(fn_out)
    print()

    # open input and output files
    ifor = [open(fn, "r") for fn in fns_for]
    iOut = open(fn_out, "w")

    # Some constants
    const = Constants.hbar**2 / (12.0 * (nbeads * Constants.kb * temperature) ** 3)

    iOut.write(
        "# Atom, Cartesian components of the effective temperature, average effective temperature (in Kelvin)\n"
    )

    natoms = 0
    ifr = 0
    f, m = None, None
    while True:  # Reading input files and calculating PPI correction
        if ifr % 100 == 0:
            print("\rProcessing frame {:d}".format(ifr), end=" ")
            sys.stdout.flush()

        try:
            for i in range(nbeads):
                ret = read_file("xyz", ifor[i], dimension="force")["atoms"]
                if natoms == 0:
                    m, natoms, names = ret.m, ret.natoms, ret.names
                    f = np.zeros((nbeads, 3 * natoms))
                    f2_av = np.zeros(3 * natoms)
                f[i, :] = ret.q
        except EOFError:  # finished reading files
            break

        if ifr >= skipSteps:  # PPI correction
            f2 = np.zeros(3 * natoms)
            for i in range(natoms):
                for j in range(nbeads):
                    f2[i * 3 : i * 3 + 3] += f[j, i * 3 : i * 3 + 3] ** 2 / m[i]

            f2_av[:] += f2[:]
            ifr += 1

        else:
            ifr += 1

    dT = const * f2_av / float(ifr - skipSteps)

    temperature = unit_to_user("temperature", "kelvin", temperature)

    for i in range(natoms):
        iOut.write(
            "%s    %f     %f     %f     %f\n"
            % (
                names[i],
                temperature * (1 + dT[i * 3]),
                temperature * (1 + dT[i * 3 + 1]),
                temperature * (1 + dT[i * 3 + 2]),
                temperature * (1 + np.sum(dT[i * 3 : i * 3 + 3]) / 3.0),
            )
        )

    for f in ifor:
        f.close()

    iOut.close()


def main(*arg):
    effectiveTemperatures(*arg)


if __name__ == "__main__":
    main(*sys.argv[1:])

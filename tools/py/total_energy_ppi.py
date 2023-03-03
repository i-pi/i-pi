#!/usr/bin/python

__author__ = "Igor Poltavsky"
__version__ = "1.0"

""" total_energy_ppi.py
Reads simulation time, potential energy, positions and forces from
an i-PI run and computes a virial total energy estimator and a ppi correction
for each time frame. The output is saved to 'prefix.energy.dat' file which is created in the
folder which contains the input files. The results are printed out in the format: "time
frame", "virial total energy estimator", and "ppi correction".

The script assumes that the input files are in 'xyz' format, with prefix.out (contains simulation time and
potential energy among other output properties), prefix.pos_*.xyz (positions) and
prefix.for_*.xyz (forces) naming scheme.
This would require the following lines in input.xml file:
<properties filename='out' stride='n'> [step, time{picosecond}, potential{kelvin}] </properties>
<trajectory filename='pos' stride='n' format='xyz' cell_units='angstrom'> positions{angstrom} </trajectory>
<trajectory filename='for' stride='n' format='xyz' cell_units='angstrom'> forces{piconewton} </trajectory>
where n is the same integer number.

Syntax:
   python total_energy_ppi.py "prefix" "simulation temperature (in Kelvin)" "number of time frames to skip
   in the beginning of each file (default 0)" "units for the energy output (default choice is the units for
   the potential energy in input file)"

Currently supported energy units are: atomic_unit, electronvolt, j/mol, cal/mol, and kelvin

"""

import numpy as np
import sys
import glob
import os
import re

from ipi.utils.units import unit_to_internal, unit_to_user, Constants
from ipi.utils.io import read_file


def totalEnergy(prefix, temp, ss=0, unit=""):
    """
    Computes the virial centroid estimator for the total energy and PPI correction.
    """

    # Adding fortran functions (when exist)
    sys.path.append(os.path.abspath(os.path.dirname(sys.argv[0]))[:-2] + "f90")
    fast_code = True
    try:
        import fortran
    except ImportError:
        fast_code = False
        print(
            "WARNING: No compiled fortran module for fast calculations have been found.\n"
            "Calculations will use a slower python script."
        )

    temperature = unit_to_internal(
        "temperature", "kelvin", float(temp)
    )  # simulation temperature
    skipSteps = int(ss)  # steps to skip for thermalization

    f2_av, ePA_av, eVir_av, f2ePA_av = 0.0, 0.0, 0.0, 0.0  # some required sums

    fns_pos = sorted(glob.glob(prefix + ".pos*"))
    fns_for = sorted(glob.glob(prefix + ".for*"))
    fns_iU = glob.glob(prefix + ".out")[0]
    fn_out_en = prefix + ".total_energy.dat"

    # check that we found the same number of positions and forces files
    nbeads = len(fns_pos)
    if nbeads != len(fns_for):
        print(fns_pos)
        print(fns_for)
        raise ValueError(
            "Mismatch between number of input files for forces and positions."
        )

    # print some information
    print("temperature = {:f} K".format(float(temp)))
    print()
    print("number of beads = {:d}".format(nbeads))
    print()
    print("positions and forces file names:")
    for fn_pos, fn_for in zip(fns_pos, fns_for):
        print("{:s}   {:s}".format(fn_pos, fn_for))
    print()
    print("potential energy file: {:s}".format(fns_iU))
    print()
    print("output file name:")
    print(fn_out_en)
    print()

    # open input and output files
    ipos = [open(fn, "r") for fn in fns_pos]
    ifor = [open(fn, "r") for fn in fns_for]
    iU = open(fns_iU, "r")
    iE = open(fn_out_en, "w")

    # Some constants
    beta = 1.0 / (Constants.kb * temperature)
    const_1 = 0.5 * nbeads / (beta * Constants.hbar) ** 2
    const_2 = 1.5 * nbeads / beta
    const_3 = 1.5 / beta
    const_4 = Constants.kb**2 / Constants.hbar**2
    const_5 = Constants.hbar**2 * beta**3 / (24.0 * nbeads**3)

    timeUnit, potentialEnergyUnit, potentialEnergy_index, time_index = extractUnits(
        iU
    )  # extracting simulation time
    # and potential energy units

    # Defining the output energy unit
    if unit == "":
        unit = potentialEnergyUnit

    iE.write(
        "# Simulation time (in %s), virial total energy and PPI energy correction (in %s)\n"
        % (timeUnit, unit)
    )

    natoms = 0
    ifr = 0
    time0 = 0
    q, f, m = None, None, None
    while True:  # Reading input files and calculating PPI correction
        if ifr % 100 == 0:
            print("\rProcessing frame {:d}".format(ifr), end=" ")
            sys.stdout.flush()

        try:
            for i in range(nbeads):
                ret = read_file("xyz", ipos[i], output="arrays")
                if natoms == 0:
                    m, natoms = ret["masses"], ret["natoms"]
                    q = np.zeros((nbeads, 3 * natoms))
                    f = np.zeros((nbeads, 3 * natoms))
                q[i, :] = ret["data"]
                f[i, :] = read_file("xyz", ifor[i], output="arrays")["data"]
            U, time = read_U(iU, potentialEnergyUnit, potentialEnergy_index, time_index)
        except EOFError:  # finished reading files
            sys.exit(0)

        if ifr < skipSteps:
            time0 = time

        if ifr >= skipSteps:  # PPI correction
            time -= time0

            ePA, f2, f2ePA, eVir = 0.0, 0.0, 0.0, 0.0

            if not fast_code:
                for j in range(nbeads):
                    for i in range(natoms):
                        f2 += (
                            np.dot(f[j, i * 3 : i * 3 + 3], f[j, i * 3 : i * 3 + 3])
                            / m[i]
                        )
                for i in range(natoms):
                    ePA -= (
                        np.dot(
                            q[0, i * 3 : i * 3 + 3] - q[nbeads - 1, i * 3 : i * 3 + 3],
                            q[0, i * 3 : i * 3 + 3] - q[nbeads - 1, i * 3 : i * 3 + 3],
                        )
                        * m[i]
                    )
                for j in range(nbeads - 1):
                    for i in range(natoms):
                        ePA -= (
                            np.dot(
                                q[j + 1, i * 3 : i * 3 + 3] - q[j, i * 3 : i * 3 + 3],
                                q[j + 1, i * 3 : i * 3 + 3] - q[j, i * 3 : i * 3 + 3],
                            )
                            * m[i]
                        )
                rc = np.zeros(3)
                for i in range(natoms):
                    rc[:] = 0.0
                    for j in range(nbeads):
                        rc[:] += q[j, i * 3 : i * 3 + 3]
                    rc[:] /= nbeads
                    for j in range(nbeads):
                        eVir += np.dot(
                            rc[:] - q[j, i * 3 : i * 3 + 3], f[j, i * 3 : i * 3 + 3]
                        )

                ePA *= const_1
                ePA += const_2 * natoms + U
                f2ePA = f2 * ePA
                eVir /= 2.0 * nbeads
                eVir += const_3 * natoms + U

            else:
                f2 = fortran.f2divm(
                    np.array(f, order="F"), np.array(m, order="F"), natoms, nbeads
                )
                ePA = fortran.findcoupling(
                    np.array(q, order="F"),
                    np.array(m, order="F"),
                    temperature,
                    natoms,
                    nbeads,
                )
                eVir = fortran.findcentroidvirialkineticenergy(
                    np.array(f, order="F"), np.array(q, order="F"), natoms, nbeads
                )

                ePA *= const_4
                ePA += const_2 * natoms + U
                f2ePA = f2 * ePA
                eVir += const_3 * natoms + U

            ePA_av += ePA
            f2_av += f2
            f2ePA_av += f2ePA
            eVir_av += eVir
            ifr += 1

            norm = float(ifr - skipSteps)

            dE = (
                3.0 * Constants.kb * temperature + ePA_av / norm
            ) * f2_av / norm - f2ePA_av / norm
            dE *= const_5

            dE = unit_to_user("energy", unit, dE)
            eVir = unit_to_user("energy", unit, eVir_av / norm)

            iE.write("%f    %f     %f\n" % (time, eVir, dE))

        else:
            ifr += 1


def extractUnits(filedescU):
    """
    Extracting potential energy and time step units.
    Also this function looking for the simulation time
    and potential energy position in the prefix.out
    file. Thus, this file can contain any number of
    output properties in arbitrary ordering.

    Args:
       filedesc: An open readable file object from a xyz formatted file.

    Returns:
       Simulation time and potential energy units.
    """

    text = []
    read = True
    while read:  # the loop reads all lines which have # as a first word
        position = filedescU.tell()
        line = filedescU.readline()
        if line == "":
            raise EOFError("The file descriptor hit EOF.")
        elif line.split()[0] == "#":
            text.append(line)
        else:
            filedescU.seek(position)
            read = False

    timeUnit, potentialEnergyUnit = None, None
    potentialEnergy_index, time_index = 0, 0

    potential_re = re.compile(r"potential\{[a-z]*\}")
    time_re = re.compile(r"time\{[a-z]*\}")

    line_index = 0
    for line in text:
        pot = potential_re.search(line)
        time = time_re.search(line)
        if pot is not None:
            potentialEnergyUnit = pot.group()[:-1].split("{")[1]
            potentialEnergy_index = line_index
        if time is not None:
            timeUnit = time.group()[:-1].split("{")[1]
            time_index = line_index
        line_index += 1

    if timeUnit is None or potentialEnergyUnit is None:
        print("Cannot read time and potential energy units")
        sys.exit(1)

    return timeUnit, potentialEnergyUnit, potentialEnergy_index, time_index


def read_U(filedesc, potentialEnergyUnit, potentialEnergy_index, time_index):
    """Takes a file which contains simulation time and potential energy information and
       returns these data. Potential energy is transformed into internal units.

    Args:
       filedesc: An open readable file object.

    Returns:
       The simulation time and potential energy of the system.
    """

    line = filedesc.readline()
    if line == "":
        raise EOFError("The file descriptor hit EOF.")

    line = line.strip().split()
    time = float(line[time_index])
    U = float(line[potentialEnergy_index])
    U = unit_to_internal("energy", potentialEnergyUnit, U)

    return U, time


def main(*arg):
    totalEnergy(*arg)


if __name__ == "__main__":
    main(*sys.argv[1:])

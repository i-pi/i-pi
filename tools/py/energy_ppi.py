#!/usr/bin/env python3

__author__ = "Igor Poltavsky"
__version__ = "1.0"

""" energy_ppi.py
Reads simulation time, potential energy, positions and forces from
an i-PI run and computes a virial total energy estimator and a ppi correction
for each time frame. The output is saved to 'prefix.energy.dat' file which is created in the
folder which contains the input files. The results are printed out in the format: "time
frame", "virial total energy estimator", and "ppi correction".

The script assumes that the input files are in 'xyz' format, with prefix.out (contains simulation time and
potential energy among other output properties), prefix.pos_*.xyz (positions) and
prefix.for_*.xyz (forces) naming scheme.
This would require the following lines in input.xml file:
<properties filename='out' stride='1'> [step, time{picosecond}, potential{kelvin}] </properties>
<trajectory filename='pos' stride='1' format='xyz' cell_units='angstrom'> positions{angstrom} </trajectory>
<trajectory filename='force' stride='1' format='xyz' cell_units='angstrom'> forces{piconewton} </trajectory>

Syntax:
   python energy_ppi.py "prefix" "simulation temperature (in Kelvin)" "number of time frames to skip
   in the beginning of each file (default 0)"
"""

import numpy as np
import sys
import glob
from ipi.utils.messages import verbosity

from ipi.utils.io import read_file
from ipi.utils.units import unit_to_internal, unit_to_user, Constants

verbosity.low = "low"
time_index, potentialEnergy_index = (
    0,
    0,
)  # global variables for time step and potential energy units
potentialEnergyUnit = None  # potential energy unit in input file prefix.out
temperature = None  # simulation temperature
skipSteps = 0  # steps to skip for thermalization


def totalEnergy(prefix, temp, ss=0):
    """
    Computes the virial centroid estimator for the total energy and PPI correction.
    """

    global temperature, skipSteps

    temperature = unit_to_internal("temperature", "kelvin", float(temp))
    skipSteps = int(ss)

    f2_av, ePA_av, eVir_av, f2ePA_av = (
        0.0,
        0.0,
        0.0,
        0.0,
    )  # average square forces, virial and primitive energy
    # estimators, and the product of the primitive energy estimator and square forces

    ipos = []  # input positions files
    for filename in sorted(glob.glob(prefix + ".pos*")):
        ipos.append(open(filename, "r"))

    ifor = []  # input forces files
    for filename in sorted(glob.glob(prefix + ".for*")):
        ifor.append(open(filename, "r"))

    iU = None  # input potential energy and simulation time file
    for filename in sorted(glob.glob(prefix + ".out")):
        iU = open(filename, "r")

    global potentialEnergyUnit, timeUnit
    timeUnit, potentialEnergyUnit = extractUnits(
        iU
    )  # extracting simulation time and potential energy units

    iE = open(prefix + ".energy" + ".dat", "w")
    iE.write(
        "# Simulation time (in %s), virial total energy and PPI energy correction (in %s)\n"
        % (timeUnit, potentialEnergyUnit)
    )

    nbeads = len(ipos)
    if nbeads != len(ifor):
        raise ValueError(
            "Mismatch between number of output files for forces and positions"
        )
    natoms = 0
    ifr = 0
    time0 = 0
    q, f, m = None, None, None
    while True:  # Reading input files and calculating PPI correction
        try:
            for i in range(nbeads):
                ret = read_file("xyz", ipos[i])
                m, n = ret["masses"], ret["atoms"].natoms
                pos = unit_to_internal(ret["units"][0], ret["units"][1], ret["atoms"].q)
                ret = read_file("xyz", ifor[i])
                force = unit_to_internal(
                    ret["units"][0], ret["units"][1], ret["atoms"].q
                )
                if natoms == 0:
                    natoms = n
                    q = np.zeros((nbeads, 3 * natoms))
                    f = np.zeros((nbeads, 3 * natoms))
                q[i, :] = pos
                f[i, :] = force
            time, U = read_U(iU)
        except EOFError:  # finished reading files
            sys.exit(0)

        if ifr < skipSteps:
            time0 = time

        if ifr >= skipSteps:  # PPI correction
            time -= time0

            ePA, f2, f2ePA = 0.0, 0.0, 0.0
            eVir, rc = 0.0, np.zeros(3)

            for j in range(nbeads):
                for i in range(natoms):
                    f2 += (
                        np.dot(f[j, i * 3 : i * 3 + 3], f[j, i * 3 : i * 3 + 3]) / m[i]
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
            for i in range(natoms):
                rc[:] = 0.0
                for j in range(nbeads):
                    rc[:] += q[j, i * 3 : i * 3 + 3]
                rc[:] /= nbeads
                for j in range(nbeads):
                    eVir += np.dot(
                        rc[:] - q[j, i * 3 : i * 3 + 3], f[j, i * 3 : i * 3 + 3]
                    )

            ePA *= 0.5 * nbeads * (Constants.kb * temperature) ** 2 / Constants.hbar**2
            ePA += 0.5 * nbeads * (3 * natoms) * Constants.kb * temperature + U
            f2ePA = f2 * ePA
            eVir /= 2.0 * nbeads
            eVir += 0.5 * (3 * natoms) * Constants.kb * temperature + U

            ePA_av += ePA
            f2_av += f2
            f2ePA_av += f2ePA
            eVir_av += eVir
            ifr += 1
            print(
                (ifr - skipSteps)
            )  # Printing current time frame (excluding thermalization)

            dE = (
                3.0 * Constants.kb * temperature + ePA_av / float(ifr - skipSteps)
            ) * f2_av / float(ifr - skipSteps) - f2ePA_av / float(ifr - skipSteps)
            dE *= Constants.hbar**2 / (
                24.0 * (nbeads * Constants.kb * temperature) ** 3
            )

            dE = unit_to_user(
                "energy", potentialEnergyUnit, dE
            )  # Output in the same unit as potential energy
            eVir = unit_to_user(
                "energy", potentialEnergyUnit, eVir_av / float(ifr - skipSteps)
            )  # Output in the same unit
            # as potential energy

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

    global potentialEnergy_index, time_index

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

    line_index = 0
    for line in text:
        ind = 0
        try:
            ind = line.find("potential{")
            if ind != -1:
                potentialEnergy_index = int(line_index)
                line = line[ind + 10 :]
                ind = line.find("}")
                unit = line[:ind]
                potentialEnergyUnit = unit
        except ValueError:
            pass
        try:
            ind = line.find("time{")
            if ind != -1:
                time_index = int(line_index)
                line = line[ind + 5 :]
                ind = line.find("}")
                unit = line[:ind]
                timeUnit = unit
        except ValueError:
            pass
        line_index += 1

    if timeUnit is None or potentialEnergyUnit is None:
        print("Cannot read time and potential energy units")
        sys.exit(1)

    return timeUnit, potentialEnergyUnit


def read_U(filedesc):
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

    line = line.strip()
    line = line.split()
    time = float(line[time_index])
    U = float(line[potentialEnergy_index])
    U = unit_to_internal("energy", potentialEnergyUnit, U)

    return time, U


def main(*arg):
    totalEnergy(*arg)


if __name__ == "__main__":
    main(*sys.argv[1:])

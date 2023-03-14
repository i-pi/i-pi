#!/usr/local/bin/python

__author__ = "Igor Poltavsky"
__version__ = "1.0"

""" rdf_ppi.py
The script reads simulation time, potential energy, positions and forces from
standard i-PI output files and computes a conventional and PPI radial distribution function (RDF) estimators.
The output is saved to two files which are created in the folder which contains the input files.
The results are printed out in the format: "distance", "RDF".

The script assumes that the input files are in 'xyz' format, with prefix.pos_*.xyz (positions) and
prefix.for_*.xyz (forces) naming scheme.
This would require the following lines in input.xml file:
<trajectory filename='pos' stride='n' format='xyz' cell_units='angstrom'> positions </trajectory>
<trajectory filename='for' stride='n' format='xyz' cell_units='angstrom'> forces </trajectory>
where n is the same integer number.

Syntax:
   python rdf_ppi.py "prefix" "simulation temperature (in Kelvin)" "element A" "element B" "number of bins for RDF"
   "minimum distance (in Angstroms)" "maximum distance (in Angstroms) "number of time frames to skip in the beginning
   of each file (default 0)"

WARNING:
   Since the python code for computing RDF is inefficient the fortran function in f90 folder must be compiled to
   compute RDFs.
"""

import numpy as np
import sys
import glob
import os
from ipi.utils.units import unit_to_internal, unit_to_user, Constants, Elements
from ipi.utils.io import read_file


def RDF(prefix, temp, A, B, nbins, r_min, r_max, ss=0, unit="angstrom"):
    # Adding fortran functions (when exist)
    sys.path.append(os.path.abspath(os.path.dirname(sys.argv[0]))[:-2] + "f90")
    try:
        import fortran
    except ImportError:
        print(
            "WARNING: No compiled fortran module for fast calculations have been found.\n"
            "Proceeding the calculations is not possible."
        )
        sys.exit(0)

    temperature = unit_to_internal(
        "temperature", "kelvin", float(temp)
    )  # simulation temperature
    skipSteps = int(ss)  # steps to skip
    nbins = int(nbins)

    fns_pos = sorted(glob.glob(prefix + ".pos*"))
    fns_for = sorted(glob.glob(prefix + ".for*"))
    fn_out_rdf, fn_out_rdf_q = (
        prefix + "." + A + B + ".rdf.dat",
        prefix + "." + A + B + ".rdf-ppi.dat",
    )

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
    print("output file names:")
    print(fn_out_rdf)
    print(fn_out_rdf_q)
    print()

    # open input and output files
    ipos = [open(fn, "r") for fn in fns_pos]
    ifor = [open(fn, "r") for fn in fns_for]
    # iRDF, iRDFq = open(fn_out_rdf, "w"), open(fn_out_rdf_q, "w")

    # Species for RDF
    species = (A, B)
    speciesMass = np.array(
        [Elements.mass(species[0]), Elements.mass(species[1])], order="F"
    )

    r_min = unit_to_internal("length", unit, float(r_min))  # Minimal distance for RDF
    r_max = unit_to_internal("length", unit, float(r_max))  # Maximal distance for RDF
    dr = (r_max - r_min) / nbins  # RDF step
    rdf = np.array(
        [[r_min + (0.5 + i) * dr, 0] for i in range(nbins)], order="F"
    )  # conventional RDF

    # RDF auxiliary variables
    cell = None  # simulation cell matrix
    inverseCell = None  # inverse simulation sell matrix
    cellVolume = None  # simulation cell volume
    natomsA = 0  # the total number of A type particles in the system
    natomsB = 0  # the total number of B type particles in the system
    # Here A and B are the types of elements used for RDF calculations

    # temp variables
    f2 = 0.0  # the sum of square forces divided by corresponding masses (f2/m)
    f2rdf = np.zeros(
        nbins, order="F"
    )  # the sum f2/m multiplied by the rdf at each time step
    frdf = np.zeros(
        nbins, order="F"
    )  # the sum of f/m multiplied by the rdf derivative at each time step

    natoms = 0  # total number of atoms
    ifr = 0  # time frame number
    pos, force, mass = None, None, None  # positions, forces, and mass arrays
    noteof = True  # end of file test variable

    while noteof:  # Reading input files and calculating PPI correction
        if ifr % 100 == 0:
            print("\rProcessing frame {:d}".format(ifr), end=" ")
            sys.stdout.flush()

        try:
            for i in range(nbeads):
                ret = read_file("xyz", ipos[i], dimension="length")
                if natoms == 0:
                    mass, natoms = ret["atoms"].m, ret["atoms"].natoms
                    pos = np.zeros((nbeads, 3 * natoms), order="F")
                    force = np.zeros((nbeads, 3 * natoms), order="F")
                cell = ret["cell"].h
                inverseCell = ret["cell"].get_ih()
                cellVolume = ret["cell"].get_volume()
                pos[i, :] = ret["atoms"].q
                force[i, :] = read_file("xyz", ifor[i], dimension="force")["atoms"].q
        except EOFError:  # finished reading files
            noteof = False

        if noteof:
            if ifr >= skipSteps:  # RDF calculations
                species_A = [
                    3 * i + j
                    for i in np.where(mass == speciesMass[0])[0]
                    for j in range(3)
                ]
                species_B = [
                    3 * i + j
                    for i in np.where(mass == speciesMass[1])[0]
                    for j in range(3)
                ]
                natomsA = len(species_A)
                natomsB = len(species_B)
                posA = np.zeros((nbeads, natomsA), order="F")
                posB = np.zeros((nbeads, natomsB), order="F")
                forA = np.zeros((nbeads, natomsA), order="F")
                forB = np.zeros((nbeads, natomsB), order="F")
                for bead in range(nbeads):
                    posA[bead, :] = pos[bead, species_A]
                    forA[bead, :] = force[bead, species_A]
                    posB[bead, :] = pos[bead, species_B]
                    forB[bead, :] = force[bead, species_B]

                # RDF amd PPI RDF calculations
                f2temp = fortran.f2divm(force, mass, natoms, nbeads)
                f2 += f2temp
                fortran.updateqrdf(
                    rdf,
                    f2rdf,
                    frdf,
                    posA,
                    posB,
                    forA,
                    forB,
                    natomsA / 3,
                    natomsB / 3,
                    nbins,
                    r_min,
                    r_max,
                    cell,
                    inverseCell,
                    nbeads,
                    f2temp,
                    speciesMass[0],
                    speciesMass[1],
                )

                ifr += 1

            else:
                ifr += 1

            if ifr > skipSteps and ifr % 100 == 0:
                # Some constants
                const = 1.0 / float(ifr - skipSteps)
                alpha = Constants.hbar**2 / (
                    24.0 * nbeads**3 * (temperature * Constants.kb) ** 3
                )
                beta = Constants.hbar**2 / (
                    12.0 * nbeads**3 * (temperature * Constants.kb) ** 2
                )

                # Normalization
                _rdf = np.copy(rdf)
                _f2rdf = np.copy(f2rdf)
                _frdf = np.copy(frdf)
                _rdf[:, 1] *= const / nbeads
                _f2 = f2 * alpha * const
                _f2rdf[:] *= alpha * const
                _frdf[:] *= beta * const

                # PPI correction
                rdfQ = np.copy(_rdf)
                for bin in range(nbins):
                    rdfQ[bin, 1] += _rdf[bin, 1] * _f2 - _f2rdf[bin] / nbeads
                    rdfQ[bin, 1] -= _frdf[bin] / 2.0

                # Creating RDF from N(r)
                const, dr = cellVolume / (4 * np.pi / 3.0), _rdf[1, 0] - _rdf[0, 0]
                for bin in range(nbins):
                    _rdf[bin, 1] = (
                        const
                        * _rdf[bin, 1]
                        / (
                            (_rdf[bin, 0] + 0.5 * dr) ** 3
                            - (_rdf[bin, 0] - 0.5 * dr) ** 3
                        )
                    )
                    rdfQ[bin, 1] = (
                        const
                        * rdfQ[bin, 1]
                        / (
                            (rdfQ[bin, 0] + 0.5 * dr) ** 3
                            - (rdfQ[bin, 0] - 0.5 * dr) ** 3
                        )
                    )
                for bin in range(nbins):
                    _rdf[bin, 0] = unit_to_user("length", unit, _rdf[bin, 0])
                    rdfQ[bin, 0] = unit_to_user("length", unit, rdfQ[bin, 0])

                # Writing the results into files
                np.savetxt(fn_out_rdf, _rdf)
                np.savetxt(fn_out_rdf_q, rdfQ)


def main(*arg):
    RDF(*arg)


if __name__ == "__main__":
    main(*sys.argv[1:])

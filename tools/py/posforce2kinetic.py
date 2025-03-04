#!/usr/bin/env python3

"""posforce2kinetic.py

Reads positions and forces from an i-PI run and computes the centroid-virial
kinetic energy estimator for each particle and time frame, printing it out in
the same format as it would have been obtained had the run included the
kinetic_cv and kinetic_od output trajectories.

Assumes the input files are in xyz format and atomic units, with
prefix.pos_*.xyz and prefix.for_*.xyz naming scheme.

Syntax:
   posforce2kinetic.py prefix temperature[K]
"""


import numpy as np
import sys
import glob
from ipi.utils.io import read_file
from ipi.engine.beads import Beads
from ipi.utils.depend import dstrip
from ipi.utils.units import unit_to_internal, Constants
from ipi.utils.messages import verbosity

verbosity.level = "low"


def main(prefix, temp):
    T = float(temp)
    fns_pos = sorted(glob.glob(prefix + ".pos*"))
    fns_for = sorted(glob.glob(prefix + ".for*"))
    fn_out_kin = prefix + ".kin.xyz"
    fn_out_kod = prefix + ".kod.xyz"

    # check that we found the same number of positions and forces files
    nbeads = len(fns_pos)
    if nbeads != len(fns_for):
        print(fns_pos)
        print(fns_for)
        raise ValueError(
            "Mismatch between number of input files for forces and positions."
        )

    # print some information
    print("temperature = {:f} K".format(T))
    print()
    print("number of beads = {:d}".format(nbeads))
    print()
    print("positions and forces file names:")
    for fn_pos, fn_for in zip(fns_pos, fns_for):
        print("{:s}   {:s}".format(fn_pos, fn_for))
    print()
    print("output file names:")
    print(fn_out_kin)
    print(fn_out_kod)
    print()

    temp = unit_to_internal("energy", "kelvin", T)

    # open input and output files
    ipos = [open(fn, "r") for fn in fns_pos]
    ifor = [open(fn, "r") for fn in fns_for]
    ikin = open(fn_out_kin, "w")
    ikod = open(fn_out_kod, "w")

    natoms = 0
    ifr = 0
    while True:
        # print progress
        if ifr % 100 == 0:
            print("\rProcessing frame {:d}".format(ifr), end=" ")
            sys.stdout.flush()

        # load one frame
        try:
            for i in range(nbeads):
                ret = read_file("xyz", ipos[i])
                pos = ret["atoms"]
                ret = read_file("xyz", ifor[i])
                force = ret["atoms"]
                if natoms == 0:
                    natoms = pos.natoms
                    beads = Beads(natoms, nbeads)
                    forces = Beads(natoms, nbeads)
                    kcv = np.zeros((natoms, 6), float)
                beads[i].q = pos.q
                forces[i].q = force.q
        except EOFError:
            # finished reading files
            break

        # calculate kinetic energies
        q = dstrip(beads.q)
        f = dstrip(forces.q)
        qc = dstrip(beads.qc)
        kcv[:] = 0
        for j in range(nbeads):
            for i in range(natoms):
                kcv[i, 0] += (q[j, i * 3 + 0] - qc[i * 3 + 0]) * f[j, i * 3 + 0]
                kcv[i, 1] += (q[j, i * 3 + 1] - qc[i * 3 + 1]) * f[j, i * 3 + 1]
                kcv[i, 2] += (q[j, i * 3 + 2] - qc[i * 3 + 2]) * f[j, i * 3 + 2]
                kcv[i, 3] += (q[j, i * 3 + 0] - qc[i * 3 + 0]) * f[j, i * 3 + 1] + (
                    q[j, i * 3 + 1] - qc[i * 3 + 1]
                ) * f[j, i * 3 + 0]
                kcv[i, 4] += (q[j, i * 3 + 0] - qc[i * 3 + 0]) * f[j, i * 3 + 2] + (
                    q[j, i * 3 + 2] - qc[i * 3 + 2]
                ) * f[j, i * 3 + 0]
                kcv[i, 5] += (q[j, i * 3 + 1] - qc[i * 3 + 1]) * f[j, i * 3 + 2] + (
                    q[j, i * 3 + 2] - qc[i * 3 + 2]
                ) * f[j, i * 3 + 1]
        kcv *= -0.5 / nbeads
        kcv[:, 0:3] += 0.5 * Constants.kb * temp
        kcv[:, 3:6] *= 0.5

        # write output
        ikin.write(
            "%d\n# Centroid-virial kinetic energy estimator [a.u.] - diagonal terms: xx yy zz\n"
            % natoms
        )
        ikod.write(
            "%d\n# Centroid-virial kinetic energy estimator [a.u.] - off-diag terms: xy xz yz\n"
            % natoms
        )
        for i in range(natoms):
            ikin.write(
                "%8s %12.5e %12.5e %12.5e\n"
                % (pos.names[i], kcv[i, 0], kcv[i, 1], kcv[i, 2])
            )
            ikod.write(
                "%8s %12.5e %12.5e %12.5e\n"
                % (pos.names[i], kcv[i, 3], kcv[i, 4], kcv[i, 5])
            )

        ifr += 1

    print("\rProcessed {:d} frames.".format(ifr))

    ikin.close()
    ikod.close()


if __name__ == "__main__":
    main(*sys.argv[1:])

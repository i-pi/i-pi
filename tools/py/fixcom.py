#!/usr/bin/env python3

"""pos2centroid.py

Reads positions of individual beads from an i-PI run and
assemles them in a pdb describing the ring polymer connectivity.

Assumes the input files are in pdb format names prefix.pos_*.pdb.

Syntax:
   mergebeadspdb.py prefix
"""


import numpy as np
import sys
import glob
from ipi.utils.io import *
from ipi.engine.atoms import Atoms
from ipi.utils.depend import *
from ipi.utils.units import *


def main(prefix, suffix="pos", outprefix="fixcom"):
    ipos = []
    imode = []
    for filename in sorted(glob.glob(prefix + "." + suffix + "*")):
        imode.append(filename.split(".")[-1])
        ipos.append(open(filename, "r"))

    nbeads = len(ipos)
    natoms = 0
    ifr = 0

    lout = []
    for b in range(nbeads):
        # zero-padded bead number
        padb = ("%0" + str(int(1 + np.floor(np.log(nbeads) / np.log(10)))) + "d") % (b)
        lout.append(
            open(
                prefix + "." + outprefix + "." + suffix + "_" + padb + "." + imode[b],
                "a",
            )
        )

    while True:
        allbeads = []
        for i in range(nbeads):
            try:
                poscell = read_file(imode[i], ipos[i])
                cell = poscell["cell"]
                pos = poscell["atoms"]
                allbeads.append(pos)
                if natoms == 0:
                    natoms = pos.natoms
                    atoms = Atoms(natoms)

                atoms.q += pos.q
                atoms.names = pos.names
                atoms.m += pos.m
            except EOFError:  # finished reading files
                for ib in range(nbeads):
                    lout[ib].close()
                sys.exit(0)
        atoms.q /= nbeads
        atoms.m /= nbeads
        com = np.zeros(3)
        tm = 0
        for i in range(atoms.natoms):
            com += atoms.m[i] * atoms.q[3 * i : 3 * (i + 1)]
            tm += atoms.m[i]
        com /= tm
        for ib in range(nbeads):
            for i in range(allbeads[ib].natoms):
                allbeads[ib].q[3 * i : 3 * (i + 1)] -= com
            print_file(imode[ib], allbeads[ib], cell, filedesc=lout[ib])
        atoms.q[:] = 0.0
        atoms.m[:] = 0.0
        ifr += 1


if __name__ == "__main__":
    main(*sys.argv[1:])

#!/usr/bin/env python3

"""mergebeadspdb.py

Reads positions of individual beads from an i-PI run and
assemles them in a pdb describing the ring polymer connectivity.

Assumes the input files are in named following the convention prefix.pos_*.*.

Syntax:
   mergebeadspdb.py prefix
"""


import sys
import glob
from ipi.utils.io import read_file, print_file_path
from ipi.engine.beads import Beads
from ipi.utils.messages import verbosity

verbosity.level = "low"


def main(prefix, suffix="pos", unitconv="1.0"):
    ipos = []
    imode = []
    for filename in sorted(glob.glob(prefix + "." + suffix + "*")):
        imode.append(filename.split(".")[-1])
        ipos.append(open(filename, "r"))

    nbeads = len(ipos)
    natoms = 0
    ifr = 0
    while True:
        try:
            for i in range(nbeads):
                ret = read_file(imode[i], ipos[i])
                pos = ret["atoms"]
                cell = ret["cell"]
                if natoms == 0:
                    natoms = pos.natoms
                    beads = Beads(natoms, nbeads)
                cell.h *= float(unitconv)
                beads[i].q = pos.q * float(unitconv)
                beads.names = pos.names
        except EOFError:  # finished reading files
            sys.exit(0)

        print_file_path("pdb", beads, cell)
        ifr += 1


if __name__ == "__main__":
    main(*sys.argv[1:])

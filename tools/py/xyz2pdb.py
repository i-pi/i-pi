#!/usr/bin/env python3

"""xyz2pdb.py

Reads positions of a system in xyz format and returns a
pdb file. Also reads the comment line and returns the
positions and the cell in Angstrom.

Assumes the input file is in a pdb format.

Syntax:
   xyz2pdb.py filename

"""


import sys
from ipi.utils.io import read_file, print_file


def main(filename, wrap=True):
    ipos = open(filename, "r")
    if wrap == "False":
        wrap = False

    ifr = 0
    while True:
        try:
            ret = read_file("xyz", ipos)
            pos = ret["atoms"]
            cell = ret["cell"]
            if wrap:
                cell.array_pbc(pos.q)
        except EOFError:  # finished reading files
            sys.exit(0)

        print_file("pdb", pos, cell)
        ifr += 1


if __name__ == "__main__":
    main(*sys.argv[1:])

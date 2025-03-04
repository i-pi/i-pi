#!/usr/bin/env python3

"""xyz2bin.py

Reads positions of a system in xyz format and returns a
binary file.

Assumes the input file is in a xyz format.

Syntax:
   xyz2bin.py filename

"""


import sys
from ipi.utils.io import read_file, print_file


def main(filename):
    ipos = open(filename, "r")

    ifr = 0
    while True:
        try:
            ret = read_file("xyz", ipos)
            pos = ret["atoms"]
            cell = ret["cell"]
            cell.array_pbc(pos.q)
        except EOFError:  # finished reading files
            sys.exit(0)

        print_file("bin", pos, cell)
        ifr += 1


if __name__ == "__main__":
    main(*sys.argv[1:])

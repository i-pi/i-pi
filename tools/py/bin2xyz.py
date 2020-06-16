#!/usr/bin/env python2

""" bin2xyz.py

Reads positions of a system in binary format and returns a
xyz file.

Assumes the input file is in a binary format.

Syntax:
   bin2xyz.py filename
"""


import sys
from ipi.utils.io import read_file, print_file
from ipi.utils.depend import *
from ipi.utils.units import *


def main(filename):

    ipos = open(filename, "rb")

    natoms = 0
    ifr = 0
    while True:
        try:
            ret = read_file("bin", ipos)
            pos = ret["atoms"]
            cell = ret["cell"]
            cell.array_pbc(pos.q)
        except EOFError:  # finished reading files
            sys.exit(0)

        print_file("xyz", pos, cell)
        ifr += 1


if __name__ == '__main__':
    main(*sys.argv[1:])

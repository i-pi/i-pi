"""Functions used to read input configurations and print trajectories
in the XYZ format.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import sys

import numpy as np

import ipi.utils.mathtools as mt
from ipi.utils.depend import dstrip

__all__ = ["print_txt"]


def print_txt(
    atoms, cell, filedesc=sys.stdout, title="", cell_conv=1.0, atoms_conv=1.0
):
    """Prints an array into an TXT formatted file.

    Args:
        atoms: An atoms object giving the centroid positions.
        cell: A cell object giving the system box.
        filedesc: An open writable file object. Defaults to standard output.
        title: This gives a string to be appended to the comment line.
    """

    a, b, c, alpha, beta, gamma = mt.h2abc_deg(cell.h * cell_conv)

    natoms = atoms.natoms

    string = "# n. rows: {:d}\n".format(natoms)
    filedesc.write(string)
    data = dstrip(atoms.q) * atoms_conv
    np.savetxt(filedesc, data)
    filedesc.write("\n")
    pass

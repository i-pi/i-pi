"""Functions used to read input configurations and print trajectories
in the XYZ format.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import sys

__all__ = ["print_text"]



def print_text(
    atoms, cell, filedesc=sys.stdout, title="", cell_conv=1.0, atoms_conv=1.0
):
    """Prints an array into an TXT formatted file.

    Args:
        atoms: An atoms object giving the centroid positions.
        cell: A cell object giving the system box.
        filedesc: An open writable file object. Defaults to standard output.
        title: This gives a string to be appended to the comment line.
    """

    string = atoms + "\n"
    filedesc.write(string)

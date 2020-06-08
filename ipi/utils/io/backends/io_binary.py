"""Functions used to print the trajectories and read input configurations
(or even full status dump) as unformatted binary.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import sys
import numpy as np


__all__ = ["print_binary", "read_binary"]


def print_binary(
    atoms, cell, filedesc=sys.stdout, title="", cell_conv=1.0, atoms_conv=1.0
):
    """Prints an atomic configuration into a binary file.

    Args:
        beads: An atoms object giving the centroid positions.
        cell: A cell object giving the system box.
        filedesc: An open writable file object. Defaults to standard output.
        title: This gives a string to be appended to the comment line.
    """

    buff = filedesc
    # applies conversion of units to the cell and saves to file.
    (cell.h.flatten() * cell_conv).tofile(buff)
    nat = np.asarray([atoms.natoms])
    nat.tofile(buff)
    # applies conversion of units to the atoms and saves to file.
    (atoms.q * atoms_conv).tofile(buff)
    # concatenates names, assuming none contains '|'
    names = "|".join(atoms.names)
    nat[0] = len(names)
    nat.tofile(buff)
    np.asarray([names]).tofile(buff)
    # also saves the title to the file.
    np.asarray([title]).tofile(buff)


def read_binary(filedesc):
    try:
        cell = np.fromfile(filedesc, dtype=float, count=9)
        cell.shape = (3, 3)
        nat = np.fromfile(filedesc, dtype=int, count=1)[0]
        qatoms = np.fromfile(filedesc, dtype=float, count=3 * nat)
        nat = np.fromfile(filedesc, dtype=int, count=1)[0]
        names = np.fromfile(filedesc, dtype="|U1", count=nat)
        names = "".join(names)
        names = names.split("|")
        masses = np.zeros(len(names))
        title = "".join(np.fromfile(filedesc, dtype="|U1", count=-1))
    except (StopIteration, ValueError):
        raise EOFError
    return (title, cell, qatoms, names, masses)

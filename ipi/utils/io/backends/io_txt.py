"""Functions used to read input configurations and print trajectories
in the XYZ format.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import sys
import re

import numpy as np

import ipi.utils.mathtools as mt
from ipi.utils.depend import dstrip
from ipi.utils.units import Elements
from ipi import ipi_global_settings


__all__ = ["print_txt_path", "print_xyz"]

deg2rad = np.pi / 180.0


def print_txt_path(beads, cell, filedesc=sys.stdout, cell_conv=1.0, atoms_conv=1.0):
    """Prints all the bead configurations into a TXT formatted file.

    Prints all the replicas for each time step separately, rather than all at
    once.

    Args:
        beads: A beads object giving the bead positions.
        cell: A cell object giving the system box.
        filedesc: An open writable file object. Defaults to standard output.
    """

    a, b, c, alpha, beta, gamma = mt.h2abc_deg(cell.h * cell_conv)

    fmt_header = (
        "%d\n# bead: %d CELL(abcABC): %10.5f  %10.5f  %10.5f  %10.5f  %10.5f  %10.5f \n"
    )
    natoms = beads.natoms
    nbeads = beads.nbeads
    for j in range(nbeads):
        filedesc.write(fmt_header % (natoms, j, a, b, c, alpha, beta, gamma))
        for i in range(natoms):
            qs = dstrip(beads.q) * atoms_conv
            lab = dstrip(beads.names)
            filedesc.write(
                "%8s %12.5e %12.5e %12.5e\n"
                % (lab[i], qs[j][3 * i], qs[j][3 * i + 1], qs[j][3 * i + 2])
            )


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
    # fmt_header = (
    #     "%d\n# CELL(abcABC): %10.5f  %10.5f  %10.5f  %10.5f  %10.5f  %10.5f  %s\n"
    # )
    # filedesc.write(fmt_header % (natoms, a, b, c, alpha, beta, gamma, title))
    # direct access to avoid unnecessary slow-down

    string = "# n. rows: {:d}\n".format(natoms)
    filedesc.write(string)
    data = dstrip(atoms.q) * atoms_conv
    np.savetxt(filedesc, data)
    filedesc.write("\n")
    pass

    # qs = dstrip(atoms.q) * atoms_conv
    # lab = dstrip(atoms.names)
    # string="# n. rows: {:d}\n\n".format(natoms)
    # for i in range(natoms):
    #     for j in range(qs.shape[1]):
    #         string += ipi_global_settings["floatformat"]%(qs[i,j])
    #     string += "\n"
    # string += "\n"
    # filedesc.write(string)

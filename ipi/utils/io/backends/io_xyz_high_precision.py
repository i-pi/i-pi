"""Read and write XYZ trajectories with round-trip double precision."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import sys

import ipi.utils.mathtools as mt
from ipi.utils.depend import dstrip
from ipi.utils.io.backends.io_xyz import read_xyz as read_xyz_high_precision

__all__ = [
    "print_xyz_high_precision_path",
    "print_xyz_high_precision",
    "read_xyz_high_precision",
]


def _write_atoms(names, coordinates, filedesc):
    for i, name in enumerate(names):
        offset = 3 * i
        filedesc.write(
            "%8s %.16e %.16e %.16e\n"
            % (
                name,
                coordinates[offset],
                coordinates[offset + 1],
                coordinates[offset + 2],
            )
        )


def print_xyz_high_precision_path(
    beads, cell, filedesc=sys.stdout, title="", cell_conv=1.0, atoms_conv=1.0
):
    """Print every bead as an XYZ frame with round-trip double precision."""

    a, b, c, alpha, beta, gamma = mt.h2abc_deg(cell.h * cell_conv)
    header = "%d\n# bead: %d CELL(abcABC): %.16e %.16e %.16e " "%.16e %.16e %.16e %s\n"
    coordinates = dstrip(beads.q) * atoms_conv
    names = dstrip(beads.names)
    for bead in range(beads.nbeads):
        filedesc.write(
            header
            % (
                beads.natoms,
                bead,
                a,
                b,
                c,
                alpha,
                beta,
                gamma,
                title,
            )
        )
        _write_atoms(names, coordinates[bead], filedesc)


def print_xyz_high_precision(
    atoms, cell, filedesc=sys.stdout, title="", cell_conv=1.0, atoms_conv=1.0
):
    """Print one XYZ frame with round-trip double precision."""

    a, b, c, alpha, beta, gamma = mt.h2abc_deg(cell.h * cell_conv)
    header = "%d\n# CELL(abcABC): %.16e %.16e %.16e %.16e %.16e %.16e %s\n"
    filedesc.write(header % (atoms.natoms, a, b, c, alpha, beta, gamma, title))
    _write_atoms(dstrip(atoms.names), dstrip(atoms.q) * atoms_conv, filedesc)

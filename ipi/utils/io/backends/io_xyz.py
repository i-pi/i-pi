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
from ipi.utils.array_backend import to_numpy
from ipi.utils.depend import dstrip
from ipi.utils.units import Elements

__all__ = ["print_xyz_path", "print_xyz", "read_xyz"]

deg2rad = np.pi / 180.0


def print_xyz_path(beads, cell, filedesc=sys.stdout, cell_conv=1.0, atoms_conv=1.0):
    """Prints all the bead configurations into a XYZ formatted file.

    Prints all the replicas for each time step separately, rather than all at
    once.

    Args:
        beads: A beads object giving the bead positions.
        cell: A cell object giving the system box.
        filedesc: An open writable file object. Defaults to standard output.
    """

    # I/O boundary: pull into numpy once, reshape to (nbeads, natoms, 3)
    # so the inner loop is a plain row access.
    h = to_numpy(cell.h) * cell_conv
    a, b, c, alpha, beta, gamma = mt.h2abc_deg(h)

    natoms = beads.natoms
    nbeads = beads.nbeads
    qs_all = (to_numpy(beads.q) * atoms_conv).reshape(nbeads, natoms, 3)
    lab = dstrip(beads.names)
    for j in range(nbeads):
        filedesc.write(
            f"{natoms}\n# bead: {j} CELL(abcABC): "
            f"{a:10.5f}  {b:10.5f}  {c:10.5f}  "
            f"{alpha:10.5f}  {beta:10.5f}  {gamma:10.5f} \n"
        )
        qs = qs_all[j]
        for i in range(natoms):
            filedesc.write(
                f"{lab[i]:>8s} {qs[i, 0]:12.5e} {qs[i, 1]:12.5e} {qs[i, 2]:12.5e}\n"
            )


def print_xyz(
    atoms, cell, filedesc=sys.stdout, title="", cell_conv=1.0, atoms_conv=1.0
):
    """Prints an atomic configuration into an XYZ formatted file.

    Args:
        atoms: An atoms object giving the centroid positions.
        cell: A cell object giving the system box.
        filedesc: An open writable file object. Defaults to standard output.
        title: This gives a string to be appended to the comment line.
    """

    # I/O boundary: pull into numpy once, reshape to (natoms, 3) so the
    # inner loop is a plain row access.
    h = to_numpy(cell.h) * cell_conv
    a, b, c, alpha, beta, gamma = mt.h2abc_deg(h)

    natoms = atoms.natoms
    filedesc.write(
        f"{natoms}\n# CELL(abcABC): "
        f"{a:10.5f}  {b:10.5f}  {c:10.5f}  "
        f"{alpha:10.5f}  {beta:10.5f}  {gamma:10.5f}  {title}\n"
    )
    qs = (to_numpy(atoms.q) * atoms_conv).reshape(natoms, 3)
    lab = dstrip(atoms.names)
    for i in range(natoms):
        filedesc.write(
            f"{lab[i]:>8s} {qs[i, 0]:12.5e} {qs[i, 1]:12.5e} {qs[i, 2]:12.5e}\n"
        )


# Cell type patterns
cell_re = [
    re.compile(r"CELL[\(\[\{]abcABC[\)\]\}]: ([-+0-9\.Ee ]*)\s*"),
    re.compile(r"CELL[\(\[\{]H[\)\]\}]: ([-+0-9\.?Ee ]*)\s*"),
    re.compile(r"CELL[\(\[\{]GENH[\)\]\}]: ([-+0-9\.?Ee ]*)\s*"),
]


def read_xyz(filedesc):
    """Reads an XYZ-style file with i-PI style comments and returns data in raw format for further units transformation
    and other post processing.

    Args:
        filedesc: An open readable file object from a xyz formatted file with i-PI header comments.

    Returns:
        i-Pi comment line, cell array, data (positions, forces, etc.), atoms names and masses
    """

    try:
        natoms = int(next(filedesc))
    except (StopIteration, ValueError):
        raise EOFError

    comment = next(filedesc)

    # Extracting cell
    cell = [key.search(comment) for key in cell_re]
    usegenh = False
    if cell[0] is not None:  # abcABC
        a, b, c = [float(x) for x in cell[0].group(1).split()[:3]]
        alpha, beta, gamma = [float(x) * deg2rad for x in cell[0].group(1).split()[3:6]]
        h = mt.abc2h(a, b, c, alpha, beta, gamma)
    elif cell[1] is not None:  # H
        h = np.array(cell[1].group(1).split()[:9], float)
        h.resize((3, 3))
    elif cell[2] is not None:  # GENH
        genh = np.array(cell[2].group(1).split()[:9], float)
        genh.resize((3, 3))
        invgenh = np.linalg.inv(genh)
        # convert back & forth from abcABC representation to get an upper triangular h
        h = mt.abc2h(*mt.genh2abc(genh))
        usegenh = True
    else:  # defaults to unit box
        h = np.array([[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]])
    cell = h

    qatoms = np.zeros(3 * natoms)
    names = np.zeros(natoms, dtype="|U4")
    masses = np.zeros(natoms)

    # Extracting a time-frame information
    atom_counter = 0
    for iat, line in enumerate(filedesc):
        body = line.split()
        names[iat], masses[iat] = body[0], Elements.mass(body[0])
        x, y, z = float(body[1]), float(body[2]), float(body[3])

        if usegenh:
            # must convert from the input cell parameters to the internal convention
            u = np.array([x, y, z])
            us = np.dot(u, invgenh)
            u = np.dot(h, us)
            x, y, z = u

        qatoms[3 * iat], qatoms[3 * iat + 1], qatoms[3 * iat + 2] = x, y, z
        atom_counter += 1
        if atom_counter == natoms:
            break

    if natoms != len(names):
        raise ValueError(
            "The number of atom records does not match the header of the xyz file."
        )
    return comment, cell, qatoms, names, masses

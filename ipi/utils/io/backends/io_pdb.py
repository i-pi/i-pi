"""Functions used to read input configurations and print trajectories
in the PDB format.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import sys
import copy

import numpy as np

import ipi.utils.mathtools as mt
from ipi.utils.depend import dstrip
from ipi.utils.units import Elements


__all__ = ["print_pdb_path", "print_pdb", "read_pdb"]


def print_pdb_path(beads, cell, filedesc=sys.stdout, cell_conv=1.0, atoms_conv=1.0):
    """Prints all the bead configurations, into a pdb formatted file.

    Prints the ring polymer springs as well as the bead positions using the
    CONECT command. Also prints the cell parameters in standard pdb form. Note
    that the angles are in degrees.

    Args:
        beads: A beads object giving the bead positions.
        cell: A cell object giving the system box.
        filedesc: An open writable file object. Defaults to standard output.
    """

    fmt_cryst = "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%s%4i\n"
    fmt_atom = (
        "ATOM  %5i %4s%1s%3s %1s%4i%1s  %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2i\n"
    )
    fmt_conect = "CONECT%5i%5i\n"

    a, b, c, alpha, beta, gamma = mt.h2abc_deg(cell.h * cell_conv)

    z = 1  # What even is this parameter?
    filedesc.write(fmt_cryst % (a, b, c, alpha, beta, gamma, " P 1        ", z))

    natoms = beads.natoms
    nbeads = beads.nbeads
    for j in range(nbeads):
        for i in range(natoms):
            qs = dstrip(beads.q) * atoms_conv
            lab = dstrip(beads.names)
            data = (
                j * natoms + i + 1,
                lab[i],
                " ",
                "  1",
                " ",
                1,
                " ",
                qs[j][3 * i],
                qs[j][3 * i + 1],
                qs[j][3 * i + 2],
                0.0,
                0.0,
                "  ",
                0,
            )
            filedesc.write(fmt_atom % data)

    if nbeads > 1:
        for i in range(natoms):
            filedesc.write(fmt_conect % (i + 1, (nbeads - 1) * natoms + i + 1))
        for j in range(nbeads - 1):
            for i in range(natoms):
                filedesc.write(
                    fmt_conect % (j * natoms + i + 1, (j + 1) * natoms + i + 1)
                )

    filedesc.write("END\n")


def print_pdb(
    atoms, cell, filedesc=sys.stdout, title="", cell_conv=1.0, atoms_conv=1.0
):
    """Prints an atomic configuration into a pdb formatted file.

    Also prints the cell parameters in standard pdb form. Note
    that the angles are in degrees.

    Args:
        atoms: An atoms object giving the atom positions.
        cell: A cell object giving the system box.
        filedesc: An open writable file object. Defaults to standard output.
        title: An optional string of max. 70 characters.
    """

    fmt_cryst = "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%s%4i\n"
    fmt_atom = (
        "ATOM  %5i %4s%1s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2i\n"
    )

    if title != "":
        filedesc.write("TITLE   %70s\n" % (title))

    a, b, c, alpha, beta, gamma = mt.h2abc_deg(cell.h * cell_conv)

    z = 1
    filedesc.write(fmt_cryst % (a, b, c, alpha, beta, gamma, " P 1        ", z))

    natoms = atoms.natoms
    qs = dstrip(atoms.q) * atoms_conv
    lab = dstrip(atoms.names)
    for i in range(natoms):
        data = (
            i + 1,
            lab[i],
            " ",
            "  1",
            " ",
            1,
            " ",
            qs[3 * i],
            qs[3 * i + 1],
            qs[3 * i + 2],
            0.0,
            0.0,
            "  ",
            0,
        )
        filedesc.write(fmt_atom % data)

    filedesc.write("END\n")


def read_pdb(filedesc):
    """Reads a PDB-style file and creates an Atoms and Cell object.

    Args:
        filedesc: An open readable file object from a pdb formatted file.

    Returns:
        An Atoms object with the appropriate atom labels, masses and positions,
        and a Cell object with the appropriate cell dimensions and an estimate
        of a reasonable cell mass.
    """

    header = filedesc.readline()
    comment = "# comment line might contain special i-PI keywords "
    if "TITLE" in header:
        # skip the comment field
        comment = copy.copy(header)
        header = filedesc.readline()
    # PDB defaults to Angstrom, because so says the standard
    if "positions{" not in comment:
        comment = comment.strip()
        comment += " positions{angstrom}\n"
    if "cell{" not in comment:
        comment = comment.strip()
        comment += " cell{angstrom}\n"
    if header == "":
        raise EOFError("End of file or empty header in PDB file")

    a = float(header[6:15])
    b = float(header[15:24])
    c = float(header[24:33])
    alpha = float(header[33:40])
    beta = float(header[40:47])
    gamma = float(header[47:54])
    alpha *= np.pi / 180.0
    beta *= np.pi / 180.0
    gamma *= np.pi / 180.0
    h = mt.abc2h(a, b, c, alpha, beta, gamma)
    cell = h

    natoms = 0
    body = filedesc.readline()
    qatoms = []
    names = []
    masses = []
    while body.strip() != "" and body.strip() != "END":
        natoms += 1
        name = body[12:16].strip()
        names.append(name)
        masses.append(Elements.mass(name))
        x = float(body[31:39])
        y = float(body[39:47])
        z = float(body[47:55])
        qatoms.append(x)
        qatoms.append(y)
        qatoms.append(z)

        body = filedesc.readline()

    return (
        comment,
        cell,
        np.asarray(qatoms),
        np.asarray(names, dtype="|U4"),
        np.asarray(masses),
    )

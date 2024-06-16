"""Functions used to read input configurations and print trajectories
using the Atomic Simulation Environment
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

# import os
import sys
import re
import numpy as np
from ipi.utils.units import Constants


try:
    import ase
except ImportError:
    ase = None

__all__ = ["print_ase", "read_ase"]


def _asecheck():
    if ase is None:
        raise RuntimeError(
            """The Atomic Simulation Environment must be installed to use the mode 'ase'. Please run

$~ pip install -U ase

to install it."""
        )

    from packaging import version

    assert version.parse(ase.__version__) >= version.parse(
        "3.22.1"
    ), "Please get a newer version of ASE"


def print_ase(
    atoms, cell, filedesc=sys.stdout, title="", cell_conv=1.0, atoms_conv=1.0
):
    """Prints an atomic configuration into an ASE extended XYZ formatted file.

    Args:
        atoms: An atoms object giving the centroid positions.
        cell: A cell object giving the system box.
        filedesc: An open writable file object. Defaults to standard output.
        title: This gives a string to be appended to the comment line.
    """

    _asecheck()

    from ase import Atoms

    # Reads the quantity from the ipi_comment
    quantity = next(
        (
            (keyword, value)
            for keyword, value in re.findall(r"(\w+)\{(\w+)\}", title)
            if keyword != "cell"
        ),
        (None, None),
    )[0]

    ase_atoms = Atoms(
        atoms.names,
        cell=cell.h.T * cell_conv,
        pbc=True,
    )
    ase_atoms.arrays[quantity] = atoms.q.reshape((-1, 3)) * atoms_conv
    ase_atoms.info["ipi_comment"] = title
    ase_atoms.write(filedesc, format="extxyz")


# this is a ugly hack to support reading iteratively from a file descriptor while coping with the quirks of both ASE and i-PI
_ase_open_files = {}


def read_ase(filedesc):
    """Reads an ASE-compatible file and returns data in raw format for further units transformation
    and other post processing.

    Args:
        filedesc: An open readable file object from an ase-readable file.

    Returns:
        i-Pi comment line, cell array, data (positions, forces, etc.), atoms names and masses
    """

    _asecheck()

    from ase.io import iread

    # keep a list of open files to iterate on using iread
    if filedesc not in _ase_open_files:
        _ase_open_files[filedesc] = iread(filedesc, ":")

    try:
        try:
            atoms = next(_ase_open_files[filedesc])
        except ase.io.formats.UnknownFileTypeError:  # defaults to extxyz
            _ase_open_files[filedesc] = iread(filedesc, ":", format="extxyz")
            atoms = next(_ase_open_files[filedesc])
    except StopIteration:
        _ase_open_files.pop(filedesc)
        raise EOFError()

    if all(atoms.get_pbc()):
        # We want to make the cell conform
        a = atoms.cell[0]
        atoms.rotate(a, "x", rotate_cell=True)  # a along x
        b = atoms.cell[1]
        b = b.copy() / np.linalg.norm(b)
        ang = -np.arctan2(b[2], b[1]) * 180 / np.pi
        atoms.rotate(ang, "x", rotate_cell=True)  # b in xy

    comment = "Structure read with ASE with composition %s" % atoms.symbols.formula
    cell = atoms.cell.array.T
    qatoms = atoms.positions.reshape((-1))
    names = list(atoms.symbols)
    masses = atoms.get_masses() * Constants.amu

    return comment, cell, qatoms, names, masses

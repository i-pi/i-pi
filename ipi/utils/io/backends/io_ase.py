"""Functions used to read input configurations and print trajectories
using the Atomic Simulation Environment
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

import os
import sys

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


def print_ase(
    atoms, cell, filedesc=sys.stdout, title="", cell_conv=1.0, atoms_conv=1.0
):
    """Prints an atomic configuration into an XYZ formatted file.

    Args:
        atoms: An atoms object giving the centroid positions.
        cell: A cell object giving the system box.
        filedesc: An open writable file object. Defaults to standard output.
        title: This gives a string to be appended to the comment line.
    """

    _asecheck()

    from ase import Atoms

    atoms = Atoms(atoms.names, positions=atoms.q * atoms_conv, cell=cell.h * cell_conv)
    atoms.write(filedesc)


def read_ase(filedesc):
    """Reads an ASE-compatible file and returns data in raw format for further units transformation
    and other post processing.

    Args:
        filedesc: An open readable file object from an ase-readable file.

    Returns:
        i-Pi comment line, cell array, data (positions, forces, etc.), atoms names and masses
    """

    _asecheck()

    from ase.io import read
    from ase.build import niggli_reduce

    # Extract the format
    fmt = os.path.splitext(filedesc.name)[1].strip('.')

    try:
        a = read(filedesc, format=fmt)
    except ValueError:
        raise EOFError()

    if all(a.get_pbc()):
        niggli_reduce(a)

    comment = "Structure read with ASE with composition %s" % a.symbols.formula
    cell = a.cell.array
    qatoms = a.positions.reshape((-1))
    names = list(a.symbols)
    masses = a.get_masses()

    return comment, cell, qatoms, names, masses

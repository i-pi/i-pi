"""Functions used to read input configurations and print trajectories
in the JSON format.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2015 i-PI developers
# See the "licenses" directory for full license information.


import sys
import json

import numpy as np

import ipi.utils.mathtools as mt
from ipi.utils.depend import dstrip


__all__ = ["print_json_path", "print_json", "read_json", "iter_json"]


def print_json_path(beads, cell, filedesc=sys.stdout):
    """Prints all the bead configurations into a JSON formatted file.

    Prints all the replicas for each time step separately, rather than all at
    once.

    Args:
        beads: A beads object giving the bead positions.
        cell: A cell object giving the system box.
        filedesc: An open writable file object. Defaults to standard output.
    """
    raise NotImplementedError("print_json_path is not implemented yet.")


def print_json(
    atoms, cell, filedesc=sys.stdout, title="", cell_conv=1.0, atoms_conv=1.0
):
    """Prints an atomic configuration into an XYZ formatted file.

    Args:
        atoms: An atoms object giving the centroid positions.
        cell: A cell object giving the system box.
        filedesc: An open writable file object. Defaults to standard output.
        title: This gives a string to be appended to the comment line.
    """

    a, b, c, alpha, beta, gamma = mt.h2abc_deg(cell.h * cell_conv)

    natoms = atoms.natoms
    # direct access to avoid unnecessary slow-down
    qs = dstrip(atoms.q) * atoms_conv
    lab = dstrip(atoms.names)

    data = {}
    data["natoms"] = natoms
    data["cell"] = [a, b, c, alpha, beta, gamma]
    data["title"] = title
    data["q"] = qs.tolist()
    data["labels"] = lab.tolist()

    filedesc.write(json.dumps(data))
    filedesc.write(" \n")


def read_json(filedesc):
    """Reads a JSON-style file with i-pi style comments and creates an Atoms and Cell object.

    Args:
        filedesc: An open readable file object from a json formatted file with i-PI header comments.

    Returns:
        An Atoms object with the appropriate atom labels, masses and positions.
        A Cell object.
    """

    try:
        line = filedesc.readline()
        data = json.loads(line)
    except ValueError:
        raise EOFError("The file descriptor hit EOF.")

    line.strip("\n")
    qatoms = np.asarray(data["q"], dtype=float)
    names = np.asarray(data["labels"], dtype="|S4")
    title = data["title"]
    masses = np.zeros(len(names))

    a, b, c, alpha, beta, gamma = data["cell"]
    alpha *= np.pi / 180.0
    beta *= np.pi / 180.0
    gamma *= np.pi / 180.0
    cell = mt.abc2h(a, b, c, alpha, beta, gamma)

    return (title, cell, qatoms, names, masses)


def iter_json(filedesc, **kwargs):
    """Takes a json-style file and yields one Atoms object after another.

    Args:
        filedesc: An open readable file object from a json formatted file.

    Returns:
        Generator over the json trajectory, that yields
        Atoms objects with the appropriate atom labels, masses and positions.
    """

    try:
        while True:
            yield read_json(filedesc, **kwargs)
    except EOFError:
        pass

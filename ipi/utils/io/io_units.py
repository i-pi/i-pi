"""Functions used to transform units in input files into default atomic system of units.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import re
from ipi.utils.messages import verbosity, info
from ipi.engine.atoms import Atoms
from ipi.engine.cell import Cell
from ipi.utils.units import unit_to_internal
from ipi.engine.properties import Trajectories as Traj

__all__ = ["auto_units", "process_units"]
# Regular expressions initialization for read_xyz function
cell_unit_re = re.compile(r"cell\{([A-Za-z_]*)\}")  # cell unit pattern
traj_dict = Traj().traj_dict  # trajectory dictionary
traj_re = [
    re.compile("%s%s" % (key, r"\{[A-Za-z_]*\}")) for key in list(traj_dict.keys())
]  # trajectory patterns


def auto_units(
    comment="",
    dimension="automatic",
    units="automatic",
    cell_units="automatic",
    mode="xyz",
):
    """Processes comment line and requested units to determine how to interpret the I/O conversion."""
    # heuristics to detect units
    if mode == "pdb":  # these are the default units
        auto_cell = "angstrom"
        auto_units = "angstrom"
        auto_dimension = "length"
    elif mode == "ase":
        auto_cell = "ase"
        auto_units = "ase"
        auto_dimension = "undefined"
    else:
        auto_cell = "atomic_unit"
        auto_units = "atomic_unit"
        auto_dimension = "undefined"

    is_comment_useful = []
    if comment != "":  # but they can be overridden by a special comment line
        # tries to guess units from the input
        # Extracting trajectory units
        is_comment_useful = [
            _f for _f in [key.search(comment.strip()) for key in traj_re] if _f
        ]
        if len(is_comment_useful) > 0:
            traj = is_comment_useful[0].group()[:-1].split("{")
            auto_dimension, auto_units = traj_dict[traj[0]]["dimension"], traj[1]

        # Extracting cell units
        tmp = cell_unit_re.search(comment)
        if tmp is not None:
            auto_cell = tmp.group(1)
    if dimension == "automatic":
        dimension = auto_dimension
    elif dimension != auto_dimension and len(is_comment_useful) > 0:
        raise ValueError(
            "Requested dimension mismatch with property indicated in the comment string"
        )

    if units == "automatic":
        units = auto_units
    elif units != auto_units and len(is_comment_useful) > 0:
        raise ValueError(
            "Requested units mismatch with units indicated in the comment string"
        )

    if cell_units == "automatic":
        cell_units = auto_cell
    elif cell_units != auto_cell and len(is_comment_useful) > 0:
        raise ValueError(
            "Requested cell units mismatch with units indicated in the comment string"
        )

    return dimension, units, cell_units


def process_units(
    comment,
    cell,
    data,
    names,
    masses,
    natoms,
    dimension="automatic",
    units="automatic",
    cell_units="automatic",
    mode="xyz",
):
    """Convert the data in the file according to the units written in the i-PI format.

    Args:
        comment:
        cell:
        data:
        names:
        masses:
        output:

    Returns:

    """
    dimension, units, cell_units = auto_units(
        comment, dimension, units, cell_units, mode
    )

    info(
        " # Interpreting input with dimension %s, units %s and cell units %s"
        % (dimension, units, cell_units),
        verbosity.high,
    )

    # Units transformation
    cell *= unit_to_internal("length", cell_units, 1)  # cell units transformation
    data *= unit_to_internal(dimension, units, 1)  # units transformation

    # Return data as i-PI structures
    cell = Cell(cell)
    atoms = Atoms(natoms)
    atoms.q[:] = data
    atoms.names[:] = names
    atoms.m[:] = masses

    return {"atoms": atoms, "cell": cell, "comment": comment}

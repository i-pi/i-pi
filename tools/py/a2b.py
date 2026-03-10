#!/usr/bin/env python3

"""a2b.py

Reads positions of a system in format 'a' and returns a
file written in format 'b'.

Assumes the input file is in a format 'a'.

Syntax:
   a2b.py filename format_a format_b

"""


import sys
import re
from ipi.utils.io import read_file, print_file, read_file_raw
from ipi.engine.properties import Trajectories as Traj
from ipi.utils.messages import verbosity

cell_unit_re = re.compile(r"cell\{([A-Za-z_]*)\}")  # cell unit pattern
traj_dict = Traj().traj_dict  # trajectory dictionary
traj_re = [
    re.compile("%s%s" % (key, r"\{[A-Za-z_]*\}")) for key in list(traj_dict.keys())
]  # trajectory patterns
verbosity.level = "low"


def get_cell_units(comment, mode):
    # default units
    if mode == "pdb":
        auto_cell = "angstrom"
    else:
        auto_cell = "atomic_unit"

    cunit = cell_unit_re.search(comment)
    if cunit is not None:
        auto_cell = cunit.group(1)
    return auto_cell


def get_key_dim_units(comment, mode):
    # if mode is pdb return standard units and dimensions.
    if mode == "pdb":
        auto_units = "angstrom"
        auto_dimension = "length"

    # if mode is not pdb looks for the units in the comment line.
    else:
        auto_units = "atomic_unit"
        auto_dimension = "undefined"

    is_comment_useful = []
    if comment != "":
        is_comment_useful = [
            _f for _f in [key.search(comment.strip()) for key in traj_re] if _f
        ]
        if len(is_comment_useful) > 0:
            traj = is_comment_useful[0].group()[:-1].split("{")
            key, auto_dimension, auto_units = (
                traj[0],
                traj_dict[traj[0]]["dimension"],
                traj[1],
            )

    if mode == "pdb" and auto_dimension != "length":
        raise ValueError(
            "PDB Standard is only designed for atomic positions with units in Angstroms"
        )

    return key, auto_dimension, auto_units


def main(filename, imode, omode):
    ipos = open(filename, "r")
    ifr = 0

    # extracts the dimension, its units and the cell_units from the first frame
    ret = read_file_raw(imode, ipos)
    ipos.close()
    comment = ret["comment"]
    cell_units = get_cell_units(comment, imode)
    key, dim, dim_units = get_key_dim_units(comment, imode)

    ipos = open(filename, "r")
    while True:
        try:
            ret = read_file(imode, ipos)
            pos = ret["atoms"]
            cell = ret["cell"]

        except EOFError:  # finished reading files
            sys.exit(0)
        print_file(
            omode,
            pos,
            cell,
            filedesc=sys.stdout,
            title="",
            key=key,
            dimension=dim,
            units=dim_units,
            cell_units=cell_units,
        )
        ifr += 1


if __name__ == "__main__":
    main(*sys.argv[1:])

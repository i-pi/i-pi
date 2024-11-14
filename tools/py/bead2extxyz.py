"""
Convert i-PI bead trajectories to ASE-readable extxyz files.

This script reads bead trajectories from an i-PI XML file and creates separate extxyz files for each bead. 
The output files can be used with ASE for further analysis.

Usage
-----
Basic conversion:
    $ python bead2extxyz.py -i input.xml -o output_prefix

Include additional arrays (e.g., velocities):
    $ python bead2extxyz.py -i input.xml -o output_prefix -t velocities

Add specific properties (e.g., bead-level potential):
    $ python bead2extxyz.py -i input.xml -o output_prefix -p bead_potentials

Arguments
---------
-i, --input : str
    Path to the XML input file (required).
-o, --output : str
    Prefix for output extxyz files (required).
-t, --trajectories : list of str (optional)
    Arrays to include in extxyz files. Default: `["forces"]`.
-p, --properties : list of str (optional)
    Properties to include as `info` fields. Default: `["bead_potentials"]`.

Requirements
------------
- Python with `argparse`.
- `ase` package (install via `$ pip install ase`).

Examples
--------
Basic usage:
    $ python bead2extxyz.py -i trajectory.xml -o output_prefix

Include velocities:
    $ python bead2extxyz.py -i trajectory.xml -o output_prefix -t velocities
"""

import argparse
from ipi.utils.parsing import create_bead_trajectories

try:
    import ase
    from ase.io import write
except ImportError:
    ase = None

description = "Convert a classical i-PI trajectory to an extxyz file readable by ASE."


# ---------------------------------------#
def prepare_args():
    parser = argparse.ArgumentParser(description=description)
    argv = {"metavar": "\b"}
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        **argv,
        help="xml input file",
    )
    parser.add_argument(
        "-t",
        "--trajectories",
        type=str,
        nargs="*",  # Allows zero or more arguments to be passed, creating a list
        required=False,
        **argv,
        help="trajectories to be added as `arrays` to the extxyz file (default: %(default)s)",
        default=["forces"],
    )
    parser.add_argument(
        "-p",
        "--properties",
        type=str,
        nargs="*",  # Allows zero or more arguments to be passed, creating a list
        required=False,
        **argv,
        help="properties to be added as `info` to the extxyz file (default: %(default)s)",
        default=["bead_potentials"],
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        **argv,
        help="prefix of the extxyz output file",
    )
    return parser.parse_args()


# ---------------------------------------#
def main():

    if ase is None:
        raise ImportError(
            "'bead2extxyz.py' requires the `ase` package to return the structure in ASE format"
        )

    args = prepare_args()
    print(f"\n\t {description}")
    print("\n\t Input arguments:")
    for k in args.__dict__.keys():
        print("\t {:>20s}:".format(k), getattr(args, k))
    print()

    print("\t Constructing the classical trajectory.")
    list_atoms = create_bead_trajectories(
        args.input, args.trajectories, args.properties
    )
    print(
        f"\n\t The trajectories have {len(list_atoms[0])} snapshots.   | The least common multiple of all the strides has been used."
    )
    print(
        "\t The trajectories have the following arrays: ",
        list(list_atoms[0][0].arrays.keys()),
    )
    print(
        "\t The trajectories have the following info: ",
        list(list_atoms[0][0].info.keys()),
    )

    print(f"\n\t The trajectories will be saved to file with prefix '{args.output}':")
    for n, atoms in enumerate(list_atoms):
        file = f"{args.output}.bead={n}.extxyz"
        print(f"\t - bead {n} to file {file}")
        write(file, atoms, format="extxyz")

    print("\t Job done :)\n")

    return


# ---------------------------------------#
if __name__ == "__main__":
    main()

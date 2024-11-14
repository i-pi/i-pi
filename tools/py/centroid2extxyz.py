#!/usr/bin/env python
"""
This script converts the i-PI centroid trajectory into an extxyz file format readable by ASE.
It can be used to include various trajectory arrays and saves the results in a format that is compatible with ASE for further processing or analysis.

Usage
-----
To use this script, run the following command in the terminal:

1. Basic usage to convert an XML trajectory to an extxyz file with default options:
    $ python centroid2extxyz.py -i input.xml -o output.extxyz

2. Including additional trajectory arrays (e.g., forces and velocities):
    $ python centroid2extxyz.py -i input.xml -o output.extxyz -t forces velocities

Arguments
---------
-i, --input : str
    The XML input file containing the centroid trajectory. This argument is required.

-o, --output : str
    The name of the output extxyz file. This argument is required.

-t, --trajectories : list of str (optional)
    A list of trajectory data to be added as arrays to the extxyz file. Default is `["forces"]` (and `positions` will be added automatically).

Requirements
------------
- Python with `argparse` for argument parsing.
- The `ase` package for handling ASE trajectory objects. Install it via:
    $ pip install ase

Exceptions
----------
- Raises `ImportError` if the `ase` package is not installed.
- Prints relevant error messages for missing input or output arguments.

"""

# author: Elia Stocco

import argparse
from ipi.utils.parsing import create_centroid_trajectory

try:
    import ase
    from ase.io import write
except ImportError:
    ase = None

description = (
    "Convert the i-PI trajectory for the centroid to an extxyz file readable by ASE."
)


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
        "-o", "--output", type=str, required=True, **argv, help="output extxyz file"
    )
    return parser.parse_args()


# ---------------------------------------#
def main():

    if ase is None:
        raise ImportError(
            "'centroid2extxyz.py' requires the `ase` package to return the structure in ASE format"
        )

    args = prepare_args()
    print(f"\n\t {description}")
    print("\n\t Input arguments:")
    for k in args.__dict__.keys():
        print("\t {:>20s}:".format(k), getattr(args, k))
    print()

    print("\t Constructing the centroid trajectory.")
    atoms = create_centroid_trajectory(args.input, args.trajectories)
    print(
        f"\t The trajectory has {len(atoms)} snapshots.   | The least common multiple of all the strides has been used."
    )
    print(
        "\t The trajectory has the following arrays: ",
        list(atoms[0].arrays.keys()),
    )
    print(f"\t The trajectory will be saved to '{args.output}'.")
    write(args.output, atoms, format="extxyz")

    print("\t Job done :)\n")

    return


# ---------------------------------------#
if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""
This script converts a classical i-PI trajectory into an extxyz file format readable by ASE.
The script must be run in the same folder where the i-PI simulation was done.

Usage
-----
To use this script, run the following command in the terminal:

1. Basic usage to convert an XML trajectory to an extxyz file with default options:
    $ python classical_md2extxyz.py -i input.xml -o output.extxyz

2. Including additional trajectory arrays (e.g., forces and velocities):
    $ python classical_md2extxyz.py -i input.xml -o output.extxyz -t forces velocities

3. Including additional properties (e.g., potential energy and temperature):
    $ python classical_md2extxyz.py -i input.xml -o output.extxyz -p potential temperature

By default `positions` and `forces` are always read, as well as `potential`.
To change this behavor, just modify the default parameters in `prepare_args`.

Arguments
---------
-i, --input : str
    The XML input file containing the classical i-PI trajectory. This argument is required.

-o, --output : str
    The name of the output extxyz file. This argument is required.

-t, --trajectories : list of str (optional)
    A list of trajectory data to be added as arrays to the extxyz file. Default is `["forces"]` (and `positions` will be added automatically).

-p, --properties : list of str (optional)
    A list of properties to be added as `info` fields in the extxyz file. Default is `["potential"]`.

Requirements
------------
- Python with `argparse` for argument parsing.
- The `ase` package for working with ASE trajectory objects. Install it via:
    $ pip install ase

Exceptions
----------
- Raises `ImportError` if the `ase` package is not installed.
- Raises `ValueError` if the input or output arguments are missing.

"""

import argparse
from ipi.utils.parsing import create_classical_trajectory

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
        default=["potential"],
    )
    parser.add_argument(
        "-o", "--output", type=str, required=True, **argv, help="output extxyz file"
    )
    return parser.parse_args()


# ---------------------------------------#
def main():

    if ase is None:
        raise ImportError(
            "classical_md2extxyz.py requires the `ase` package to return the structure in ASE format"
        )

    args = prepare_args()
    print(f"\n\t {description}")
    print("\n\t Input arguments:")
    for k in args.__dict__.keys():
        print("\t {:>20s}:".format(k), getattr(args, k))
    print()

    print("\t Constructing the classical trajectory.")
    atoms = create_classical_trajectory(args.input, args.trajectories, args.properties)
    print(
        f"\t The trajectory has {len(atoms)} snapshots.   | The least common multiple of all the strides has been used."
    )
    print(
        "\t The trajectory has the following arrays: ",
        list(atoms[0].arrays.keys()),
    )
    print(
        "\t The trajectory has the following info: ",
        list(atoms[0].info.keys()),
    )

    print(f"\t The trajectory will be saved to '{args.output}'.")
    write(args.output, atoms, format="extxyz")

    print("\t Job done :)\n")

    return


# ---------------------------------------#
if __name__ == "__main__":
    main()

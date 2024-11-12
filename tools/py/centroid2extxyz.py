#!/usr/bin/env python
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
            "centroid2extxyz.py requires the `ase` package to return the structure in ASE format"
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

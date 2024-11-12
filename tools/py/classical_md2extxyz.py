#!/usr/bin/env python
import argparse
from ipi.utils.parsing import create_classical_trajectory

try:
    import ase
    from ase.io import write
except ImportError:
    ase = None

description = "Convert a classical i-PI trajectory to an extxyz file readable from ASE."


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
        help="trajectories to be added to the extxyz file (default: %(default)s)",
        default=["forces"],
    )
    # parser.add_argument(
    #     "-f",
    #     "--folder",
    #     type=str,
    #     required=False,
    #     **argv,
    #     help="folder containing the i-PI files (default: %(default)s)",
    #     default=".",
    # )
    # parser.add_argument(
    #     "-n",
    #     "--nbeads",
    #     type=int,
    #     required=True,
    #     **argv,
    #     help="number of beads (default: %(default)s)",
    #     default=":",
    # )
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

    print("\t Constructing the classical trajectory")
    atoms = create_classical_trajectory(args.input, args.trajectories)
    print(
        "\t The trajectory will have the following arrays: ",
        list(atoms[0].arrays.keys()),
    )

    print(f"\t The trajectory will be saved to '{args.output}'.")
    write(args.output, atoms)

    print("\t Job dobe :)")

    return


# ---------------------------------------#
if __name__ == "__main__":
    main()

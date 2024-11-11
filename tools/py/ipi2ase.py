#!/usr/bin/env python
import argparse
from ipi.utils.parsing import merge_beads

description = "Convert i-PI trajectories to ASE files."


# ---------------------------------------#
def prepare_args():
    parser = argparse.ArgumentParser(description=description)
    argv = {"metavar": "\b"}
    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        required=True,
        **argv,
        help="prefix of the i-PI files",
    )
    parser.add_argument(
        "-f",
        "--folder",
        type=str,
        required=False,
        **argv,
        help="folder containing the i-PI files (default: %(default)s)",
        default=".",
    )
    parser.add_argument(
        "-n",
        "--nbeads",
        type=int,
        required=True,
        **argv,
        help="number of beads (default: %(default)s)",
        default=":",
    )
    parser.add_argument(
        "-o", "--output", type=str, required=True, **argv, help="output extxyz file"
    )
    return parser.parse_args()


# ---------------------------------------#
def main():

    # ------------------#
    args = prepare_args()
    print(f"\n\t {description}")
    print("\n\t Input arguments:")
    for k in args.__dict__.keys():
        print("\t {:>20s}:".format(k), getattr(args, k))
    print()

    # ------------------#
    merge_beads(args.prefix, args.folder, args.nbeads, args.output)


# ---------------------------------------#
if __name__ == "__main__":
    main()

#!/usr/bin/env python3


import os
import sys
import numpy as np
import argparse
from ipi.utils.messages import verbosity
from ipi.utils.io.io_units import auto_units, process_units

from ipi.utils.io import open_backup, iter_file_name_raw, print_file


description = """
Read positions of individual beads from a set of trajectory files and
a) multiplex them into a single output trajectory
b) wrap/unwrap them.
Trajectory file formats are inferred from file extensions, the number
of beads is given by the number of input files.
"""


def main(fns_in, fn_out, begin, end, stride, wrap, unwrap):
    verbosity.level = "low"
    print("Multiplexing {:d} beads into one trajectory.".format(len(fns_in)))
    print()

    print("input file names:")
    for fn in fns_in:
        print(fn)
    print()

    print("output file name:", fn_out)
    print()

    # Open input trajectory iterators.
    trjs_in = [iter_file_name_raw(fn) for fn in fns_in]
    mode = os.path.splitext(fns_in[0])[-1]

    # Open output file.
    f_out = open_backup(fn_out, "w")
    mode_out = os.path.splitext(fn_out)[-1]

    # Loop over all frames.
    i_frame = 0
    i_frame_saved = 0

    # There can be multiple trajectories, so we store a frame_last for each trajectory
    frame_last = [None] * len(fns_in)

    while True:
        # Check the endpoint index, exit if we're done.
        if (end > -1) and (i_frame >= end):
            break

        # Should we save output from this frame?
        do_output = (i_frame >= begin) and ((i_frame % stride) == 0)

        try:
            # Get the frames from all trajectories...
            for idx, trj in enumerate(trjs_in):
                frame = next(trj)

                # gets units from first frame
                dimension, units, cell_units = auto_units(comment=frame["comment"])
                frame = process_units(
                    dimension=dimension,
                    units=units,
                    cell_units=cell_units,
                    mode=mode,
                    **frame
                )

                if wrap:
                    frame = wrap_positions(frame)
                if unwrap:
                    frame = unwrap_positions(frame, frame_last[idx])

                frame_last[idx] = frame.copy()

                # ... and possibly save them in the output trajectory.
                if do_output:
                    print_file(
                        mode_out,
                        frame["atoms"],
                        frame["cell"],
                        f_out,
                        dimension=dimension,
                        units=units,
                        cell_units=cell_units,
                    )
            if do_output:
                i_frame_saved += 1
        except StopIteration:
            # Stop when any of the trajectories runs out of frames.
            break

        # Count frames and print information on progress.
        i_frame += 1
        if i_frame % 100 == 0:
            print("\rframe {:d}".format(i_frame), end="")
        sys.stdout.flush()

    f_out.close()

    print()
    print()
    print("Loaded {:d} frames.".format(i_frame))
    print("Saved {:d} frames.".format(i_frame_saved))


def wrap_positions(frame):
    pos = frame["atoms"].q.copy()
    pos.shape = (len(pos) / 3, 3)
    cell = frame["cell"].h
    cell_inv = np.linalg.inv(cell)
    s = np.dot(cell_inv, pos.T)
    s -= np.round(s)
    frame["atoms"].q = np.dot(cell, s).T.flatten()
    return frame


def unwrap_positions(frame, framelast):
    if framelast is None:
        return frame

    poslast_uwr = framelast["atoms"].q.copy()
    poslast_uwr.shape = (len(poslast_uwr) / 3, 3)

    pos = frame["atoms"].q.copy()
    pos.shape = (len(pos) / 3, 3)

    d = pos - poslast_uwr

    cell = frame["cell"].h
    cell_inv = np.linalg.inv(cell)
    s = np.dot(cell_inv, d.T)
    s -= np.round(s)
    d = np.dot(cell, s).T
    frame["atoms"].q = (poslast_uwr + d).flatten()
    return frame


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "filenames", type=str, nargs="+", help="Bead trajectories to process."
    )
    parser.add_argument("--filename-out", type=str, help="Output file name.")
    parser.add_argument("--begin", type=int, default=0, help="Step to begin.")
    parser.add_argument("--end", type=int, default=-1, help="Step to end.")
    parser.add_argument("--stride", type=int, default=1, help="Stride in steps.")
    parser.add_argument(
        "--wrap", action="store_true", default=False, help="Wrap atomic positions."
    )
    parser.add_argument(
        "--unwrap", action="store_true", default=False, help="Unwrap atomic positions."
    )

    args = parser.parse_args()

    # Process everything.
    main(
        args.filenames,
        args.filename_out,
        args.begin,
        args.end,
        args.stride,
        args.wrap,
        args.unwrap,
    )

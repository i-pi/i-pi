#!/usr/bin/env python2

from __future__ import print_function

import os
import sys
import argparse
from ipi.utils.messages import verbosity

from ipi.utils.io import open_backup, iter_file_name, print_file


description = """
Read positions of individual beads from a set of trajectory files and
multiplex them into a single output trajectory. Trajectory file formats are
inferred from file extensions, the number of beads is given by the number of
input files.
"""


def main(fns_in, fn_out, begin, end, stride):

    verbosity.level = "low"
    print('Multiplexing {:d} beads into one trajectory.'.format(len(fns_in)))
    print()

    print('input file names:')
    for fn in fns_in:
        print(fn)
    print()

    print('output file name:', fn_out)
    print()

    # Open input trajectory iterators.
    trjs_in = [iter_file_name(fn) for fn in fns_in]

    # Open output file.
    f_out = open_backup(fn_out, 'w')
    mode_out = os.path.splitext(fn_out)[1]

    # Loop over all frames.
    i_frame = 0
    i_frame_saved = 0
    while True:

        # Check the endpoint index, exit if we're done.
        if (end > -1) and (i_frame >= end):
            break

        # Should we save output from this frame?
        do_output = (i_frame >= begin) and ((i_frame % stride) == 0)

        try:
            # Get the frames from all trajectories...
            for trj in trjs_in:
                frame = trj.next()
                # ... and possibly save them in the output trajectory.
                if do_output:
                    print_file(mode_out, frame['atoms'], frame['cell'], f_out)
            if do_output:
                i_frame_saved += 1
        except StopIteration:
            # Stop when any of the trajectories runs out of frames.
            break

        # Count frames and print information on progress.
        i_frame += 1
        if i_frame % 100 == 0:
            print('\rframe {:d}'.format(i_frame), end='')
        sys.stdout.flush()

    f_out.close()

    print()
    print()
    print('Loaded {:d} frames.'.format(i_frame))
    print('Saved {:d} frames.'.format(i_frame_saved))


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('filenames', type=str, nargs='+',
                        help='Bead trajectories to process.')
    parser.add_argument('--filename-out', type=str,
                        help='Output file name.')
    parser.add_argument('--begin', type=int, default=0,
                        help='Step to begin.')
    parser.add_argument('--end', type=int, default=-1,
                        help='Step to end.')
    parser.add_argument('--stride', type=int, default=1,
                        help='Stride in steps.')

    args = parser.parse_args()

    # Process everything.
    main(args.filenames, args.filename_out, args.begin, args.end, args.stride)

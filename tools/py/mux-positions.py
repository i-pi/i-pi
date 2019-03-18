#!/usr/bin/env python2

from __future__ import print_function

import os
import sys
import numpy as np
import argparse
from ipi.utils.messages import verbosity
import copy

from ipi.utils.io import open_backup, iter_file_name, print_file


description = """
Read positions of individual beads from a set of trajectory files and
a) multiplex them into a single output trajectory
b) wrap/unwrap them.
Trajectory file formats are inferred from file extensions, the number
of beads is given by the number of input files.
"""

def cell_dot_pos(cell,pos):
    s = np.array([])
    for row in pos:
        s = np.append(s,np.dot(cell,row))
    s.shape = (len(pos),3)
    #s = np.round(s,12)
    return s

def wrap_positions(pos, cell):
    if len(pos.shape) == 1:
        pos.shape = (len(pos)/3,3)
    cell_inv = np.linalg.inv(cell)
    s = cell_dot_pos(cell_inv,pos)
    sc = s - np.round(s)
    wrapped = cell_dot_pos(cell,sc)
    return wrapped

def unwrap_positions(pos_wr, poslast_wr, poslast_uwr, cell):
    if len(pos_wr.shape) == 1:
        pos_wr.shape = (len(pos_wr)/3,3)
    if len(poslast_wr.shape) == 1:
        poslast_wr.shape = (len(poslast_wr)/3,3)
    cell_inv = np.linalg.inv(cell)
    d = pos_wr - poslast_wr
    s = cell_dot_pos(cell_inv,d)
    sc = s-np.round(s)
    q = cell_dot_pos(cell,sc)
    return q+poslast_uwr

def unwrap_or_wrap(frame,pos_last_wr,pos_last_uwr,idx,args):
    frc = frame.copy()
    pos = frc['atoms'].q
    pos_= pos.copy()
    pos_wr = wrap_positions(pos_, frc['cell'].h)

    pos_uwr = False
    if args.unwrap and idx == 0:
        pos_uwr = pos_wr
    elif args.unwrap and idx > 0:
        pos_uwr = unwrap_positions(pos_wr, pos_last_wr, pos_last_uwr, frc['cell'].h)

    if args.wrap == True:
        frc['atoms'].q = pos_wr.flatten()
    elif args.unwrap == True:
        frc['atoms'].q = pos_uwr.flatten()
    frame = frc.copy()
    return frame, pos_wr, pos_uwr



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

                if args.wrap or args.unwrap:
                    if i_frame == 0:
                        pos_last_wr = False
                        pos_last_uwr = False
                    frame,pos_last_wr,pos_last_uwr = unwrap_or_wrap(frame,pos_last_wr,pos_last_uwr,i_frame,args)

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
    parser.add_argument('--wrap',  action='store_true', default=False,
                        help='Wrap atomic positions.')
    parser.add_argument('--unwrap',  action='store_true', default=False,
                        help='Unwrap atomic positions.')

    args = parser.parse_args()

    # Process everything.
    main(args.filenames, args.filename_out, args.begin, args.end, args.stride)

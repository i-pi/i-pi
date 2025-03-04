#!/usr/bin/env python3
"""Generate fake parameters useful to test the softare."""

import tempfile as tmp
import numpy as np
from ipi.utils.units import Elements
from copy import copy

all_elem = list(Elements.mass_list.keys())
at_names = None


def xyz_rand(natoms, comment, names=None):
    """Generate a fake xyz file. Atoms and coordinates are random."""

    global at_names
    xyz = np.random.random((natoms, 3))
    xyz = xyz * 10 - 5  # To have both, positive and negative numbers
    if not names:
        # if True:
        at_names = [all_elem[i] for i in np.random.randint(0, len(all_elem), natoms)]

    output = str(natoms) + "\n"

    if not comment.endswith("\n"):
        comment += "\n"
    output += comment
    for i in range(natoms):
        output += "%8s %12.5e %12.5e %12.5e\n" % (
            at_names[i],
            xyz[i, 0],
            xyz[i, 1],
            xyz[i, 2],
        )
    return (output, xyz.flatten(), copy(at_names))


def xyz_traj(natoms, nframe, comment):
    """Generate a fake xyz trajectory. Atoms and coordinates are random."""
    output, xyz, all_names = xyz_rand(natoms, comment)
    for _ in range(nframe - 1):
        raw_output = xyz_rand(natoms, comment, names=True)
        output += raw_output[0]
        xyz = np.concatenate([xyz, raw_output[1]])
        all_names += raw_output[2]

    return (output, xyz, all_names)


def xyz_traj_filedesc(natoms, nframe, comment):
    """Generate a file descriptor containing a fake xyz trajectory."""
    contents, xyz, all_names = xyz_traj(natoms, nframe, comment)
    filedesc = tmp.NamedTemporaryFile(mode="w")
    filedesc.write(contents)
    filedesc.seek(0)
    # This works only on unix!
    return (open(filedesc.name), xyz, all_names)


# if __name__ == '__main__':

# Fast autocheck... if the test is wrong itself... it is bad ;)
#     natoms = 100
#     raw = xyz_traj_filedesc(natoms, 200, '')

#     raw[0].seek(0)

#     frames = -1
#     atoms = 0
#     for line in raw[0].readlines():
#         fields = line.split()
#         if len(fields) == 1 and fields[0] != '':
#             print 'Check number of atoms', fields[0], '== ' + str(natoms)
#             assert fields[0] == str(natoms)
#             frames += 1
#             atoms = 0
#         if len(fields) > 1:
# print 'Atom #', atoms, 'Frame #', frames+1
# print 'Check name of atoms #' + str(fields[0]) +\
# '# == #' + str(raw[2][atoms]) + '#'
#             assert fields[0] == raw[2][atoms]
#             for _i in xrange(3):
# print 'Check coordinate of atoms #'+str(frames* natoms *3 + atoms*3+_i) + \
#                     ' ==> ' + str(float(fields[_i+1]) -
#                                   float(raw[1][frames * atoms*3+_i]))+\
#                     ' < 1E-9 ('+str(fields[_i+1])  + '  ' +\
#                     str(raw[1][frames * atoms*3+_i]) +')'

#                 assert abs(float(fields[_i+1]) - \
#                            float(raw[1][frames* natoms *3 + atoms*3+_i])) < 1E-9
#             atoms += 1

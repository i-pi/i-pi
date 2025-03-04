#!/usr/bin/env python3
"""Generate fake parameters useful to test the softare."""

import tempfile as tmp
import numpy as np
from ipi.utils.units import Elements
from copy import copy

all_elem = list(Elements.mass_list.keys())
at_names = None


def xyz_rand(natoms, comment, names=None):
    """Generate a fake pdb pseudo-file. Atoms and coordinates are random."""

    global at_names
    xyz = np.random.random((natoms, 3))
    xyz = xyz * 10 - 5  # To have both, positive and negative numbers
    if not names:
        # if True:
        at_names = [all_elem[i] for i in np.random.randint(0, len(all_elem), natoms)]

    if not comment.endswith("\n"):
        comment += "\n"
    output = comment
    for i in range(natoms):
        output += (
            "%-6s"  # Record name                       1- 6  1
            "%5i"  # Atom serial number                7-11  2
            " "  # Space                               12
            "%3s"  # Atom name                        13-16  3
            " "  # Alternate location indicator (?)    17
            "%3s"  # Residue name                     18-20  4
            " "  # Space                               21
            "%1s"  # Chain identifier                    22  5
            "%3i"  # Residue sequence number          23-26  6
            "%1s"  # Code for insertion of residues(?)   27  7
            "   "  # Spacese                          28-30
            "%8.3f"  # X coordinate                     31-38  8
            "%8.3f"  # Y coordinate                     39-46  9
            "%8.3f"  # Z coordinate                     47-54 10
            "%6.2f"  # Occupancy                        55-60 11
            "%6.2f"  # Temp Factor                      61-66 12
            "          "  # space                        67-76
            "%3s"  # Element symbol, right-justif     77-78 13
            "%3s\n"  # Charge on atom.                  79-80 14
            % (
                "ATOM",
                i,
                at_names[i],
                "UNK",
                " ",
                i,
                " ",
                xyz[i, 0],
                xyz[i, 1],
                xyz[i, 2],
                1.0,
                0.0,
                at_names[i].rjust(3),
                " ",
            )
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


if __name__ == "__main__":
    # Fast autocheck... if the test is wrong itself... it is bad ;)
    natoms = 100
    print(xyz_rand(natoms, "")[0])

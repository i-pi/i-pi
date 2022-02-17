#!/usr/bin/env python3
"""neb_interpolate.py
Creates interpolated path between 2 or more geometries. Beads are placed
equidistantly to a piecewise-linear path or a cubic spline formed by 2 or more
provided beads.

Arguments:
    --mode          - "endpoints" or "xyz" or "chk"
    -n, -N, --nbeads    - number of NEB replicas incl. endpoints
    EITHER
      --ini         - initial geometry (.xyz) if mode "endoints"
      --fin         - final geometry (.xyz) if mode "endpoints"
    OR
      --xyz         - .xyz with all geometries if mode "xyz"
    OR
      --chk         - i-PI checkpoint with >=2 beads.

Output:
    Writes a file "init.xyz" with N beads.
"""
# (c) Karen Fidanyan 2022

import os
import numpy as np
import sys
import argparse

from ipi.utils.io import read_file, read_file_raw, print_file
from ipi.utils.io.io_units import auto_units
from ipi.utils.units import unit_to_internal

try:
    import scipy
    from scipy.interpolate import make_interp_spline, splev
except Exception as e:
    scipy = None
    scipy_exception = e


def spline_resample(q, nbeads_old, nbeads_new, k=3):
    """Resamples the intermediate points along the spline so that
    all points are equidistant by the spline arc length.

    Arguments:
        q           - beads.q[:]
        k           - 1 or 3: order of the spline
        nbeads_old  - number of beads in the input path
        nbeads_new      - desired number of beads
    Returns:
        new_q - resampled coordinates
    """
    from numpy.linalg import norm as npnorm

    if scipy is None:
        print("spline interpolation requires scipy module.")
        raise (scipy_exception)

    if nbeads_new <= 2:
        softexit.trigger(status="bad", message="nbeads_new < 3 in string optimization.")

    # First, we calculate the current parameterization of the path
    # according to 3N-D Euclidean distances between adjacent beads.
    t = [0.0]
    current_t = 0.0
    print("Cartesian 3N-D distances between old beads:")
    for i in range(1, nbeads_old):
        dist = npnorm(q[i] - q[i - 1])
        print("\tfrom %3d to %3d : %.6f bohr" % (i - 1, i, dist))
        if dist <= 1e-2:
            print("Warning: two adjacent beads are very close.")
        current_t += dist
        t.append(current_t)

    t = np.array(t)
    t /= t[-1]
    print("Spline parameter 't' before resampling:")
    print(t)
    print("t[i] - t[i-1]:")
    print(t - np.roll(t, 1))

    # New parameterization, equidistant by spline parameter
    new_t = np.arange(nbeads_new) / (nbeads_new - 1)
    print("New t:")
    print(new_t)

    # This reshapes q to a list of 3*Natoms sets consisting of nbeads_old numbers each,
    # describing trajectory of one Cartesian component of an atom over the path.
    atoms = q.T

    new_atoms = []
    for component in atoms:
        # Interpolate the trajectory of one cartesian component by a spline.
        # bc_type='natural' imposes zero second derivative at the endpoints.
        if k == 3:
            spl = make_interp_spline(t, component, k=k, bc_type="natural")
        elif k == 1:
            spl = make_interp_spline(t, component, k=k)
        # Sample new parameter values uniformly from the old range.
        # Resample positions of this atom.
        new_atoms.append(splev(new_t, spl))

    # Reshape the path back to the beads shape
    new_q = np.array(new_atoms).T

    print("Cartesian 3N-D distances between NEW beads:")
    for i in range(1, nbeads):
        dist = npnorm(new_q[i] - new_q[i - 1])
        print("\tfrom %3d to %3d : %.6f bohr" % (i - 1, i, dist))

    return new_q


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""A script for linear interpolation of a NEB path"""
    )
    parser.add_argument(
        "-m",
        "--mode",
        type=str,
        default="endpoints",
        help="Mode: 'endpoints' to interpolate between 2 geometries provided separately,\n"
        "   'xyz' to interpolate between multiple geometries in a single file,\n"
        "   or 'chk' to load an existing checkpoint and resample path from it.",
    )
    parser.add_argument(
        "--ini",
        type=str,
        default="ini.xyz",
        help="Filename of the initial geometry when interpolating between 2 points. XYZ expected.",
    )
    parser.add_argument(
        "--fin",
        type=str,
        default="fin.xyz",
        help="Filename of the final geometry when interpolating between 2 points. XYZ expected.",
    )
    parser.add_argument(
        "--xyz",
        type=str,
        default="None",
        help="Name of a XYZ file with 2 or more replicas",
    )
    parser.add_argument(
        "--chk",
        type=str,
        default="None",
        help="Name of i-PI checkpoint file.",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default="interpolated_path.xyz",
        help="Name of the output xyz file (including extension).",
    )
    parser.add_argument(
        "-n",
        "-N",
        "--nbeads",
        required=True,
        help="Desired number of beads (including the endpoints)",
        type=int,
    )
    parser.add_argument(
        "-k",
        "--order",
        type=int,
        default=3,
        help="The order of interpolation. 1 and 3 are supported",
    )
    parser.add_argument(
        "--units",
        type=str,
        default="None",
        help="Units to output cell and positions.",
    )

    args = parser.parse_args()
    mode = args.mode
    k = args.order
    assert k in [1, 3]  # I don't know whether other spline orders work well.
    nbeads = args.nbeads
    assert nbeads >= 2  # 2 can be helpful also - to extract the endpoints quickly.

    if mode == "endpoints":
        input_ini = args.ini
        input_fin = args.fin
    elif mode == "xyz":
        input_xyz = args.xyz
    elif mode == "chk":
        input_chk = args.chk
    else:
        print("Error: cannot recognize mode %s." % mode)
        exit(-1)

    if mode == "endpoints":
        # Reading data from 2 geometries
        if os.path.exists(input_ini):
            with open(input_ini, "r") as fd_ini:
                inipos = []
                nframes = 0
                while True:
                    try:
                        ret = read_file("xyz", fd_ini)
                    except EOFError:  # finished reading file
                        break
                    inipos.append(ret["atoms"])
                    cell = ret["cell"]
                    if np.all(
                        cell.h == -np.eye(3)
                    ):  # default when cell not found in xyz file.
                        print("Error: no cell found in %s." % input_ini)
                        exit(-1)
                    nframes += 1
                    if nframes > 1:
                        print(
                            "Error: in 'endpoints' mode, xyz files should contain only one geometry each."
                        )
                        exit(-1)
            inipos = inipos[0]
            natoms = inipos.natoms
            with open(input_ini, "r") as fd_ini:
                rr = read_file_raw("xyz", fd_ini)
                (_, units, cell_units) = auto_units(rr["comment"])
        else:
            print("Error: cannot find {}.".format(input_ini))
            exit(-1)

        if os.path.exists(input_fin):
            with open(input_fin, "r") as fd_fin:
                finpos = []
                nframes = 0
                while True:
                    try:
                        ret = read_file("xyz", fd_fin)
                        finpos.append(ret["atoms"])
                        cell2 = ret["cell"]
                        if np.any(cell2.h != cell.h):
                            print(
                                "Error: initial and final structures have different cells:"
                            )
                            print(cell.h, cell2.h)
                            exit(-1)
                        nframes += 1
                        if nframes > 1:
                            print(
                                "Error: in 'endpoints' mode, xyz files should contain only one geometry each."
                            )
                            exit(-1)
                    except EOFError:  # finished reading file
                        break
            finpos = finpos[0]
            if natoms != finpos.natoms:
                print(
                    "Error: the number of atoms is different in ini and fin geometries."
                )
                exit(-1)
        else:
            print("Error: cannot find {}.".format(input_fin))
            exit(-1)

        prototype = inipos
        q = np.concatenate((inipos.q[np.newaxis, :], finpos.q[np.newaxis, :]), axis=0)
        nbeads_old = 2

    elif mode == "xyz":
        # Reading beads from a xyz file
        if os.path.exists(input_xyz):
            with open(input_xyz, "r") as fd_xyz:
                path = []
                nbeads_old = 0
                while True:
                    try:
                        ret = read_file("xyz", fd_xyz)
                        path.append(ret["atoms"])
                        cell = ret["cell"]
                        if np.all(
                            cell.h == -np.eye(3)
                        ):  # default when cell not found in xyz file.
                            print("Error: no cell found in %s." % input_xyz)
                            exit(-1)
                        nbeads_old += 1
                    except EOFError:  # fxyzshed reading files
                        break
                if nbeads_old < 2:
                    print("Error: at least 2 geometries are needed.")
                    exit(-1)

                natoms = path[0].natoms
            # Extract units information to write output with the same units as input
            with open(input_xyz, "r") as fd_xyz:
                rr = read_file_raw("xyz", fd_xyz)
                (_, units, cell_units) = auto_units(rr["comment"])
        else:
            print("Error: cannot find {}.".format(input_xyz))
            exit(-1)

        prototype = path[0]
        q = np.concatenate([x.q[np.newaxis, :] for x in path], axis=0)

    elif mode == "chk":
        # Reading beads from a checkpoint
        from ipi.engine.simulation import Simulation

        if os.path.exists(input_chk):
            simulation = Simulation.load_from_xml(
                input_chk, custom_verbosity="low", request_banner=False, read_only=True
            )
        else:
            print("Error: cannot find {}.".format(input_chk))
            sys.exit()
        cell = simulation.syslist[0].cell
        beads = simulation.syslist[0].motion.beads.copy()
        natoms = simulation.syslist[0].motion.beads.natoms
        nbeads_old = beads.nbeads
        if nbeads_old < 2:
            print("Error: at least 2 geometries are needed.")
            exit(-1)
        q = beads.q

        prototype = beads._blist[0]
        units = cell_units = "automatic"

    # Below this point, it shouldn't matter where the inputs came from.

    # Create new Atoms objects from the prorotype
    newpath = []
    for i in range(nbeads):
        newpath.append(prototype)

    new_q = spline_resample(q, nbeads_old, nbeads, k)

    if args.units != "None":
        units = cell_units = args.units

    out_fname = args.output
    with open(out_fname, "w") as fdout:
        for i in range(nbeads):
            newpath[i].q = new_q[i]
            print_file(
                mode="xyz",
                atoms=newpath[i],
                cell=cell,
                filedesc=fdout,
                title="Bead: %d " % i,
                key="positions",
                dimension="length",
                units=units,
                cell_units=cell_units,
            )

    print("The new path written to %s." % out_fname)

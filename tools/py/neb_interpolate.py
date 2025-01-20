#!/usr/bin/env python3
"""neb_interpolate.py
Creates interpolated path between 2 or more geometries. Beads are placed
equidistantly to a piecewise-linear path or a cubic spline formed by 2 or more
provided beads.

Arguments:
    --mode              - "endpoints" or "xyz" or "chk"
    -n, -N, --nbeads    - number of NEB replicas incl. endpoints
    -k                  - degree of spline: 1 or 3
    --units             - i-pi-supported distance units to write output geometries.
    -al, --activelist   - whitespace-separated list of active atoms to read from the command line
    -af, --activefile   - a file with the space or linebreak separated list of active atoms.
    -o, --output        - Name of the output xyz file (including extension). Default is "interpolated_path.xyz".
    --dry               - a boolean flag to switch off output into file(s).
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
from numpy.linalg import norm as npnorm

from ipi.utils.io import read_file, read_file_raw, print_file
from ipi.utils.io.io_units import auto_units
from ipi.utils.units import Constants

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

    if scipy is None:
        print("@spline_resample: spline interpolation requires scipy module.")
        raise (scipy_exception)

    if nbeads_new < 2:
        raise RuntimeError("@spline_resample: nbeads_new < 2.")

    # First, we calculate the current parameterization of the path
    # according to 3N-D Euclidean distances between adjacent beads.
    t = [0.0]
    current_t = 0.0
    print("@spline_resample: Cartesian 3N-dimensional distances between old beads:")
    for i in range(1, nbeads_old):
        dist = npnorm(q[i] - q[i - 1])
        print("\tfrom %3d to %3d : %.6f bohr" % (i - 1, i, dist))
        if dist <= 1e-2:
            print("Warning: two adjacent beads are very close.")
        current_t += dist
        t.append(current_t)

    t = np.array(t)
    t /= t[-1]
    print("@spline_resample: Spline parameter 't' before resampling:")
    print(t)
    print("t[i] - t[i-1]:")
    print(t - np.roll(t, 1))

    # New parameterization, equidistant by spline parameter
    new_t = np.arange(nbeads_new) / (nbeads_new - 1)
    print("New t:")
    print(new_t)

    # This reshapes q to a list of 3*Natoms sets consisting of nbeads_old numbers each,
    # describing trajectory of one Cartesian component of an atom over the path.
    components = q.T

    new_components = []
    for comp in components:
        # Interpolate the trajectory of one cartesian component by a spline.
        # bc_type='natural' imposes zero second derivative at the endpoints.
        if k == 3:
            spl = make_interp_spline(t, comp, k=k, bc_type="natural")
        elif k == 1:
            spl = make_interp_spline(t, comp, k=k)
        # Sample new parameter values uniformly from the old range.
        # Resample positions of this atom.
        new_components.append(splev(new_t, spl))

    # Reshape the path back to the beads shape
    new_q = np.array(new_components).T

    print("@spline_resample: Cartesian 3N-dimensional distances between NEW beads:")
    for i in range(1, nbeads_new):
        dist = npnorm(new_q[i] - new_q[i - 1])
        print("\tfrom %3d to %3d : %.6f bohr" % (i - 1, i, dist))

    return new_q


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""A script for interpolating a NEB path"""
    )
    parser.add_argument(
        "-m",
        "--mode",
        type=str,
        choices=["endpoints", "xyz", "chk"],
        default="endpoints",
        help="Mode: 'endpoints' to interpolate between 2 geometries provided separately,\n"
        "   'xyz' to interpolate between multiple geometries in a single file,\n"
        "   or 'chk' to load an existing checkpoint and resample path from it.",
    )
    parser.add_argument(
        "--ini",
        type=str,
        default="ini.xyz",
        help="Filename of the initial geometry when interpolating between 2 points (-m endpoints). XYZ expected.",
    )
    parser.add_argument(
        "--fin",
        type=str,
        default="fin.xyz",
        help="Filename of the final geometry when interpolating between 2 points (-m endpoints). XYZ expected.",
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
        help="Desired number of beads (including the endpoints).",
        type=int,
    )
    parser.add_argument(
        "-k",
        "--order",
        type=int,
        default=3,
        help="The order of interpolation. 1 (piecewise-linear) and 3 (cubic spline) are supported.",
    )
    parser.add_argument(
        "--units",
        type=str,
        default="None",
        help="Units to output cell and positions. Default is keeping units from inputs.",
    )
    parser.add_argument(
        "--activelist",
        "-al",
        nargs="+",
        type=int,
        default=[],
        help="List of atom indices to use. Positions of the others will be kept from bead 0, "
        "except one particular case when old and new number of beads coincide. "
        "Then, old positions will be kept. Atoms are counted from 0.",
    )
    parser.add_argument(
        "--activefile",
        "-af",
        type=str,
        default=None,
        help="A file with atom indices to use, the others will be kept from bead 0. "
        "Atoms are counted from 0.",
    )
    parser.add_argument(
        "--dry",
        type=bool,
        default=False,
        help="Flag to switch off all output to file(s)",
    )

    args = parser.parse_args()
    mode = args.mode
    k = args.order
    assert k in [1, 3]  # I don't know whether other spline orders work well.
    nbeads_new = args.nbeads
    assert nbeads_new >= 2  # 2 can be helpful also - to extract the endpoints quickly.
    assert (
        args.activelist == [] or args.activefile is None
    ), "Only one active list is allowed."

    if mode == "endpoints":
        # Reading data from 2 geometries
        input_ini = args.ini
        input_fin = args.fin
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
                    inimasses = ret["atoms"].m
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
                        finmasses = ret["atoms"].m
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
            print("Error: cannot find {} file.".format(input_fin))
            exit(-1)

        if np.all(inimasses == finmasses):
            masses = inimasses
        else:
            print(
                "Error: initial and final masses differ. "
                "Check that the input files have correct geometries and the same order of atoms."
            )
            exit(-1)

        prototype = inipos
        q = np.concatenate((inipos.q[np.newaxis, :], finpos.q[np.newaxis, :]), axis=0)
        nbeads_old = 2

    elif mode == "xyz":
        # Reading beads from a xyz file
        input_xyz = args.xyz
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
                masses = path[0].m
                # Check that the masses in all frames are the same
                for atoms in path:
                    if np.any(atoms.m != masses):
                        print(
                            "Error: masses differ across XYZ frames. "
                            "All frames should represent the same atomic system in the same order."
                        )
                        exit(-1)
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

        input_chk = args.chk
        if os.path.exists(input_chk):
            simulation = Simulation.load_from_xml(
                open(input_chk),
                custom_verbosity="low",
                request_banner=False,
                read_only=True,
            )
        else:
            print("Error: cannot find {}.".format(input_chk))
            sys.exit()
        cell = simulation.syslist[0].cell
        beads = simulation.syslist[0].motion.beads.clone()
        natoms = simulation.syslist[0].motion.beads.natoms
        masses = simulation.syslist[0].motion.beads.m  # KF: not sure about this line
        nbeads_old = beads.nbeads
        if nbeads_old < 2:
            print("Error: at least 2 geometries are needed.")
            exit(-1)
        q = beads.q

        prototype = beads._blist[0]
        units = cell_units = "automatic"

    # =================================================================
    # Below this point, it shouldn't matter where the inputs came from.
    # =================================================================
    print("Imported structures: %d geometries, %d atoms each." % (nbeads_old, natoms))
    print("Atomic masses:")
    masses /= Constants.amu
    for mm in masses:
        print(mm)

    # Parse the arguments for the list of active atoms
    activelist = None
    if args.activelist != []:
        activelist = args.activelist
    elif args.activefile is not None:
        with open(args.activefile, "r") as af:
            al = af.read().split()
            activelist = [int(w) for w in al]
    if activelist is not None:
        if any([a < 0 or a >= natoms for a in activelist]):
            print("Error: index < 0 or > natoms found in activelist.")
            print("Activelist: %s" % activelist)
            exit(-1)
        else:
            activelist = np.array(activelist)
            print("Activelist: %s" % activelist)
            # Mask to exclude fixed atoms from 3N-array. 0 are frozen, 1 are active.
            fixmask = np.zeros(3 * natoms, dtype=bool)
            fixmask[3 * activelist] = 1
            fixmask[3 * activelist + 1] = 1
            fixmask[3 * activelist + 2] = 1
    else:
        fixmask = np.ones(3 * natoms, dtype=bool)

    # Calculate distances in the mass-scaled space
    print("q.shape:")
    print(q.shape)
    sqrtm = np.sqrt(masses)
    sqrtm = np.column_stack((sqrtm, sqrtm, sqrtm)).flatten()  # We need the shape (3N)
    mass_scaled_q = q[:, fixmask] * sqrtm[fixmask]
    mass_scaled_t = 0
    print("Distances between the old beads in the mass-scaled space:")
    for i in range(1, nbeads_old):
        dist = npnorm(mass_scaled_q[i] - mass_scaled_q[i - 1])
        mass_scaled_t += dist
        print(
            "\tfrom %3d to %3d : %.6f bohr×√m, from start:  %6f"
            % (i - 1, i, dist, mass_scaled_t)
        )

    # Resample the path
    masked_new_q = spline_resample(q[:, fixmask], nbeads_old, nbeads_new, k)

    # Create new Atoms objects from the prorotype, and put there resampled coordinates
    # If nbeads remains the same, we keep fixed atoms as they were in the input beads.
    newpath = []
    if nbeads_new == nbeads_old:
        if activelist is not None:
            print(
                "nbeads_new = nbeads_old, "
                "therefore we keep fixed atoms in their original positions."
            )
        for i in range(nbeads_new):
            # copy() is necessary, otherwise they all point to the same obj.
            newpath.append(prototype.copy())
            newpath[i].q[:] = q[i]
            newpath[i].q[fixmask] = masked_new_q[i]
    else:  # If nbeads changes, we put fixed atoms from the bead 0 to all beads.
        if activelist is not None:
            print(
                "nbeads_new != nbeads_old, therefore we put fixed atoms from bead 0 to all beads."
            )
        for i in range(nbeads_new):
            newpath.append(prototype.copy())
            newpath[i].q[fixmask] = masked_new_q[i]

    if args.units != "None":
        units = cell_units = args.units

    out_fname = args.output
    if not args.dry:
        with open(out_fname, "w") as fdout:
            for i in range(nbeads_new):
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

        print("The new path is written to %s." % out_fname)
    else:
        print("Dry run, no files were written.")

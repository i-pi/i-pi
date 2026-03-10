"""Instanton_interpolation.py
Reads a hessian file  and/or a positions file (xyz format) and creates an interpolation
that can be used in a further instanton optimization with more beads

Syntax manual:    python  Instanton_interpolation.py -m -xyz <geometry file> -h <hessian file> -n <new-beads(half-polymer)>
Syntax chk:       python  Instanton_interpolation.py -chk  <checkpoint_file>  -n <new-beads(half-polymer)>

Example:   python  Instanton_interpolation.py  -xyz INSTANTON.xyz  -hess  INSTANTON.hess -n 30
           python  Instanton_interpolation.py  -chk RESTART -n 30

Relies on the infrastructure of i-pi, so the ipi package should
be installed in the Python module directory, or the i-pi
main directory must be added to the PYTHONPATH environment variable.

"""

# Y. Litman 2017

import os
import numpy as np
import sys
import argparse


# You can insert the i-pi path with the following lines.
# Uncomment them and adjust the ipi_path variable

# ipi_path='/home/litman/Yair/Instanton/I-PI-mc/i-pi-mc'

# if not (os.path.exists(ipi_path)):
#   print 'We can not find ipi in %s' %ipi_path
#   print 'Please correct the path'
#   sys.exit()
# sys.path.insert(0, ipi_path)

from ipi.utils.io import read_file, print_file
from ipi.utils.nmtransform import nm_rescale
from ipi.utils.units import unit_to_internal

# INPUT
parser = argparse.ArgumentParser(
    description="""Script for interpolate hessian and/or instanton geometry"""
)
parser.add_argument(
    "-m",
    "--manual",
    action="store_true",
    default=False,
    help="Boolean which decides between a checkpoint file or a manual entry.",
)
parser.add_argument(
    "-chk",
    "--checkpoint",
    type=str,
    default="None",
    help="Name of the instanton checkpoint file.",
)
parser.add_argument(
    "-xyz",
    "--xyz",
    type=str,
    default="None",
    help="Name of the instanton geometry file.",
)
parser.add_argument(
    "-hess", "--hessian", type=str, default="None", help="Name of the hessian file."
)
parser.add_argument(
    "-n",
    "--nbeadsNew",
    required=True,
    default=0,
    help="New number of beads (half polymer)",
    type=int,
)

args = parser.parse_args()
chk = args.checkpoint
input_geo = args.xyz
input_hess = args.hessian
nbeadsNew = args.nbeadsNew
manual = args.manual

if not manual:
    if chk == "None":
        print("Manual mode not specified and checkpoint file name not provided")
        sys.exit()
else:
    if input_geo == "None":
        print("Manual mode  specified and geometry file name not provided")
        sys.exit()

# OPEN AND READ   ###########################################################3


if input_geo != "None" or chk != "None":
    if manual:
        if os.path.exists(input_geo):
            ipos = open(input_geo, "r")
        else:
            print("We can't find {}".format(input_geo))
            sys.exit()

        pos = list()
        nbeads = 0
        while True:
            try:
                ret = read_file("xyz", ipos)
                pos.append(ret["atoms"])
                cell = ret["cell"]
                nbeads += 1
            except EOFError:  # finished reading files
                break
        ipos.close()

        natoms = pos[0].natoms
        atom = pos[0]
        # Compose the half ring polymer.
        q = np.vstack([i.q for i in pos])
    else:
        from ipi.engine.simulation import Simulation

        if os.path.exists(chk):
            simulation = Simulation.load_from_xml(
                open(chk), custom_verbosity="low", request_banner=False, read_only=True
            )
        else:
            print("We can't find {}".format(chk))
            sys.exit()
        cell = simulation.syslist[0].cell
        beads = simulation.syslist[0].motion.beads.clone()
        natoms = simulation.syslist[0].motion.beads.natoms
        nbeads = beads.nbeads
        q = beads.q
        atom = beads._blist[0]

    print(" ")
    print(
        "We have a half ring polymer made of {} beads and {} atoms.".format(
            nbeads, natoms
        )
    )
    print(
        "We will expand the ring polymer to get a half polymer of {} beads.".format(
            nbeadsNew
        )
    )

    # Make the rpc step (standard). It is better that open path in some corner cases.
    q2 = np.concatenate((q, np.flipud(q)), axis=0)  # Compose the full ring polymer.
    rpc = nm_rescale(2 * nbeads, 2 * nbeadsNew)
    new_q = rpc.b1tob2(q2)[0:nbeadsNew]

    # Make the rpc step (open path)
    # rpc = nm_rescale(nbeads, nbeadsNew, np.asarray(range(natoms)))
    # new_q = rpc.b1tob2(q)

    # Print
    out = open("new_instanton.xyz", "w")
    for i in range(nbeadsNew):
        atom.q = new_q[i] / unit_to_internal(
            "length", "angstrom", 1.0
        )  # Go back to angstrom
        print_file(
            "xyz", atom, cell, out, title="cell{atomic_unit}  Traj: positions{angstrom}"
        )
        # print_file("xyz",pos[0],cell,out,title='cell  }')
    out.close()

    print("The new Instanton geometry (half polymer) was generated")
    print("Check new_instanton.xyz")
    print("")
    print(
        "Don't forget to change the number of beads to the new value ({}) in your input file".format(
            nbeadsNew
        )
    )
    print("when starting your new simulation with an increased number of beads.")
    print("")

if input_hess != "None" or chk != "None":
    if manual:
        try:
            hess = open(input_hess, "r")
        except:
            print("We can't find {}".format(input_hess))
            sys.exit()
        h = np.zeros((natoms * 3) ** 2 * nbeads)
        aux = hess.readline().split()

        for i in range((natoms * 3) ** 2 * nbeads):
            h[i] = float(aux[i])
        h = h.reshape((natoms * 3, natoms * 3 * nbeads))
        hess.close()

    else:
        from ipi.engine.simulation import Simulation

        try:
            h = simulation.syslist[0].motion.optarrays["hessian"].copy()
        except:
            print("We don't have a hessian so there is nothing more to do")
            sys.exit()
        if np.linalg.norm(h) < 1e-13:
            print("We don't have a hessian so there is nothing more to do")
            sys.exit()

    print("The new hessian is {} x {}.".format(3 * natoms, natoms * 3 * nbeadsNew))
    out = open("new_hessian.dat", "w")

    print("Creating matrix... ")
    #    hessian = get_double_h(nbeads, natoms, h)

    hessian = h
    size0 = natoms * 3

    # # We use open path RPC
    # size1 = size0 * nbeads
    # size2 = size0 * nbeadsNew
    # new_h = np.zeros([size0, size2])
    # rpc = nm_rescale(nbeads, nbeadsNew, np.asarray(range(1)))
    # new_q = rpc.b1tob2(q)

    # # Compose the full ring polymer.
    size1 = size0 * (2 * nbeads)
    size2 = size0 * (2 * nbeadsNew)
    new_h = np.zeros([size0, size2])
    q2 = np.concatenate((q, np.flipud(q)), axis=0)  # Compose the full ring polymer.
    rpc = nm_rescale(2 * nbeads, 2 * nbeadsNew)
    new_q = rpc.b1tob2(q2)[0:nbeadsNew]

    for i in range(size0):
        for j in range(size0):
            h = np.array([])
            for n in range(nbeads):
                h = np.append(h, hessian[i, j + size0 * n])
            #           h3 = np.concatenate((h, h, h), axis=0).reshape((h.size, 3), order='F') # Open path expect three coordinates per atom
            #           diag = rpc.b1tob2(h3)[:, 0] # Open path
            h2 = np.concatenate(
                (h, np.flipud(h)), axis=0
            )  # Compose the full ring polymer.
            diag = rpc.b1tob2(h2)
            new_h[i, j:size2:size0] += diag

    new_h_half = new_h[:, 0 : size2 // 2]
    np.savetxt(out, new_h_half.reshape(1, new_h_half.size))
    # new_h_half = new_h[:, 0:size2 / 2]
    # np.savetxt(out, new_h.reshape(1, new_h.size))

    print("The new physical Hessian (half polymer) was generated")
    print("Check new_hessian.dat")
    print("")
    print("Remeber to adapt/add the following line in your input file:")
    print("")
    print(
        " <hessian mode='file' shape='({}, {})' >hessian.dat</hessian>".format(
            3 * natoms, natoms * 3 * nbeadsNew
        )
    )
    print("")

sys.exit()

#!/usr/bin/env python3

"""remdsort.py

Relies on the infrastructure of i-pi, so the ipi package should
be installed in the Python module directory, or the i-pi
main directory must be added to the PYTHONPATH environment variable.

Post-processes the output of a replica exchange simulation and
re-orders the outputs so that they correspond to the different
ensembles rather than to the time series of one of
the replicas exchanging ensembles over time.

It should be run in the same directory as where i-pi was (or is being)
run, and simply fetches all information from the simulation input file.
Will create a series of PTindex_* files, each corresponding to the
data for replica 'index'.

Syntax:
   remdsort.py inputfile.xml
"""


import sys
from copy import deepcopy
import numpy as np
from ipi.utils.messages import verbosity
from ipi.engine.outputs import *
from ipi.engine.properties import getkey
from ipi.inputs.simulation import InputSimulation
from ipi.utils.io.inputs import io_xml


def main(inputfile, prefix="SRT_"):
    verbosity.level = "low"
    # opens & parses the input file
    ifile = open(inputfile, "r")
    xmlrestart = io_xml.xml_parse_file(ifile)  # Parses the file.
    ifile.close()

    # ugly hack to remove ffplumed objects to avoid messing up with plumed output files
    newfields = [f for f in xmlrestart.fields[0][1].fields if f[0] != "ffplumed"]
    xmlrestart.fields[0][1].fields = newfields

    isimul = InputSimulation()
    isimul.parse(xmlrestart.fields[0][1])

    simul = isimul.fetch()
    swapfile = ""
    if simul.smotion is None or (
        simul.smotion.mode != "remd" and simul.smotion.mode != "multi"
    ):
        raise ValueError("Simulation does not look like a parallel tempering one.")
    else:
        if simul.smotion.mode == "remd":
            swapfile = simul.smotion.swapfile
        else:
            for sm in simul.smotion.mlist:
                if sm.mode == "remd":
                    if swapfile != "":
                        raise ValueError(
                            "I'm not equipped to deal with multiple REMD outputs, sorry"
                        )
                    swapfile = sm.swapfile
        if swapfile == "":
            raise ValueError(
                "Could not determine the REMD swapfile name. \
                 Sorry, you'll have to look carefully at your inputs."
            )

    # reconstructs the list of the property and trajectory files that have been output
    # and that should be re-ordered
    lprop = []  # list of property files
    ltraj = []  # list of trajectory files
    nsys = len(simul.syslist)
    for o in simul.outtemplate:
        o = deepcopy(o)  # avoids overwriting the actual filename
        if simul.outtemplate.prefix != "":
            o.filename = simul.outtemplate.prefix + "." + o.filename
        if (
            type(o) is CheckpointOutput
        ):  # properties and trajectories are output per system
            pass
        elif type(o) is PropertyOutput:
            nprop = []
            isys = 0
            for s in simul.syslist:  # create multiple copies
                if s.prefix != "":
                    filename = s.prefix + "_" + o.filename
                else:
                    filename = o.filename
                ofilename = prefix + filename
                nprop.append(
                    {
                        "filename": filename,
                        "ofilename": ofilename,
                        "stride": o.stride,
                        "ifile": open(filename, "r"),
                        "ofile": open(ofilename, "w"),
                    }
                )
                isys += 1
            lprop.append(nprop)
        elif (
            type(o) is TrajectoryOutput
        ):  # trajectories are more complex, as some have per-bead output
            if getkey(o.what) in [
                "positions",
                "velocities",
                "forces",
                "extras",
            ]:  # multiple beads
                nbeads = simul.syslist[0].beads.nbeads
                for b in range(nbeads):
                    ntraj = []
                    isys = 0
                    # zero-padded bead number
                    padb = (
                        "%0" + str(int(1 + np.floor(np.log(nbeads) / np.log(10)))) + "d"
                    ) % (b)
                    for s in simul.syslist:
                        if s.prefix != "":
                            filename = s.prefix + "_" + o.filename
                        else:
                            filename = o.filename
                        ofilename = prefix + filename
                        if o.ibead < 0 or o.ibead == b:
                            if getkey(o.what) == "extras":
                                filename = filename + "_" + padb
                                ofilename = ofilename + "_" + padb
                                # Sets format of extras as None
                                ntraj.append(
                                    {
                                        "filename": filename,
                                        "format": None,
                                        "ofilename": ofilename,
                                        "stride": o.stride,
                                        "ifile": open(filename, "r"),
                                        "ofile": open(ofilename, "w"),
                                    }
                                )
                            else:
                                extension = o.format
                                if extension == "ase":  # special case
                                    extension = "extxyz"
                                filename = filename + "_" + padb + "." + extension
                                ofilename = ofilename + "_" + padb + "." + extension
                                ntraj.append(
                                    {
                                        "filename": filename,
                                        "format": o.format,
                                        "ofilename": ofilename,
                                        "stride": o.stride,
                                        "ifile": open(filename, "r"),
                                        "ofile": open(ofilename, "w"),
                                    }
                                )
                        isys += 1
                    if ntraj != []:
                        ltraj.append(ntraj)

            else:
                ntraj = []
                isys = 0
                for s in simul.syslist:  # create multiple copies
                    if s.prefix != "":
                        filename = s.prefix + "_" + o.filename
                    else:
                        filename = o.filename
                    extension = o.format
                    if extension == "ase":  # special case
                        extension = "extxyz"
                    filename = filename + "." + extension
                    ofilename = prefix + filename
                    ntraj.append(
                        {
                            "filename": filename,
                            "format": o.format,
                            "ofilename": ofilename,
                            "stride": o.stride,
                            "ifile": open(filename, "r"),
                            "ofile": open(ofilename, "w"),
                        }
                    )

                    isys += 1
                ltraj.append(ntraj)

    ptfile = open(simul.outtemplate.prefix + "." + swapfile, "r")

    # now reads files one frame at a time,
    # and re-direct output to the appropriate location

    line = ptfile.readline().split()
    irep = list(range(nsys))  # Could this be harmful?
    step = 0
    while True:
        # reads one line from index file
        try:
            for prop in lprop:
                for isys in range(nsys):
                    sprop = prop[isys]
                    if step % sprop["stride"] == 0:  # property transfer
                        iline = sprop["ifile"].readline()
                        if len(iline) == 0:
                            raise EOFError  # useful if line is blank
                        while iline[0] == "#":  # fast forward if line is a comment
                            prop[irep[isys]]["ofile"].write(iline)
                            iline = sprop["ifile"].readline()
                        prop[irep[isys]]["ofile"].write(iline)
            for traj in ltraj:
                for isys in range(nsys):
                    straj = traj[isys]
                    if step % straj["stride"] == 0:  # property transfer
                        # reads one frame from the input file
                        ibuffer = []
                        if straj["format"] is None:
                            ibuffer.append(straj["ifile"].readline())
                            ibuffer.append(straj["ifile"].readline())
                            traj[irep[isys]]["ofile"].write("".join(ibuffer))
                        if straj["format"] in ["xyz", "ase"]:
                            iline = straj["ifile"].readline()
                            nat = int(iline)
                            ibuffer.append(iline)
                            ibuffer.append(straj["ifile"].readline())
                            for i in range(nat):
                                ibuffer.append(straj["ifile"].readline())
                            traj[irep[isys]]["ofile"].write("".join(ibuffer))
                        elif straj["format"] == "pdb":
                            iline = straj["ifile"].readline()
                            while iline.strip() != "" and iline.strip() != "END":
                                ibuffer.append(iline)
                                iline = straj["ifile"].readline()
                            ibuffer.append(iline)
                            traj[irep[isys]]["ofile"].write("".join(ibuffer))
        except EOFError:
            break

        if len(line) > 0 and step == int(line[0]):
            irep = [int(i) for i in line[1:]]
            line = ptfile.readline()
            line = line.split()

        step += 1


if __name__ == "__main__":
    main(*sys.argv[1:])

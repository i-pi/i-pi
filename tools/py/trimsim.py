#!/usr/bin/env python3

"""trimsim.py

Relies on the infrastructure of i-pi, so the ipi package should
be installed in the Python module directory, or the i-pi
main directory must be added to the PYTHONPATH environment variable.

Cuts short the output of a previous i-pi simulation, up to the
step indicated in the <step> field of the input file.
This is useful to restart a simulation that crashed.

It should be run in the same dyrectory as where i-pi was (or is being)
run, and simply fetches all information from the simulation input file.
One should also specify a directory name in which the trimmed files
will be output.

Syntax:
   trimsim.py inputfile.xml
"""


import sys
import os
import numpy as np
from copy import deepcopy
from ipi.engine.outputs import *
from ipi.engine.properties import getkey
from ipi.inputs.simulation import InputSimulation
from ipi.utils.io.inputs import io_xml


def main(inputfile, outdir="trim"):
    # opens & parses the input file
    ifile = open(inputfile, "r")
    xmlrestart = io_xml.xml_parse_file(ifile)  # Parses the file.
    ifile.close()

    isimul = InputSimulation()
    isimul.parse(xmlrestart.fields[0][1])

    simul = isimul.fetch()
    trimstep = isimul.step.fetch()

    os.makedirs(outdir)

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
                ofilename = outdir + "/" + filename
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
                        ofilename = outdir + "/" + filename
                        if o.ibead < 0 or o.ibead == b:
                            if getkey(o.what) == "extras":
                                filename = filename + "_" + padb
                                ofilename = ofilename + "_" + padb
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
                                filename = filename + "_" + padb + "." + o.format
                                ofilename = ofilename + "_" + padb + "." + o.format
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
                    filename = filename + "." + o.format
                    ofilename = outdir + "/" + filename
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

    ptfile = None
    # wtefile = None
    if hasattr(simul.smotion, "swapfile") and os.path.isfile(
        simul.outtemplate.prefix + "." + simul.smotion.swapfile
    ):
        ptfile = open(simul.outtemplate.prefix + "." + simul.smotion.swapfile, "r")
        optfile = open(
            outdir + "/" + simul.outtemplate.prefix + "." + simul.smotion.swapfile, "w"
        )
    # do not know if this is redudant, please uncomment if it is not
    # if os.path.isfile("PARAWTE"):
    #    wtefile = open("PARAWTE", "r")
    #    owtefile = open(outdir + "/PARAWTE", "w")

    # First reads the swap file
    if ptfile is not None:
        while True:
            try:
                line = ptfile.readline()
                step = int(line.split()[0])
                if step < trimstep:
                    optfile.write(line)
                else:
                    break
            except IndexError:
                break

    # now reads files one frame at a time, and re-direct output to the appropriate location
    for step in range(trimstep + 1):
        try:
            for prop in lprop:
                for isys in range(nsys):
                    sprop = prop[isys]
                    if step % sprop["stride"] == 0:  # property transfer
                        iline = sprop["ifile"].readline()
                        while iline[0] == "#":  # fast forward if line is a comment
                            prop[isys]["ofile"].write(iline)
                            iline = sprop["ifile"].readline()
                        prop[isys]["ofile"].write(iline)

            for traj in ltraj:
                for isys in range(nsys):
                    straj = traj[isys]
                    if step % straj["stride"] == 0:  # property transfer
                        # reads one frame from the input file
                        ibuffer = []
                        if straj["format"] is None:
                            ibuffer.append(straj["ifile"].readline())
                            ibuffer.append(straj["ifile"].readline())
                            traj[isys]["ofile"].write("".join(ibuffer))
                        if straj["format"] in ["xyz", "ase"]:
                            iline = straj["ifile"].readline()
                            nat = int(iline)
                            ibuffer.append(iline)
                            ibuffer.append(straj["ifile"].readline())
                            for i in range(nat):
                                ibuffer.append(straj["ifile"].readline())
                            traj[isys]["ofile"].write("".join(ibuffer))
                        elif straj["format"] == "pdb":
                            iline = straj["ifile"].readline()
                            while iline.strip() != "" and iline.strip() != "END":
                                ibuffer.append(iline)
                                iline = straj["ifile"].readline()
                            ibuffer.append(iline)
                            traj[isys]["ofile"].write("".join(ibuffer))
        except EOFError:
            break


if __name__ == "__main__":
    main(*sys.argv[1:])

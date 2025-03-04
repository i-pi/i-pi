#!/usr/bin/env python3

"""paraweights.py

Relies on the infrastructure of i-pi, so the ipi package should
be installed in the Python module directory, or the i-pi
main directory must be added to the PYTHONPATH environment variable.

Post-processes the output of a parallel-tempering simulation and
computes - for each output file that contains potential energy
information - the log-weights that should be given to each
line to re-weight the ensemble at a given target temperature.
For each file, it also prints out the "global" weighing factor that
should be used to combine different PT replicas, based on a simple
estimate of the error based on
Ceriotti, Brain, Riordan, Manolopoulos, Proc. Royal Soc. A 468 (2011)
and the assumption that the observable and the weight are un-correlated.
In computing the global weight, by default a certain number of steps
will be ignored.

It should be run in the same dyrectory as where i-pi was (or is being)
run, and simply fetches all information from the simulation input file.
Will create a series of (prefix)index_* files, each corresponding to the
data for replica 'index', and a (prefix)WEIGHTS file containing the
trajectory weights.

Syntax:
   paraweights.py inputfile.xml [prefix] [temperature(K)] [skip]
"""


import sys
import re
import numpy as np
from ipi.utils.messages import verbosity, banner
from ipi.engine.outputs import *
from ipi.inputs.simulation import InputSimulation
from ipi.utils.io.inputs import io_xml
from ipi.utils.units import unit_to_internal
from ipi.utils.mathtools import logsumlog


def main(inputfile, prefix="PTW-", ttemp="300.0", skip="2000"):
    txtemp = ttemp
    ttemp = unit_to_internal("energy", "kelvin", float(ttemp))
    skip = int(skip)

    # opens & parses the input file
    ifile = open(inputfile, "r")
    verbosity.level = "quiet"
    verbosity.lock = True
    xmlrestart = io_xml.xml_parse_file(ifile)  # Parses the file.
    ifile.close()

    isimul = InputSimulation()
    isimul.parse(xmlrestart.fields[0][1])
    verbosity.level = "quiet"
    banner()
    print(
        "# Printing out temperature re-weighing factors for a parallel tempering simulation"
    )
    simul = isimul.fetch()

    if simul.mode != "paratemp":
        raise ValueError("Simulation does not look like a parallel tempering one.")

    # reconstructs the list of the property and trajectory files that have been output
    # and that should be re-ordered
    lprop = []  # list of property files
    # ltraj = []  # list of trajectory files
    tlist = simul.paratemp.temp_list
    nsys = len(simul.syslist)
    for o in simul.outtemplate:
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
                ofilename = (
                    prefix
                    + (
                        (
                            "%0"
                            + str(
                                int(
                                    1
                                    + np.floor(np.log(len(simul.syslist)) / np.log(10))
                                )
                            )
                            + "d"
                        )
                        % (isys)
                    )
                    + "_"
                    + o.filename
                )
                nprop.append(
                    {
                        "filename": filename,
                        "ofilename": ofilename,
                        "stride": o.stride,
                        "ifile": open(filename, "r"),
                        "ofile": None,
                    }
                )
                isys += 1
            lprop.append(nprop)

    ptfile = open("PARATEMP", "r")

    # these are variables used to compute the weighting factors
    tprops = []
    vfields = []
    refpots = []
    vunits = []
    nw = []
    tw = []
    tw2 = []

    repot = re.compile(" ([0-9]*) *--> potential")
    reunit = re.compile("{(.*)}")

    # now reads files one frame at a time, and re-direct output to the appropriate location
    irep = np.zeros(nsys, int)
    while True:
        # reads one line from PARATEMP index file
        line = ptfile.readline()
        line = line.split()

        try:
            if len(line) == 0:
                raise EOFError

            step = int(line[0])
            irep[:] = line[1:]

            wk = 0
            for prop in lprop:
                for isys in range(nsys):
                    sprop = prop[isys]
                    if step % sprop["stride"] == 0:  # property transfer
                        iline = sprop["ifile"].readline()
                        if len(iline) == 0:
                            raise EOFError
                        while iline[0] == "#":  # fast forward if line is a comment
                            # checks if we have one single file with potential energies
                            rm = repot.search(iline)
                            if not (rm is None) and not (prop in tprops):
                                tprops.append(prop)
                                for p in prop:
                                    p["ofile"] = open(p["ofilename"], "w")
                                    p["ofile"].write(
                                        "# column   1     --> ptlogweight: ln of re-weighing factor with target temperature %s K\n"
                                        % (txtemp)
                                    )
                                vfields.append(int(rm.group(1)) - 1)
                                refpots.append(np.zeros(nsys))
                                nw.append(np.zeros(nsys))
                                tw.append(np.zeros(nsys))
                                tw2.append(np.zeros(nsys))
                                rm = reunit.search(iline)
                                if rm:
                                    vunits.append(rm.group(1))
                                else:
                                    vunits.append("atomic_unit")
                            iline = sprop["ifile"].readline()
                        if prop in tprops:  # do temperature weighing
                            pot = unit_to_internal(
                                "energy", vunits[wk], float(iline.split()[vfields[wk]])
                            )
                            ir = irep[isys]
                            if nw[wk][ir] == 0:
                                refpots[wk][
                                    ir
                                ] = pot  # picks the first value as a reference potential to avoid insane absolute values of the logweight
                            temp = tlist[ir]
                            lw = (pot - refpots[wk][ir]) * (1 / temp - 1 / ttemp)
                            if (
                                step > skip
                            ):  # computes trajectory weights avoiding the initial - possibly insane - values
                                if nw[wk][ir] == 0:
                                    tw[wk][ir] = lw
                                    tw2[wk][ir] = lw
                                else:
                                    tw[wk][ir] = logsumlog((tw[wk][ir], 1), (lw, 1))[0]
                                    tw2[wk][ir] = logsumlog(
                                        (tw2[wk][ir], 1), (2 * lw, 1)
                                    )[0]
                                nw[wk][ir] += 1
                            prop[ir]["ofile"].write("%15.7e\n" % (lw))
                            if isys == nsys - 1:
                                wk += 1
        except EOFError:
            # print out trajectory weights based on PRSA 2011, assuming that observables and weights are weakly correlated
            wk = 0
            fpw = open(prefix + "WEIGHTS", "w")
            fpw.write("# Global trajectory weights for temperature %s K\n" % (txtemp))
            fpw.write(
                "# Please cite M. Ceriotti, G. A. Brain, O. Riordan, D.E. Manolopoulos, "
                + "The inefficiency of re-weighted sampling and the curse of system size in high-order path integration. "
                + "Proceedings of the Royal Society A, 468(2137), 2-17  (2011) \n"
            )
            for prop in lprop:
                if prop in tprops:
                    for ir in range(nsys):
                        fpw.write(
                            "%s   %15.7e \n"
                            % (
                                prop[ir]["ofilename"],
                                1.0
                                / (np.exp(tw2[wk][ir] - 2 * tw[wk][ir]) * nw[wk][ir]),
                            )
                        )
                    wk += 1
            break


if __name__ == "__main__":
    main(*sys.argv[1:])

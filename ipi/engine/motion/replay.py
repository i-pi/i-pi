"""TODO

"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import time
import os
from fnmatch import fnmatch

from ipi.engine.motion import Motion
from ipi.utils.softexit import softexit
from ipi.utils.io import read_file, read_file_raw
from ipi.utils.io.inputs.io_xml import xml_parse_file
from ipi.utils.units import unit_to_internal
from ipi.utils.messages import verbosity, info


__all__ = ["Replay"]


class Replay(Motion):
    """Calculator object that just loads snapshots from an external file in sequence.

    Has the relevant conserved quantity and normal mode propagator for the
    constant energy ensemble. Note that a temperature of some kind must be
    defined so that the spring potential can be calculated.

    Attributes:
        intraj: The input trajectory file.
        ptime: The time taken in updating the velocities.
        qtime: The time taken in updating the positions.
        ttime: The time taken in applying the thermostat steps.

    Depend objects:
        None really meaningful.
    """

    def __init__(self, fixcom=False, fixatoms=None, intraj=None):
        """Initialises Replay.

        Args:
           dt: The simulation timestep.
           temp: The system temperature.
           fixcom: An optional boolean which decides whether the centre of mass
              motion will be constrained or not. Defaults to False.
           intraj: The input trajectory file.
        """

        super(Replay, self).__init__(fixcom=fixcom, fixatoms=fixatoms)
        if intraj is None:
            raise ValueError(
                "Must provide an initialized InitFile object to read trajectory from"
            )
        self.intraj = intraj
        if intraj.mode == "manual":
            raise ValueError(
                "Replay can only read from PDB or XYZ files -- or a single frame from a CHK file"
            )
        # Posibility to read beads from separate XYZ files by a wildcard
        if any(char in self.intraj.value for char in "*?[]"):
            infilelist = []
            for file in sorted(os.listdir(".")):
                if fnmatch(file, self.intraj.value):
                    infilelist.append(file)
            # determine bead numbers in input files
            bead_map_list = []
            for file in infilelist:
                fdin = open(file, "r")
                rr = read_file_raw(self.intraj.mode, fdin)
                metainfo = rr["comment"].split()
                for i, word in enumerate(metainfo):
                    if word == "Bead:":
                        bead_map_list.append(int(metainfo[i + 1]))
                fdin.close()
            # check that beads are continuous (no files missing)
            if sorted(bead_map_list) != list(range(len(bead_map_list))):
                info(
                    "ATTENTION: Provided trajectory files have non-sequential "
                    "range of bead indices.\n"
                    "\tIndices found: %s\n"
                    "\tMake sure that the wildcard does what it's supposed to do."
                    % str(bead_map_list),
                    verbosity.low,
                )
            # sort the list of files according to their bead indices
            infilelist_sorted, _ = zip(
                *sorted(zip(infilelist, bead_map_list), key=lambda t: t[1])
            )
            self.rfile = [open(f, "r") for f in infilelist_sorted]
        else:  # no wildcard
            self.rfile = open(self.intraj.value, "r")
        self.rstep = 0

    def step(self, step=None):
        """Does one replay time step."""

        self.ptime = 0.0
        self.ttime = 0.0
        self.qtime = -time.time()

        # If wildcard is used, check that it is consistent with Nbeads
        wildcard_used = False
        if any(char in self.intraj.value for char in "*?[]"):
            wildcard_used = True
            if len(self.rfile) != len(self.beads):
                info(
                    "Error: if a wildcard is used for replay, then "
                    "the number of files should be equal to the number of beads.",
                    verbosity.low,
                )
                softexit.trigger(status="bad", message=" # Error in replay input.")
        while True:
            self.rstep += 1
            try:
                if self.intraj.mode == "xyz":
                    for bindex, b in enumerate(self.beads):
                        if wildcard_used:
                            myframe = read_file("xyz", self.rfile[bindex])
                        else:
                            myframe = read_file("xyz", self.rfile)
                        myatoms = myframe["atoms"]
                        mycell = myframe["cell"]
                        myatoms.q *= unit_to_internal("length", self.intraj.units, 1.0)
                        mycell.h *= unit_to_internal("length", self.intraj.units, 1.0)
                        b.q[:] = myatoms.q
                elif self.intraj.mode == "pdb":
                    for bindex, b in enumerate(self.beads):
                        if wildcard_used:
                            myframe = read_file("pdb", self.rfile[bindex])
                        else:
                            myframe = read_file("pdb", self.rfile)
                        myatoms = myframe["atoms"]
                        mycell = myframe["cell"]
                        myatoms.q *= unit_to_internal("length", self.intraj.units, 1.0)
                        mycell.h *= unit_to_internal("length", self.intraj.units, 1.0)
                        b.q[:] = myatoms.q
                elif self.intraj.mode == "ase":
                    for bindex, b in enumerate(self.beads):
                        if wildcard_used:
                            myframe = read_file(
                                "ase", self.rfile[bindex], dimension="length"
                            )
                        else:
                            myframe = read_file("ase", self.rfile, dimension="length")
                        myatoms = myframe["atoms"]
                        mycell = myframe["cell"]
                        b.q[:] = myatoms.q

                elif self.intraj.mode == "chk" or self.intraj.mode == "checkpoint":
                    # TODO: Adapt the new `Simulation.load_from_xml`?
                    # reads configuration from a checkpoint file
                    xmlchk = xml_parse_file(self.rfile)  # Parses the file.

                    from ipi.inputs.simulation import InputSimulation

                    simchk = InputSimulation()
                    simchk.parse(xmlchk.fields[0][1])
                    mycell = simchk.cell.fetch()
                    mybeads = simchk.beads.fetch()
                    self.beads.q[:] = mybeads.q
                    softexit.trigger(
                        status="success", message=" # Read single checkpoint"
                    )
                # do not assign cell if it contains an invalid value (typically missing cell in the input)
                if mycell.V > 0:
                    self.cell.h[:] = mycell.h
            except EOFError:
                softexit.trigger(
                    status="success", message=" # Finished reading re-run trajectory"
                )
            if (step is None) or (self.rstep > step):
                break

        self.qtime += time.time()

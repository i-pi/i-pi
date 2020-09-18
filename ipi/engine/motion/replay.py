"""TODO

"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import time

from ipi.engine.motion import Motion
from ipi.utils.softexit import softexit
from ipi.utils.io import read_file
from ipi.utils.io.inputs.io_xml import xml_parse_file
from ipi.utils.units import unit_to_internal


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
        self.rfile = open(self.intraj.value, "r")
        self.rstep = 0

    def step(self, step=None):
        """Does one replay time step."""

        self.ptime = 0.0
        self.ttime = 0.0
        self.qtime = -time.time()

        while True:
            self.rstep += 1
            try:
                if self.intraj.mode == "xyz":
                    for b in self.beads:
                        myframe = read_file("xyz", self.rfile)
                        myatoms = myframe["atoms"]
                        mycell = myframe["cell"]
                        myatoms.q *= unit_to_internal("length", self.intraj.units, 1.0)
                        mycell.h *= unit_to_internal("length", self.intraj.units, 1.0)
                        b.q[:] = myatoms.q
                elif self.intraj.mode == "pdb":
                    for b in self.beads:
                        myatoms, mycell = read_file("pdb", self.rfile)
                        myatoms.q *= unit_to_internal("length", self.intraj.units, 1.0)
                        mycell.h *= unit_to_internal("length", self.intraj.units, 1.0)
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
                    softexit.trigger(" # Read single checkpoint")
                # do not assign cell if it contains an invalid value (typically missing cell in the input)
                if mycell.V > 0:
                    self.cell.h[:] = mycell.h
            except EOFError:
                softexit.trigger(" # Finished reading re-run trajectory")
            if (step is None) or (self.rstep > step):
                break

        self.qtime += time.time()

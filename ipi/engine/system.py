"""Contains the class that deals with storing the state of a physical system.

Contains code used to hold the information which represents the state of
a system, including the particle positions and momenta, and the
forcefields which govern the interaction potential.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import threading

from ipi.utils.depend import dpipe
from ipi.utils.messages import verbosity, info
from ipi.engine.forces import Forces
from ipi.engine.properties import Properties, Trajectories


__all__ = ["System"]


class System:
    """Physical system object.

    Contains all the physical information. Also handles stepping and output.

    Attributes:
       beads: A beads object giving the atom positions.
       cell: A cell object giving the system box.
       fcomp: A list of force components that must act on each replica
       forces: A Forces object that actually compute energy and forces
       ensemble: An ensemble object giving the objects necessary for producing
          the correct ensemble.
       outputs: A list of output objects that should be printed during the run
       nm:  A helper object dealing with normal modes transformation
       properties: A property object for dealing with property output.
       trajs: A trajectory object for dealing with trajectory output.
       init: A class to deal with initializing the system.
       simul: The parent simulation object.
    """

    def __init__(
        self, init, beads, nm, cell, fcomponents, ensemble=None, motion=None, prefix=""
    ):
        """Initialises System class.

        Args:
           init: A class to deal with initializing the system.
           beads: A beads object giving the atom positions.
           cell: A cell object giving the system box.
           fcomponents: A list of force components that are active for each
              replica of the system.
           bcomponents: A list of force components that are considered as bias, and act on each
              replica of the system.
           ensemble: An ensemble object giving the objects necessary for
              producing the correct ensemble.
           nm: A class dealing with path NM operations.
           prefix: A string used to differentiate the output files of different
              systems.
        """

        info(" # Initializing system object ", verbosity.low)
        self.prefix = prefix
        self.init = init
        self.ensemble = ensemble
        self.motion = motion
        self.beads = beads
        self.cell = cell
        self.nm = nm

        self.fcomp = fcomponents
        self.forces = Forces()

        self.properties = Properties()
        self.trajs = Trajectories()

    def bind(self, simul):
        """Calls the bind routines for all the objects in the system."""

        self.simul = simul  # keeps a handle to the parent simulation object

        # binds important computation engines
        info(" # Binding the forces ", verbosity.low)
        self.forces.bind(
            self.beads,
            self.cell,
            self.fcomp,
            self.simul.fflist,
            open_paths=self.nm.open_paths,
            output_maker=simul.output_maker,
        )
        self.nm.bind(self.ensemble, self.motion, beads=self.beads, forces=self.forces)
        self.ensemble.bind(
            self.beads,
            self.nm,
            self.cell,
            self.forces,
            self.simul.fflist,
            simul.output_maker,
        )
        self.motion.bind(
            self.ensemble,
            self.beads,
            self.nm,
            self.cell,
            self.forces,
            self.prng,
            simul.output_maker,
        )

        dpipe(self.nm._omegan2, self.forces._omegan2)

        self.init.init_stage2(self)

        # binds output management objects
        self._propertylock = threading.Lock()
        self.properties.bind(self)
        self.trajs.bind(self)

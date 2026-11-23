"""TODO"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

import numpy as np

from ipi.utils.depend import depend_value, dproperties


class Motion:
    """Base motion calculation class.

    Gives the standard methods and attributes needed in all the
    motion calculation classes.

    Attributes:
        beads: A beads object giving the atoms positions.
        cell: A cell object giving the system box.
        forces: A forces object giving the virial and the forces acting on
            each bead.
        fixcom: A boolean which decides whether the centre of mass
            motion will be constrained or not.
        fixatoms_dof: A list of degrees of freedom  that should be held fixed to their
            initial positions.
        fixcentroid: Whether every atomic centroid is held at its initial position.
        fixbeads: Zero-based replica indices held at their initial positions.

    Depend objects:
        none
    """

    def __init__(
        self, fixcom=False, fixatoms_dof=None, fixcentroid=False, fixbeads=None
    ):
        """Initialises Motion object.

        Args:
           fixcom: An optional boolean which decides whether the centre of mass
              motion will be constrained or not. Defaults to False.
           fixatoms_dof: A list of degrees of freedom  that should be held fixed to their
              initial positions.
        """

        self._dt = depend_value(name="dt", value=0.0)
        self.fixcom = fixcom
        if fixatoms_dof is None:
            self.fixatoms_dof = np.zeros(0, int)
        else:
            self.fixatoms_dof = fixatoms_dof
        self.fixcentroid = bool(fixcentroid)
        self.fixbeads = np.unique(np.asarray([] if fixbeads is None else fixbeads, int))
        if self.fixcentroid and len(self.fixbeads):
            raise ValueError("fixcentroid and fixbeads cannot be used together")
        if self.fixcentroid or len(self.fixbeads):
            self.fixcom = False

        self.beads = self.cell = self.forces = self.prng = self.nm = self.enstype = None
        self.finished = False
        self.exit_status = None
        self.exit_message = ""

    def finish(self, status="success", message=""):
        """Requests a clean simulation exit from the main simulation loop."""

        if status not in ("success", "restartable", "bad"):
            raise ValueError("Unknown motion exit status: " + str(status))
        self.finished = True
        self.exit_status = status
        self.exit_message = message

    def bind(self, ens, beads, nm, cell, bforce, prng, omaker):
        """Binds beads, cell, bforce, and prng to the calculator.

        This takes a beads object, a cell object, a forcefield object and a
        random number generator object and makes them members of the atom motion calculator.
        It also then creates the objects that will hold the data needed in the
        ensemble algorithms and the dependency network. Note that the conserved
        quantity is defined in the init, but as each ensemble has a different
        conserved quantity the dependencies are defined in bind.

        Args:
            beads: The beads object from which the bead positions are taken.
            nm: A normal modes object used to do the normal modes transformation.
            cell: The cell object from which the system box is taken.
            bforce: The forcefield object from which the force and virial are taken.
            prng: The random number generator object which controls random number
                generation.
        """

        # store local references to the different bits of the simulation
        self.beads = beads
        self.cell = cell
        self.forces = bforce
        self.prng = prng
        self.nm = nm
        self.ensemble = ens
        self.output_maker = omaker

    def step(self, step=None):
        """Dummy simulation time step which does nothing."""

        pass


dproperties(Motion, "dt")

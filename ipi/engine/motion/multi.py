"""TODO
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

from ipi.utils.depend import depend_value, dd
from ipi.engine.motion import Motion


class MultiMotion(Motion):

    """A class to hold multiple motion objects to be executed serially.
    """

    def __init__(self, motionlist=None):
        """Initialises MultiMotion.

        Args:
           fixcom: An optional boolean which decides whether the centre of mass
              motion will be constrained or not. Defaults to False.
        """

        dself = dd(self)
        dself.dt = depend_value(name="dt", func=self.get_totdt)
        self.mlist = motionlist
        for m in self.mlist:
            dd(m).dt.add_dependant(dself.dt)
        a = self.dt  # DON'T ASK WHY BUT IF YOU DON'T DO THAT WEAKREFS TO SELF.DT WILL BE INVALIDATED

        self.fixatoms = set(self.mlist[0].fixatoms)
        for m in self.mlist:
            self.fixatoms = self.fixatoms.intersection(m.fixatoms)
        self.fixatoms = list(self.fixatoms)

        self.fixcom = True  # fixcom is true only if all movers are fixed
        for m in self.mlist:
            self.fixcom = self.fixcom and m.fixcom

    def get_totdt(self):
        dt = 0.0
        for m in self.mlist:
            dt += m.dt
        return dt

    def step(self, step=None):
        for m in self.mlist:
            m.step(step)

    def bind(self, ens, beads, nm, cell, bforce, prng, omaker):
        """Binds beads, cell, bforce, and prng to the calculator.

        This takes a beads object, a cell object, a forcefield object and a
        random number generator object and makes them members of the atom motion caclulator.
        It also then creates the objects that will hold the data needed in the
        ensemble algorithms and the dependency network. Note that the conserved
        quantity is defined in the init, but as each ensemble has a different
        conserved quantity the dependencies are defined in bind.

        Args:
            beads: The beads object from whcih the bead positions are taken.
            nm: A normal modes object used to do the normal modes transformation.
            cell: The cell object from which the system box is taken.
            bforce: The forcefield object from which the force and virial are taken.
            prng: The random number generator object which controls random number
                generation.
        """

        for m in self.mlist:
            m.bind(ens, beads, nm, cell, bforce, prng, omaker)

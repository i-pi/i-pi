"""Implementation of a Smotion that can hold a list of Smotion objects that
will be executed serially at each step.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

from ipi.engine.smotion import Smotion, MetaDyn
from ipi.utils.messages import verbosity, warning

__all__ = ["MultiSmotion"]


class MultiSmotion(Smotion):
    """A class to hold multiple Smotion objects to be executed serially."""

    def __init__(self, smotionlist=None):
        """Initialises MultiSmotion.

        Args:
           smotionlist: A list of Smotion objects
        """

        self.mlist = smotionlist
        self.mode = "multi"

    def step(self, step=None):
        for m in self.mlist:
            m.step(step)

    def bind(self, syslist, prng, omaker):
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

        super(MultiSmotion, self).bind(syslist, prng, omaker)
        for k, m in enumerate(self.mlist):
            m.bind(syslist, prng, omaker)
            if (
                type(m) is MetaDyn
                and k != len(self.mlist) - 1
                and type(self.mlist[k + 1]) != MetaDyn
            ):
                warning(
                    "MetaD Smotion should come last in a multi Smotion to avoid a discrepancy between i-PI and PLUMED outputs",
                    verbosity.low,
                )

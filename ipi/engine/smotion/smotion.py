"""TODO"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


class Smotion:
    """Base smootion calculation class.

    Gives the standard methods and attributes needed in all the
    smootion calculation classes.

    Attributes:
        dummy: A dummy object giving dummy information.

    Depend objects:
        none
    """

    def __init__(self):
        """Initialises Smotion object.

        Args:

        """

        pass

    def bind(self, syslist, prng, omaker):
        """Binds systems and prng to the calculator.

        This takes a beads object, a cell object, a forcefield object and a
        random number generator object and makes them members of the atom motion caclulator.
        It also then creates the objects that will hold the data needed in the
        ensemble algorithms and the dependency network. Note that the conserved
        quantity is defined in the init, but as each ensemble has a different
        conserved quantity the dependencies are defined in bind.

        Args:
            syslist: A list of systems to operate on.
            prng: The random number generator object which controls random number
                generation.
        """

        # store local references to the different bits of the simulation
        self.syslist = syslist
        self.prng = prng
        self.output_maker = omaker

    def step(self, step=None):
        """Dummy simulation time step which does nothing."""

        pass

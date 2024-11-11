"""Contains the classes that deal with the MC exchanges of isotopes.

Holds the algorithms required for alchemical exchanges. Also calculates the
appropriate conserved energy quantity.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np

from ipi.engine.motion import Motion
from ipi.utils.depend import *
from ipi.utils.units import Constants


class AtomSwap(Motion):
    """Swap atom positions (typically useful to exchange species in
    a way that is compatible with the i-PI encapsulation paradigm.

    Attributes:
        names of the species for exchanges
    """

    def __init__(
        self, fixcom=False, fixatoms_dof=None, mode=None, names=[], nxc=1, ealc=None
    ):
        """Initialises a "alchemical exchange" motion object.

        Args:
            names : A list of isotopes
            nmc : frequency of doing exchanges

        """

        super(AtomSwap, self).__init__(fixcom=fixcom, fixatoms_dof=fixatoms_dof)

        self.names = names
        self.nxc = nxc

        self._ealc = depend_value(name="ealc")
        if ealc is not None:
            self.ealc = ealc
        else:
            self.ealc = 0.0

    def bind(self, ens, beads, cell, bforce, nm, prng, omaker):
        """Binds ensemble beads, cell, bforce, and prng to the dynamics.

        This takes a beads object, a cell object, a forcefield object and a
        random number generator object and makes them members of the ensemble.
        It also then creates the objects that will hold the data needed in the
        ensemble algorithms and the dependency network. Note that the conserved
        quantity is defined in the init, but as each ensemble has a different
        conserved quantity the dependencies are defined in bind.

        Args:
            beads: The beads object from whcih the bead positions are taken.
            prng: The random number generator object which controls random number
                generation.
        """

        super(AtomSwap, self).bind(ens, beads, cell, bforce, nm, prng, omaker)
        self.ensemble.add_econs(self._ealc)
        self.dbeads = self.beads.clone()
        self.dcell = self.cell.clone()
        self.dforces = self.forces.clone(self.dbeads, self.dcell)

    def AXlist(self, atomtype):
        """This compile a list of atoms ready for exchanges."""

        # selects the types of atoms for exchange
        atomexchangelist = []
        for i in range(self.beads.natoms):
            if self.beads.names[i] in atomtype:
                atomexchangelist.append(i)

        return np.asarray(atomexchangelist)

    def step(self, step=None):
        # picks number of attempted exchanges
        ntries = self.prng.poisson(self.nxc)
        if ntries == 0:
            return

        """Does one round of alchemical exchanges."""
        # record the spring energy (divided by mass) for each atom in the exchange chain
        nb = self.beads.nbeads
        # qswap = np.zeros((nb, 3))
        axlist = self.AXlist(self.names)
        lenlist = len(axlist)
        if lenlist == 0:
            raise ValueError("Atoms exchange list is empty in MC atom swapper.")

        # does the exchange
        betaP = 1.0 / (Constants.kb * self.ensemble.temp * nb)
        nexch = 0

        # this would be double-counting, we already have a bail-out condition above
        # if (1.0/self.nxc < self.prng.u) : return  # tries a round of exhanges with probability 1/nmc
        self.dcell.h = (
            self.cell.h
        )  # just in case the cell gets updated in the other motion classes
        for x in range(ntries):
            i = self.prng.integers(0, lenlist)
            j = self.prng.integers(0, lenlist)
            while self.beads.names[axlist[i]] == self.beads.names[axlist[j]]:
                j = self.prng.integers(0, lenlist)  # makes sure we pick a real exchange

            # map the "subset" indices back to the "absolute" atom indices
            i = axlist[i]
            j = axlist[j]

            old_energy = self.forces.pot
            # swap the atom positions
            self.dbeads.q[:] = self.beads.q[:]
            self.dbeads.q[:, 3 * i : 3 * i + 3] = self.beads.q[:, 3 * j : 3 * j + 3]
            self.dbeads.q[:, 3 * j : 3 * j + 3] = self.beads.q[:, 3 * i : 3 * i + 3]
            new_energy = self.dforces.pot
            pexchange = np.exp(-betaP * (new_energy - old_energy))

            # attemps the exchange, and actually propagate the exchange if something has happened
            if pexchange > self.prng.u:
                nexch += 1

                # copy the exchanged beads position
                self.beads.q[:] = self.dbeads.q[:]
                # transfers the (already computed) status of the force, so we don't need to recompute
                self.forces.transfer_forces(self.dforces)

                self.ealc += -(new_energy - old_energy)


dproperties(AtomSwap, ["ealc"])

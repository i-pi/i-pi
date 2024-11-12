"""Contains the classes that deal with the MC exchanges of isotopes.

Holds the algorithms required for alchemical exchanges. Also calculates the
appropriate conserved energy quantity.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np

from ipi.engine.motion import Motion
from ipi.utils.depend import dproperties, depend_value, dstrip
from ipi.utils.units import Constants


# __all__ = ['AlchemyMC']


class AlchemyMC(Motion):
    """Monte Carlo alchemical exchanges class.

    Attributes:
        names of the isotopes for exchanges

    Depend objects:
        The spring potential energy, names of the atoms.
    """

    def __init__(
        self, fixcom=False, fixatoms_dof=None, mode=None, names=[], nxc=1, ealc=None
    ):
        """Initialises a "alchemical exchange" motion object.

        Args:
            names : A list of isotopes
            nmc : frequency of doing exchanges

        """

        super(AlchemyMC, self).__init__(fixcom=fixcom, fixatoms_dof=fixatoms_dof)

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

        super(AlchemyMC, self).bind(ens, beads, cell, bforce, nm, prng, omaker)
        self.ensemble.add_econs(self._ealc)

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
        # q = dstrip(self.beads.q)
        nb = self.beads.nbeads
        axlist = self.AXlist(self.names)
        lenlist = len(axlist)
        if lenlist == 0:
            raise ValueError("Atoms exchange list is empty in alchemical sampler.")
        atomspring = np.zeros(lenlist)
        wk2 = dstrip(self.nm.omegak2)

        i = 0
        for atomnum in axlist:
            na3 = atomnum * 3

            # computes spring in NM representation
            spr = 0.0
            for b in range(1, nb):
                spr += wk2[b] * (self.nm.qnm[b, na3 : na3 + 3] ** 2).sum()
            spr *= 0.5

            # no mass here - just the massless spring term
            atomspring[i] = spr
            i += 1

        # does the exchange
        betaP = 1.0 / (Constants.kb * self.ensemble.temp * nb)
        nexch = 0

        # this would be double-counting, we already have a bail-out condition above
        # if (1.0/self.nxc < self.prng.u) : return  # tries a round of exhanges with probability 1/nmc

        for x in range(ntries):
            i = self.prng.integers(0, lenlist)
            j = self.prng.integers(0, lenlist)
            while self.beads.names[axlist[i]] == self.beads.names[axlist[j]]:
                j = self.prng.integers(0, lenlist)  # makes sure we pick a real exchange

            # energy change due to the swap
            difspring = (atomspring[i] - atomspring[j]) * (
                self.beads.m[axlist[j]] - self.beads.m[axlist[i]]
            )
            pexchange = np.exp(-betaP * difspring)

            # attemps the exchange
            if pexchange > self.prng.u:
                nexch += 1
                # print 'exchange atom No.  ', axlist[i], '  and  ', axlist[j]

                # swap names
                nameswap = self.beads.names[axlist[i]]
                self.beads.names[axlist[i]] = self.beads.names[axlist[j]]
                self.beads.names[axlist[j]] = nameswap

                # change masses
                massratio = self.beads.m[axlist[i]] / self.beads.m[axlist[j]]
                self.beads.m[axlist[i]] /= massratio
                self.beads.m[axlist[j]] *= massratio

                # adjust the (classical) momenta to conserve ke
                self.beads.p[:, 3 * axlist[i] : 3 * (axlist[i] + 1)] /= np.sqrt(
                    massratio
                )
                self.beads.p[:, 3 * axlist[j] : 3 * (axlist[j] + 1)] *= np.sqrt(
                    massratio
                )

                # adjusts the conserved quantity counter based on the change in spring energy
                self.ealc += -difspring


dproperties(AlchemyMC, ["ealc"])

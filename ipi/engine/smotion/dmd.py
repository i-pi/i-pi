"""Holds an algorithm to update an external time dependent drive.
   Should work with PBC.

"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2016 i-PI developers
# See the "licenses" directory for full license information.


from ipi.engine.smotion import Smotion


__all__ = ["DMD"]


class DMD(Smotion):
    """Class that updates the time step for a time-dependent forcefield.
    Solves issues when a i-PI native time-dependent forcefield needs to be
    called several times (several beads or replicas) at the same time step
    and only when all evaluations have been performed the "time" of the forcecield
    should be updated. Depending on how the actual forcefield is implemented, the
    update step can be made more intricate than just a time increment.


    Attributes:

    """

    def __init__(self, dmdff=""):
        """Initialises DMD smotion.

        Args:
        """

        super(DMD, self).__init__()
        self.dmdff = dmdff

    def step(self, step=None):
        """Updates driven md time step."""

        for s in self.syslist:
            for ik, bc in enumerate(s.ensemble.bcomp):
                k = bc.ffield
                if k not in self.dmdff:
                    continue  # only does dmd for the indicated forcefield
                f = s.ensemble.bias.ff[k]  # dmd is sort of a bias...
                if not hasattr(
                    f, "dmd_update"
                ):  # forcefield does not expose dmd_update interface
                    raise ValueError(
                        "The forcefield '%s' associated with driven MD \
                                      does not have a dmd_update interface"
                        % (k)
                    )
                f.dmd_update()  # updates the time step

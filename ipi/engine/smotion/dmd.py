"""Holds the algorithms to perform replica exchange.

Algorithms implemented by Robert Meissner and Riccardo Petraglia, 2016
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2016 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np
import time

from ipi.engine.smotion import Smotion
from ipi.engine.ensembles import ensemble_swap
from ipi.utils.depend import *
from ipi.utils.softexit import softexit
from ipi.utils.messages import verbosity, info


__all__ = ["DMD"]


class DMD(Smotion):
    """Metadynamics routine based on a FFPlumed forcefield.

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
            #            oldf = dstrip(s.forces.f).copy()
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
                #                if s.ensemble.bweights[ik] == 0:
                #                    continue  # do not put metad bias on biases with zero weights (useful to do remd+metad!)
                f.dmd_update()  # updates the time step


#                if fmtd:  # if metadyn has updated, then we must recompute forces.
#                    # hacky but cannot think of a better way: we must manually taint *just* that component
#                    for fc in s.ensemble.bias.mforces:
#                        if fc.ffield == k:
#                            for fb in fc._forces:
#                                dd(fb).ufvx.taint()

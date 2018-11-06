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


__all__ = ['MetaDyn']


class MetaDyn(Smotion):
    """Metadynamics routine based on a FFPlumed forcefield.

    Attributes:

    """

    def __init__(self, metaff=""):
        """Initialises REMD.

        Args:
        """

        super(MetaDyn, self).__init__()
        self.metaff = metaff

    def bind(self, syslist, prng):

        super(MetaDyn, self).bind(syslist, prng)

    def step(self, step=None):
        """Updates metad bias."""

        for s in self.syslist:
            oldf = dstrip(s.forces.f).copy()
            for ik, bc in enumerate(s.ensemble.bcomp):
                k = bc.ffield
                if not k in self.metaff:
                    continue  # only does metad for the indicated forcefield
                f = s.ensemble.bias.ff[k]
                if not hasattr(f, "mtd_update"):  # forcefield does not expose mtd_update interface
                    raise ValueError("The forcefield '%s' associated with metadynamics \
                                      does not have a mtd_update interface" % (k))
                if s.ensemble.bweights[ik] == 0:
                    continue # do not put metad bias on biases with zero weights (useful to do remd+metad!)
                fmtd = f.mtd_update(pos=s.beads.qc, cell=s.cell.h)
                if fmtd:  # if metadyn has updated, then we must recompute forces.
                    # hacky but cannot think of a better way: we must manually taint *just* that component
                    for fc in s.ensemble.bias.mforces:
                        if fc.ffield == k:
                            for fb in fc._forces:
                                dd(fb).ufvx.taint()

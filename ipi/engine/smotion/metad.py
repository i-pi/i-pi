"""Holds the algorithms to perform replica exchange.

Algorithms implemented by Robert Meissner and Riccardo Petraglia, 2016
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2016 i-PI developers
# See the "licenses" directory for full license information.

from ipi.engine.smotion import Smotion

__all__ = ["MetaDyn"]


class MetaDyn(Smotion):
    """Metadynamics routine based on a FFPlumed forcefield.

    Attributes:

    """

    def __init__(self, metaff="", use_energy=False):
        """Initialises MetaDynamics.

        Args:
        """

        super(MetaDyn, self).__init__()
        self.metaff = metaff
        self.use_energy = use_energy

        self.mode = "metad"

    def bind(self, syslist, prng, omaker):
        super(MetaDyn, self).bind(syslist, prng, omaker)

        for s in self.syslist:
            for ik, bc in enumerate(s.ensemble.bcomp):
                k = bc.ffield
                if k not in self.metaff:
                    continue  # only does metad for the indicated forcefield
                f = s.ensemble.bias.ff[k]
                if not hasattr(
                    f, "mtd_update"
                ):  # forcefield does not expose mtd_update interface
                    raise ValueError("The forcefield '%s' associated with metadynamics \
                                      does not have a mtd_update interface" % (k))

                if self.use_energy and s.ensemble.bweights[ik] > 0:
                    if s.beads.nbeads > 1:
                        raise ValueError(
                            "Cannot use energy information for metadynamics with more than one bead"
                        )
                    # injects a reference to the force calculator in the plumed forcefield.
                    # this is a really ugly hack because it breaks the separation between forces and
                    # forcefield calculators, making the whole setup very fragile.
                    # for instance, this will *not* function with multiple beads (because there is no way to
                    # specify which beads' energy should be taken) or with replica exchange (because
                    # metad_update is not called when the ensembles are exchanged)
                    f.system_force = s.forces

    def step(self, step=None):
        """Updates metad bias."""

        for s in self.syslist:
            for ik, bc in enumerate(s.ensemble.bcomp):
                k = bc.ffield
                if k not in self.metaff:
                    continue  # only does metad for the indicated forcefield

                f = s.ensemble.bias.ff[k]
                if not hasattr(
                    f, "mtd_update"
                ):  # forcefield does not expose mtd_update interface
                    raise ValueError("The forcefield '%s' associated with metadynamics \
                                      does not have a mtd_update interface" % (k))
                if s.ensemble.bweights[ik] == 0:
                    continue  # do not put metad bias on biases with zero weights (useful to do remd+metad!)

                # MTD is hardcoded to be applied on the centroid variable.
                # this is the "right" thing if you need to compute kinetics based
                # on the resulting FES
                mtd_work = f.mtd_update(pos=s.beads.qc, cell=s.cell.h)

                # updates the conserved quantity with the change in bias so that
                # we remove the shift due to added hills
                s.ensemble.eens += (
                    mtd_work * s.beads.nbeads
                )  # apply ring polymer contraction!

                # Adaptive biases can change forces even when their scalar work
                # is zero. Invalidate only this component; recompute lazily on
                # the next read and leave the physical force caches untouched.
                for fc in s.ensemble.bias.mforces:
                    if fc.ffield == k:
                        for fb in fc._forces:
                            fb._ufvx.taint()

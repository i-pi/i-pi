"""Deals with creating the ensembles class.

Copyright (C) 2013, Joshua More and Michele Ceriotti

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http.//www.gnu.org/licenses/>.


Classes:
   InputAlKMC: Deals with creating the a KMC object from a file, and
      writing the checkpoints.
"""

import numpy as np
import pickle

from ipi.engine.motion import *
from ipi.utils.inputvalue import *
from ipi.inputs.thermostats import *
from ipi.inputs.initializer import *
from ipi.utils.units import *
from .geop import InputGeop

__all__ = ["InputAlKMC"]


class InputAlKMC(InputDictionary):
    """KMC for Al-6xxx.

    Contains options related with geometry optimization, such as method,
    thresholds, linear search strategy, etc.

    """

    attribs = {
        "mode": (
            InputAttribute,
            {
                "dtype": str,
                "default": "rfkmc",
                "help": "The KMC algorithm to be used",
                "options": ["rfkmc"],
            },
        )
    }

    # options of the method nstep, a0, ncell, nvac, nsi, nmg, state="",
    fields = {
        "optimizer": (
            InputGeop,
            {"default": {}, "help": "Option for geometry optimization step"},
        ),
        "nstep": (
            InputValue,
            {"dtype": int, "default": 10, "help": "The number of optimization steps."},
        ),
        "a0": (
            InputValue,
            {
                "dtype": float,
                "dimension": "length",
                "default": 1.0,
                "help": "FCC lattice parameter ",
            },
        ),
        "diffusion_barrier_al": (
            InputValue,
            {
                "dtype": float,
                "dimension": "energy",
                "default": 0.01,
                "help": "Barrier for vacancy diffusion in pure Al.",
            },
        ),
        "diffusion_prefactor_al": (
            InputValue,
            {
                "dtype": float,
                "dimension": "frequency",
                "default": 2.4188843e-05,
                "help": "Prefactor for vacancy diffusion in pure Al.",
            },
        ),
        "diffusion_barrier_mg": (
            InputValue,
            {
                "dtype": float,
                "dimension": "energy",
                "default": 0,
                "help": "Barrier for vacancy-assisted diffusion of Mg.",
            },
        ),
        "diffusion_prefactor_mg": (
            InputValue,
            {
                "dtype": float,
                "dimension": "frequency",
                "default": 0,
                "help": "Prefactor for vacancy-assisted diffusion of Mg.",
            },
        ),
        "diffusion_barrier_si": (
            InputValue,
            {
                "dtype": float,
                "dimension": "energy",
                "default": 0,
                "help": "Barrier for vacancy-assisted diffusion of Si.",
            },
        ),
        "diffusion_prefactor_si": (
            InputValue,
            {
                "dtype": float,
                "dimension": "frequency",
                "default": 0,
                "help": "Prefactor for vacancy-assisted diffusion of Si.",
            },
        ),
        "neval": (
            InputValue,
            {
                "dtype": int,
                "default": 4,
                "help": "The number of parallel force evaluators.",
            },
        ),
        "ncell": (
            InputValue,
            {
                "dtype": int,
                "default": 4,
                "help": "The number of repeat cells in each direction.",
            },
        ),
        "nvac": (
            InputValue,
            {"dtype": int, "default": 4, "help": "The number of vacancies."},
        ),
        "nsi": (
            InputValue,
            {"dtype": int, "default": 4, "help": "The number of silicon atoms."},
        ),
        "nmg": (
            InputValue,
            {"dtype": int, "default": 4, "help": "The number of magnesium atoms."},
        ),
        "idx": (
            InputArray,
            {
                "dtype": int,
                "default": input_default(factory=np.zeros, args=(0, int)),
                "help": "The position of the atoms on the lattice, relative to the canonical ordering.",
            },
        ),
        "tottime": (
            InputValue,
            {
                "dtype": float,
                "dimension": "time",
                "default": 0.0,
                "help": "Total KMC time elapsed ",
            },
        ),
        "ecache_file": (
            InputValue,
            {
                "dtype": str,
                "default": "",
                "help": "Filename for storing/loading energy cache",
            },
        ),
        "qcache_file": (
            InputValue,
            {
                "dtype": str,
                "default": "",
                "help": "Filename for storing/loading positions cache",
            },
        ),
        "max_cache_len": (
            InputValue,
            {
                "dtype": int,
                "default": 1000,
                "help": "Maximum cache length before oldest entry is deleted",
            },
        ),
    }

    STORE_STRIDE = 1.1

    dynamic = {}

    default_help = "Holds all the information for the KMC dynamics, such as timestep, rates and barriers that control it."
    default_label = "AL6XXXKMC"

    def store(self, kmc):
        """Takes a kinetic MonteCarlo instance and stores a minimal representation of it.

        Args:
            kmc: A kMC object.
        """

        if kmc == {}:
            return

        # """self.state.store(kmc.state)
        # self.cell.h.store(kmc.cell.h)
        # self.beads.q.store(kmc.beads.q)
        # self.dt(kmc.dt)"""

        # self.mode.store(kmc.mode)
        self.a0.store(kmc.a0)
        self.nvac.store(kmc.nvac)
        self.nmg.store(kmc.nmg)
        self.nsi.store(kmc.nsi)
        self.ncell.store(kmc.ncell)
        self.nstep.store(kmc.nstep)
        self.diffusion_barrier_al.store(kmc.diffusion_barrier_al)
        self.diffusion_prefactor_al.store(kmc.diffusion_prefactor_al)
        self.diffusion_barrier_mg.store(kmc.diffusion_barrier_mg)
        self.diffusion_prefactor_mg.store(kmc.diffusion_prefactor_mg)
        self.diffusion_barrier_si.store(kmc.diffusion_barrier_si)
        self.diffusion_prefactor_si.store(kmc.diffusion_prefactor_si)
        self.idx.store(kmc.idx)
        self.tottime.store(kmc.tottime)
        self.ecache_file.store(kmc.ecache_file)
        self.qcache_file.store(kmc.qcache_file)
        self.max_cache_len.store(kmc.max_cache_len)

        # only stores cache after a decent amount of new structures have been found
        if (
            kmc.struct_count - kmc.ncache_stored
        ) >= 0.1 * kmc.ncache:  # dump if new structures exceed 10% of current store
            # if (kmc.struct_count - kmc.ncache_stored) >= 100 : # Basically dump only after 100 new structures.
            if kmc.ecache_file != "":
                print(
                    "10%% new structures since last dump. Storing ECACHE in ",
                    kmc.ecache_file,
                )
                ff = open(kmc.ecache_file, "wb")
                pickle.dump(kmc.ecache, ff)
                ff.close()
            if kmc.qcache_file != "":
                print(
                    "10%% new structures since last dump. Storing QCACHE in ",
                    kmc.qcache_file,
                )
                ff = open(kmc.qcache_file, "wb")
                pickle.dump(kmc.qcache, ff)
                ff.close()
            # kmc.ncache_stored = kmc.ncache
            kmc.ncache_stored = kmc.struct_count

    def fetch(self):
        rv = super(InputAlKMC, self).fetch()
        rv["mode"] = self.mode.fetch()
        return rv

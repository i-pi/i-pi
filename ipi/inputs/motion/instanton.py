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
   InputInst: Deals with creating the Inst object from a file, and
      writing the checkpoints.
"""

import numpy as np
from ipi.utils.inputvalue import *


__all__ = ["InputInst"]


class InputInst(InputDictionary):
    """Instanton optimization options.

    Contains options related with instanton, such as method,
    thresholds, hessian update strategy, etc.

    """

    attribs = {
        "mode": (
            InputAttribute,
            {
                "dtype": str,
                "default": "rate",
                "help": "Defines whether it is an instanton rate or instanton tunneling splitting calculaion",
                "options": ["rate", "splitting"],
            },
        )
    }

    fields = {
        "tolerances": (
            InputDictionary,
            {
                "dtype": float,
                "options": ["energy", "force", "position"],
                "default": [1e-5, 1e-4, 1e-3],
                "help": "Convergence criteria for optimization.",
                "dimension": ["energy", "force", "length"],
            },
        ),
        "biggest_step": (
            InputValue,
            {
                "dtype": float,
                "default": 0.4,
                "help": "The maximum step size during the optimization.",
            },
        ),
        "old_pos": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "The previous step positions during the optimization. ",
                "dimension": "length",
            },
        ),
        "old_pot": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "The previous step potential energy during the optimization",
                "dimension": "energy",
            },
        ),
        "old_force": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "The previous step force during the optimization",
                "dimension": "force",
            },
        ),
        "opt": (
            InputValue,
            {
                "dtype": str,
                "default": "None",
                "options": ["nichols", "NR", "lbfgs", "lanczos", "None"],
                "help": """The geometry optimization algorithm to be used.
                                            For small system sizes nichols is recomended. Lanczos is tailored for big bigger than nbeads*natoms >~38*64.
                                            NR works in both cases given that the initial guess is close to the optimized geometry.
                                            Finally lbfgs is used for tunneling splitting calculations. """,
            },
        ),
        "max_e": (
            InputValue,
            {
                "dtype": float,
                "default": 0.0,
                "help": """Evaluate the forces in a reduced ring polymer such that the potential energy between consecutive replicas is smaller that the provided value.""",
                "dimension": "energy",
            },
        ),
        "max_ms": (
            InputValue,
            {
                "dtype": float,
                "default": 0.0,
                "help": """Evaluate the forces in a reduced ring polymer such that that mass-scaled distance in a.u. between consecutive replicas is  smaller that the provided value.""",
            },
        ),
        "discretization": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.ones, args=(0,)),
                "help": "Allows to specified non uniform time discretization as proposed in J. Chem. Phys. 134, 184107 (2011)",
            },
        ),
        "friction": (
            InputValue,
            {
                "dtype": bool,
                "default": False,
                "help": "Activates Friction. Add additional terms to the RP related to a position-independent frictional force. See Eq. 20 in J. Chem. Phys. 156, 194106 (2022)",
            },
        ),
        "frictionSD": (
            InputValue,
            {
                "dtype": bool,
                "default": True,
                "help": "Activates SD Friction. Add additional terms to the RP related to a position-dependent frictional force. See Eq. 32 in J. Chem. Phys. 156, 194106 (2022)",
            },
        ),
        "fric_spec_dens": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.ones, args=(0,)),
                "help": "Laplace Transform (LT) of friction. A two column data is expected. First column: w (cm^-1). Second column: LT(eta)(w). See Eq. 11 in J. Chem. Phys. 156, 194106 (2022). Note that within the separable coupling approximation the frequency dependence of the friction tensor is position independent.",
            },
        ),
        "fric_spec_dens_ener": (
            InputValue,
            {
                "dtype": float,
                "default": 0.0,
                "help": "Energy at which the LT of the friction tensor is evaluated in the client code",
                "dimension": "energy",
            },
        ),
        "eta": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.eye, args=(0,)),
                "help": "Friction Tensor. Only to be used when frictionSD is disabled.",
            },
        ),
        "alt_out": (
            InputValue,
            {
                "dtype": int,
                "default": 1,
                "help": """Alternative output:Prints different formatting of outputs for geometry, hessian and bead potential energies.
                                               All quantities are also accessible from typical i-pi output infrastructure.
                                               Default to 1, which prints every step. -1 will suppress the output (except the last one).
                                               Any other positive number will set the frequency (in steps) with which the quantities are
                                               written to file.
                                               The instanton geometry is printed in xyz format and the distances are in angrstroms
                                               The hessian is printed in one line with the following format:
                                               h1_1,h2_1,...,hN_1,   h2_2,h2_2,hN_2,   ....   ,h1_d,h2_d,...,hN_d.
                                               Where N represents the total number of replicas, d the number of dimension of each replica (3*n_atoms) and
                                               hi_j means the row j of the physical hessian corresponding to the replica i.
                                               The physical hessian uses a convention according to the positions convention used in  i-pi.
                                               Example of 2 particles, the first two rows of the physical hessian reads:
                                               'H_x1_x1, H_x1_y1, H_x1_z1, H_x1_x2, H_x1_y2,H_x1_z2'
                                               'H_x2_x1, H_x2_y1, H_x2_z1, H_x2_x2, H_x2_y2,H_x2_z2' """,
            },
        ),
        "prefix": (
            InputValue,
            {
                "dtype": str,
                "default": "instanton",
                "help": "Prefix of the output files.",
            },
        ),
        "delta": (
            InputValue,
            {"dtype": float, "default": 0.1, "help": "Initial stretch amplitude."},
        ),
        # Hessian
        "hessian_init": (
            InputValue,
            {
                "dtype": bool,
                "default": False,
                "help": "How to initialize the hessian if it is not fully provided.",
            },
        ),
        "hessian": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.eye, args=(0,)),
                "help": "(Approximate) Hessian.",
            },
        ),
        "hessian_update": (
            InputValue,
            {
                "dtype": str,
                "default": "powell",
                "options": ["powell", "recompute"],
                "help": "How to update the hessian after each step.",
            },
        ),
        "hessian_asr": (
            InputValue,
            {
                "dtype": str,
                "default": "none",
                "options": ["none", "poly", "crystal"],
                "help": "Removes the zero frequency vibrational modes depending on the symmetry of the system.",
            },
        ),
        "fric_hessian": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.eye, args=(0,)),
                "help": "(Approximate) friction second derivative from which a friction Hessian can be built.",
            },
        ),
        # L-BFGS
        "qlist_lbfgs": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "List of previous position differences for L-BFGS, if known.",
            },
        ),
        "glist_lbfgs": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "List of previous gradient differences for L-BFGS, if known.",
            },
        ),
        "old_direction": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "The previous direction in a CG or SD optimization.",
            },
        ),
        "scale_lbfgs": (
            InputValue,
            {
                "dtype": int,
                "default": 2,
                "help": """Scale choice for the initial hessian.
                                                       0 identity.
                                                       1 Use first member of position/gradient list.
                                                       2 Use last  member of position/gradient list.""",
            },
        ),
        "corrections_lbfgs": (
            InputValue,
            {
                "dtype": int,
                "default": 20,
                "help": "The number of past vectors to store for L-BFGS.",
            },
        ),
        "ls_options": (
            InputDictionary,
            {
                "dtype": [float, int],
                "help": """Options for line search methods. Includes:
                                  tolerance: stopping tolerance for the search,
                                  iter: the maximum number of iterations,
                                  step: initial step for bracketing,
                                  adaptive: whether to update initial step.
                                  """,
                "options": ["tolerance", "iter"],
                "default": [0.2, 100],
                "dimension": ["energy", "undefined"],
            },
        ),
        # Final calculations
        "energy_shift": (
            InputValue,
            {
                "dtype": float,
                "default": 0.000,
                "help": "Set the zero of energy.",
                "dimension": "energy",
            },
        ),
        "hessian_final": (
            InputValue,
            {
                "dtype": bool,
                "default": False,
                "help": "Decide if we are going to compute the final big-hessian by finite difference.",
            },
        ),
    }

    dynamic = {}

    default_help = "A class for instanton calculations"
    default_label = "INSTANTON"

    def store(self, geop):
        if geop == {}:
            return

        options = geop.options
        optarrays = geop.optarrays

        # Optimization mode
        self.mode.store(options["mode"])

        # Generic optimization
        self.tolerances.store(options["tolerances"])
        self.biggest_step.store(optarrays["big_step"])
        self.opt.store(options["opt"])

        # Generic instanton
        self.max_e.store(options["max_e"])
        self.max_ms.store(options["max_ms"])
        self.discretization.store(options["discretization"])
        self.friction.store(options["friction"])
        self.frictionSD.store(options["frictionSD"])
        if options["friction"]:
            self.fric_spec_dens.store(options["fric_spec_dens"])
            self.fric_spec_dens_ener.store(options["fric_spec_dens_ener"])
        self.alt_out.store(options["save"])
        self.prefix.store(options["prefix"])
        self.delta.store(optarrays["delta"])
        self.hessian_final.store(options["hessian_final"])
        self.old_pot.store(optarrays["old_u"])
        self.old_force.store(optarrays["old_f"])
        self.energy_shift.store(optarrays["energy_shift"])

        # Now we decide what to store depending on optimization algorithm
        if (
            geop.options["opt"] == "nichols"
            or geop.options["opt"] == "NR"
            or geop.options["opt"] == "lanczos"
        ):
            self.hessian.store(optarrays["hessian"])
            self.hessian_update.store(options["hessian_update"])
            self.hessian_asr.store(options["hessian_asr"])
            if options["friction"]:
                if options["frictionSD"]:
                    self.fric_hessian.store(optarrays["fric_hessian"])
                else:
                    self.eta.store(options["eta0"])
        elif geop.options["opt"] == "lbfgs":
            self.qlist_lbfgs.store(optarrays["qlist"])
            self.glist_lbfgs.store(optarrays["glist"])
            self.old_direction.store(optarrays["d"])
            self.scale_lbfgs.store(options["scale"])
            self.corrections_lbfgs.store(options["corrections"])
            self.ls_options.store(options["ls_options"])
            self.hessian_final.store(options["hessian_final"])
            if options["hessian_final"] == "true":
                self.hessian.store(optarrays["hessian"])

    def fetch(self):
        rv = super(InputInst, self).fetch()
        rv["mode"] = self.mode.fetch()
        return rv

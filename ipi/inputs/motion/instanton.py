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
import sys


__all__ = ['InputInst']


class InputInst(InputDictionary):

    """Instanton optimization options.

    Contains options related with instanton, such as method,
    thresholds, hessian update strategy, etc.

    """
    attribs = {"mode": (InputAttribute, {"dtype": str, "default": "rate",
                                         "help": "Defines whether it is an instanton rate or instanton tunneling splitting calculaion",
                                         "options": ['rate', 'splitting']})}

    fields = {"tolerances": (InputDictionary, {"dtype": float,
                                               "options": ["energy", "force", "position"],
                                               "default": [1e-5, 1e-4, 1e-3],
                                               "help": "Convergence criteria for optimization.",
                                               "dimension": ["energy", "force", "length"]}),
              "biggest_step": (InputValue, {"dtype": float,
                                            "default": 0.4,
                                            "help": "The maximum step size during the optimization."}),
              "old_pos": (InputArray, {"dtype": float,
                                       "default": input_default(factory=np.zeros, args=(0,)),
                                       "help": "The previous step positions during the optimization. ",
                                       "dimension": "length"}),
              "old_pot": (InputArray, {"dtype": float,
                                       "default": input_default(factory=np.zeros, args=(0,)),
                                       "help": "The previous step potential energy during the optimization",
                                       "dimension": "energy"}),
              "old_force": (InputArray, {"dtype": float,
                                         "default": input_default(factory=np.zeros, args=(0,)),
                                         "help": "The previous step force during the optimization",
                                         "dimension": "force"}),
              "opt": (InputValue, {"dtype": str,
                                   "default": 'None',
                                   "options": ["nichols", "NR", "lbfgs", "None"],
                                   "help": "The geometry optimization algorithm to be used"}),
              "alt_out": (InputValue, {"dtype": int,
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
                                               'H_x2_x1, H_x2_y1, H_x2_z1, H_x2_x2, H_x2_y2,H_x2_z2' """}),
              "prefix": (InputValue, {"dtype": str,
                                      "default": "instanton",
                                      "help": "Prefix of the output files."}),
              "delta": (InputValue, {"dtype": float,
                                     "default": 0.1,
                                     "help": "Initial stretch amplitude."}),
              # Hessian
              "hessian_init": (InputValue, {"dtype": str,
                                            "default": 'false',
                                            "options": ["true", 'false'],
                                            "help": "How to initialize the hessian if it is not fully provided."}),
              "hessian": (InputArray, {"dtype": float,
                                       "default": input_default(factory=np.eye, args=(0,)),
                                       "help": "(Approximate) Hessian."}),
              "hessian_update": (InputValue, {"dtype": str,
                                              "default": "powell",
                                              "options": ["powell", "recompute"],
                                              "help": "How to update the hessian after each step."}),
              "hessian_asr": (InputValue, {"dtype": str,
                                           "default": "none",
                                           "options": ["none", "poly", "crystal"],
                                           "help": "Removes the zero frequency vibrational modes depending on the symmerty of the system."}),
              # L-BFGS
              "qlist_lbfgs": (InputArray, {"dtype": float,
                                           "default": input_default(factory=np.zeros, args=(0,)),
                                           "help": "List of previous position differences for L-BFGS, if known."}),
              "glist_lbfgs": (InputArray, {"dtype": float,
                                           "default": input_default(factory=np.zeros, args=(0,)),
                                           "help": "List of previous gradient differences for L-BFGS, if known."}),
              "old_direction": (InputArray, {"dtype": float,
                                             "default": input_default(factory=np.zeros, args=(0,)),
                                             "help": "The previous direction in a CG or SD optimization."}),
              "scale_lbfgs": (InputValue, {"dtype": int,
                                           "default": 2,
                                           "help": """Scale choice for the initial hessian.
                                                       0 identity.
                                                       1 Use first member of position/gradient list.
                                                       2 Use last  member of position/gradient list."""}),
              "corrections_lbfgs": (InputValue, {"dtype": int,
                                                 "default": 20,
                                                 "help": "The number of past vectors to store for L-BFGS."}),
              "ls_options": (InputDictionary, {"dtype": [float, int],
                                               "help": """"Options for line search methods. Includes:
                                  tolerance: stopping tolerance for the search,
                                  iter: the maximum number of iterations,
                                  step: initial step for bracketing,
                                  adaptive: whether to update initial step.
                                  """,
                                               "options": ["tolerance", "iter"],
                                               "default": [0.2, 100],
                                               "dimension": ["energy", "undefined"]}),
              # Final calculations
              "energy_shift": (InputValue, {"dtype": float, "default": 0.000,
                                            "help": "Set the zero of energy.",
                                            "dimension": "energy"
                                            }),
              "hessian_final": (InputValue, {"dtype": str,
                                             "default": "false",
                                             "options": ["false", "true"],
                                             "help": "Decide if we are going to compute the final big-hessian by finite difference."})
              }

    dynamic = {}

    default_help = "A class for instanton calculations"
    default_label = "instanton"

    def store(self, geop):
        if geop == {}:
            return

        # Optimization mode
        self.mode.store(geop.mode)

        # Generic optimization
        self.tolerances.store(geop.tolerances)
        self.biggest_step.store(geop.big_step)
        self.opt.store(geop.opt)

        # Generic instanton
        self.alt_out.store(geop.save)
        self.prefix.store(geop.prefix)
        self.delta.store(geop.delta)
        self.hessian_final.store(geop.hessian_final)
        self.old_pot.store(geop.old_u)
        self.old_force.store(geop.old_f)
        self.energy_shift.store(geop.energy_shift)

        # Now we decide what to store depending on optimization algorithm
        if geop.opt == 'nichols' or geop.opt == 'NR':
            self.hessian.store(geop.hessian)
            self.hessian_update.store(geop.hessian_update)
            self.hessian_asr.store(geop.hessian_asr)
            self.hessian_final.store(geop.hessian_final)
        elif geop.opt == 'lbfgs':
            self.qlist_lbfgs.store(geop.qlist)
            self.glist_lbfgs.store(geop.glist)
            self.old_direction.store(geop.d)
            self.scale_lbfgs.store(geop.scale)
            self.corrections_lbfgs.store(geop.corrections)
            self.ls_options.store(geop.ls_options)
            self.hessian_final.store(geop.hessian_final)
            if geop.hessian_final == 'true':
                self.hessian.store(geop.hessian)

    def fetch(self):
        rv = super(InputInst, self).fetch()
        rv["mode"] = self.mode.fetch()
        return rv

"""Creates objects that deal with the different ensembles."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

import numpy as np
from copy import copy
import ipi.engine.thermostats
import ipi.engine.barostats
from ipi.utils.constrtools import (
    RigidBondConstraint,
    AngleConstraint,
    EckartConstraint,
    ConstraintList,
)
from ipi.utils.inputvalue import (
    InputDictionary,
    InputAttribute,
    InputValue,
    InputArray,
    Input,
    input_default,
)
from ipi.inputs.barostats import InputBaro
from ipi.inputs.thermostats import InputThermo

__all__ = ["InputConstrainedDynamics", "InputConstraint", "InputConstraintSolver"]


class InputConstraintSolver(InputDictionary):
    fields = {
        "tolerance": (
            InputValue,
            {
                "dtype": float,
                "default": 0.0001,
                "help": "Tolerance value used in the Quasi-Newton iteration scheme.",
            },
        ),
        "maxit": (
            InputValue,
            {
                "dtype": int,
                "default": 1000,
                "help": "Maximum number of steps used in the Quasi-Newton iteration scheme.",
            },
        ),
        "norm_order": (
            InputValue,
            {
                "dtype": int,
                "default": 2,
                "help": "Order of norm used to determine termination of the Quasi-newton iteration.",
            },
        ),
    }
    default_help = (
        "Holds all parameters for the numerical method used to solve the contraint."
    )
    default_label = "CSOLVER"

    def store(self, csolver):
        self.tolerance.store(csolver.tolerance)
        self.maxit.store(csolver.maxit)
        self.norm_order.store(csolver.norm_order)

    def fetch(self):
        return super(InputConstraintSolver, self).fetch()


class InputConstraintBase(Input):
    """
    An input class to define constraints. ATM built just for bonds,
    but general enough to be extended. Also incorporate a "multi"
    mode that makes it possible to define a set of constraints that
    are coupled and should hence be handled simultaneously.
    """

    attribs = {
        "mode": (
            InputAttribute,
            {
                "dtype": str,
                "default": "distance",
                "help": "The type of constraint. ",
                "options": ["distance", "angle", "eckart", "multi"],
            },
        )
    }

    fields = {
        "atoms": (
            InputArray,
            {
                "dtype": int,
                "default": np.zeros(0, int),
                "help": "List of atoms indices that are to be constrained.",
            },
        ),
        "values": (
            InputArray,
            {
                "dtype": float,
                "default": np.zeros(0, int),
                "dimension": "length",
                "help": "List of constraint lengths.",
            },
        ),
    }

    def store(self, cnstr):
        if type(cnstr) is RigidBondConstraint:
            self.mode.store("distance")
            self.atoms.store(cnstr.constrained_indices)
            self.values.store(cnstr.constraint_values)
        if type(cnstr) is AngleConstraint:
            self.mode.store("angle")
            self.atoms.store(cnstr.constrained_indices)
            self.values.store(cnstr.constraint_values)
        if type(cnstr) is EckartConstraint:
            self.mode.store("eckart")
            self.atoms.store(cnstr.constrained_indices)
            # NOTE: this is special
            self.values.store(cnstr.qref.flatten())

    def fetch(self):
        if self.mode.fetch() == "distance":
            alist = self.atoms.fetch()
            dlist = self.values.fetch()
            if len(alist.shape) == 1:
                alist.shape = (alist.shape[0] // 2, 2)
            if len(dlist) != len(alist) and len(dlist) != 0:
                raise ValueError(
                    "Length of atom indices and of distance list do not match"
                )
            robj = RigidBondConstraint(alist, dlist)
        elif self.mode.fetch() == "angle":
            alist = self.atoms.fetch()
            dlist = self.values.fetch()
            if len(alist.shape) == 1:
                alist.shape = (alist.shape[0] // 3, 3)
            if len(dlist) != len(alist) and len(dlist) != 0:
                raise ValueError(
                    "Length of atom indices and of distance list do not match"
                )
            robj = AngleConstraint(alist, dlist)
        elif self.mode.fetch() == "eckart":
            alist = self.atoms.fetch()
            dlist = self.values.fetch()
            alist.shape = -1
            if len(dlist) != 3 * len(alist) and len(dlist) != 0:
                raise ValueError(
                    "Length of atom indices and of list of coordinates do not match"
                )
            robj = EckartConstraint(alist, dlist)

        return robj


class InputConstraint(InputConstraintBase):
    attribs = copy(InputConstraintBase.attribs)

    attribs["mode"][1]["options"].append("multi")

    dynamic = {
        "constraint": (
            InputConstraintBase,
            {"help": "One or more constraints that have to be considered coupled"},
        )
    }

    def store(self, cnstr):
        if type(cnstr) is ConstraintList:
            self.mode.store("multi")
            self.extra = []
            for constr in cnstr.constraint_list:
                iobj = InputConstraint()
                iobj.store(constr)
                self.extra.append(("constraint", iobj))
        else:
            super(InputConstraint, self).store(cnstr)

    def fetch(self):
        if self.mode.fetch() == "multi":
            cnstr_list = []
            for k, v in self.extra:
                if k == "constraint":
                    cnstr_list.append(v.fetch())
                else:
                    raise ValueError("Invalid entry " + k + " in constraint(multi)")
            robj = ConstraintList(cnstr_list)
        else:
            robj = super(InputConstraint, self).fetch()
        return robj


class InputConstrainedDynamics(InputDictionary):
    """Dynamics input class.

    Handles generating the appropriate ensemble class from the xml input file,
    and generating the xml checkpoint tags and data from an instance of the
    object.

    Attributes:
        mode: An optional string giving the mode (ensemble) to be simulated.
            Defaults to 'unknown'.

    Fields:
        thermostat: The thermostat to be used for constant temperature dynamics.
        barostat: The barostat to be used for constant pressure or stress
            dynamics.
        timestep: An optional float giving the size of the timestep in atomic
            units. Defaults to 1.0.
    """

    attribs = {
        "mode": (
            InputAttribute,
            {
                "dtype": str,
                "default": "nve",
                "help": "The ensemble that will be sampled during the simulation. ",
                "options": ["nve", "nvt"],
            },
        ),
        "splitting": (
            InputAttribute,
            {
                "dtype": str,
                "default": "baoab",
                "help": "The integrator used for sampling the target ensemble. ",
                "options": ["obabo", "baoab"],
            },
        ),
    }

    fields = {
        "thermostat": (
            InputThermo,
            {
                "default": input_default(factory=ipi.engine.thermostats.Thermostat),
                "help": "The thermostat for the atoms, keeps the atom velocity distribution at the correct temperature.",
            },
        ),
        "barostat": (
            InputBaro,
            {
                "default": input_default(factory=ipi.engine.barostats.Barostat),
                "help": InputBaro.default_help,
            },
        ),
        "timestep": (
            InputValue,
            {
                "dtype": float,
                "default": 1.0,
                "help": "The time step.",
                "dimension": "time",
            },
        ),
        "nmts": (
            InputArray,
            {
                "dtype": int,
                "default": np.zeros(0, int),
                "help": "Number of iterations for each MTS level (including the outer loop, that should in most cases have just one iteration).",
            },
        ),
        "nsteps_o": (
            InputValue,
            {
                "dtype": int,
                "default": 1,
                "help": "The number of sub steps used in the evolution of the thermostat (used in function step_Oc). Relevant only for GLE thermostats",
            },
        ),
        "nsteps_geo": (
            InputValue,
            {
                "dtype": int,
                "default": 1,
                "help": "The number of sub steps used in the evolution of the geodesic flow (used in function step_Ag).",
            },
        ),
        "csolver": (
            InputConstraintSolver,
            {
                "help": "Define a numerical method for computing the projection operators associated with the constraint."
            },
        ),
    }

    dynamic = {
        "constraint": (
            InputConstraint,
            {"help": "Define a constraint to be applied onto atoms"},
        )
    }

    default_help = "Holds all the information for the MD integrator, such as timestep, the thermostats and barostats that control it."
    default_label = "CONSTRAINEDDYNAMICS"

    def store(self, dyn):
        """Takes an ensemble instance and stores a minimal representation of it.

        Args:
            dyn: An integrator object.
        """

        if dyn == {}:
            return

        self.mode.store(dyn.enstype)
        self.timestep.store(dyn.dt)
        self.thermostat.store(dyn.thermostat)
        self.barostat.store(dyn.barostat)
        self.nmts.store(dyn.nmts)
        self.splitting.store(dyn.splitting)
        self.nsteps_o.store(dyn.nsteps_o)
        self.nsteps_geo.store(dyn.nsteps_geo)
        self.csolver.store(dyn.csolver)

        self.extra = []

        for constr in dyn.constraint_list:
            iobj = InputConstraint()
            iobj.store(constr)
            self.extra.append(("constraint", iobj))

    def fetch(self):
        """Creates a ConstrainedDynamics object.

        Returns:
            An ensemble object of the appropriate mode and with the appropriate
            objects given the attributes of the InputEnsemble object.
        """

        rv = super(InputConstrainedDynamics, self).fetch()
        rv["mode"] = self.mode.fetch()
        rv["splitting"] = self.splitting.fetch()
        rv["csolver"] = self.csolver.fetch()

        cnstr_list = []
        for k, v in self.extra:
            if k == "constraint":
                cnstr_list.append(v.fetch())
            else:
                raise ValueError("Invalid entry " + k + " in constrained_dynamics")

        rv["constraint_list"] = cnstr_list

        return rv

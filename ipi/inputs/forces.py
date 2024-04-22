"""Creates objects that compose and apply forces."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


from ipi.engine.forces import *
from ipi.utils.inputvalue import *
import numpy as np

__all__ = ["InputForces", "InputForceComponent"]


class InputForceComponent(Input):
    """ForceComponent input class.

    Uses the forcefield object whose name is specified as the value of the
    field (matching one of the forcefields defined in the simulation tag)
    to compute one component of the force acting on the ring polymer.


    Attributes:
       nbeads: The number of beads that the forcefield will be evaluated on.
       weight: A scaling factor for the contribution from this forcefield.
       name: The name of the forcefield.
    """

    attribs = {
        "nbeads": (
            InputAttribute,
            {
                "dtype": int,
                "default": 0,
                "help": "If the forcefield is to be evaluated on a contracted ring polymer, this gives the number of beads that are used. If not specified, the forcefield will be evaluated on the full ring polymer.",
            },
        ),
        "weight": (
            InputAttribute,
            {
                "dtype": float,
                "default": 1.0,
                "help": "A scaling factor for this forcefield, to be applied before adding the force calculated by this forcefield to the total force.",
            },
        ),
        "fd_epsilon": (
            InputAttribute,
            {
                "dtype": float,
                "default": -0.001,
                "help": "The finite displacement to be used for calculaing the Suzuki-Chin contribution of the force. If the value is negative, a centered finite-difference scheme will be used. [in bohr]",
            },
        ),
        "name": (
            InputAttribute,
            {
                "dtype": str,
                "default": "",
                "help": "An optional name to refer to this force component.",
            },
        ),
        "forcefield": (
            InputAttribute,
            {
                "dtype": str,
                "default": "",
                "help": "Mandatory. The name of the forcefield this force is referring to.",
            },
        ),
    }

    fields = {
        "mts_weights": (
            InputArray,
            {
                "dtype": float,
                "default": np.zeros(1, float) + 1.0,
                "help": "The weight of force in each mts level starting from outer.",
                "dimension": "force",
            },
        ),
        "interpolate_extras": (
            InputArray,
            {
                "dtype": str,
                "default": np.zeros(0, str),
                "help": """A list of quantities that should be extracted from the 'extra' string returned by the forcefield,
                           that are then treated the same way as energy or forces -- that is treated as a numerical, physical
                           quantity, interpolated when changing the number of PI replicas. Same quantities from different force
                           components are summed as well. The names should correspond to entries in the JSON-formatted extra string.""",
            },
        ),
    }

    default_help = "The class that deals with how each forcefield contributes to the overall potential, force and virial calculation."
    default_label = "FORCECOMPONENT"

    def store(self, forceb):
        """Takes a ForceComponent instance and stores a minimal
        representation of it.

        Args:
           forceb: A ForceComponent object.
        """

        super(InputForceComponent, self).store(forceb.ffield)
        self.nbeads.store(forceb.nbeads)
        self.weight.store(forceb.weight)
        self.mts_weights.store(forceb.mts_weights)
        self.interpolate_extras.store(forceb.interpolate_extras)
        self.fd_epsilon.store(forceb.epsilon)
        self.name.store(forceb.name)
        self.forcefield.store(forceb.ffield)

    def fetch(self):
        """Creates a ForceComponent object.

        Returns:
           A ForceComponent object.
        """

        super(InputForceComponent, self).fetch()
        return ForceComponent(
            ffield=self.forcefield.fetch(),
            nbeads=self.nbeads.fetch(),
            weight=self.weight.fetch(),
            name=self.name.fetch(),
            mts_weights=self.mts_weights.fetch(),
            interpolate_extras=self.interpolate_extras.fetch(),
            epsilon=self.fd_epsilon.fetch(),
        )

    def check(self):
        """Checks for optional parameters."""

        super(InputForceComponent, self).check()
        if self.nbeads.fetch() < 0:
            raise ValueError(
                "The forces must be evaluated over a positive number of beads."
            )


class InputForces(Input):
    """Deals with creating all the forcefield objects required in the
    simulation.

    Dynamic fields:
       socket: Socket object to create the server socket.
    """

    # At the moment only socket driver codes implemented, other types
    # could be used in principle
    dynamic = {
        "force": (InputForceComponent, {"help": InputForceComponent.default_help})
    }

    default_help = "Deals with creating all the necessary forcefield objects."
    default_label = "FORCES"

    def fetch(self):
        """Returns a list of the output objects included in this dynamic
        container.

        Returns:
           A list of tuples, with each tuple being of the form ('type', 'object'),
           where 'type' is the type of forcefield, and 'object' is a
        """

        super(InputForces, self).fetch()
        flist = [f.fetch() for (n, f) in self.extra]
        return flist

    def store(self, flist):
        """Stores a list of the output objects, creating a sequence of
        dynamic containers.

        Args:
           flist: A list of tuples, with each tuple being of the form
           ('type', 'object') where 'type' is the type of forcefield
           and 'object' is a forcefield object of that type.
        """

        super(InputForces, self).store()

        if len(self.extra) != len(flist):
            self.extra = [0] * len(flist)

        for ii, el in enumerate(flist):
            if self.extra[ii] == 0:
                iff = InputForceComponent()
                iff.store(el)
                self.extra[ii] = ("force", iff)
            else:
                self.extra[ii][1].store(el)

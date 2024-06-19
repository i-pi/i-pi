"""Creates objects that deal with the different ensembles."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

import numpy as np
import ipi.engine.thermostats
import ipi.engine.barostats
from ipi.utils.inputvalue import (
    InputDictionary,
    InputAttribute,
    InputValue,
    InputArray,
    input_default,
)
from ipi.inputs.thermostats import InputThermo


__all__ = ["InputPlanetary"]


class InputPlanetary(InputDictionary):
    """Planetary input class.

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
                "default": "md",
                "help": "The constrained-centroid sampling mode. ",
                "options": ["md"],
            },
        )
    }

    fields = {
        "thermostat": (
            InputThermo,
            {
                "default": input_default(factory=ipi.engine.thermostats.Thermostat),
                "help": "The thermostat for the atoms, keeps the atom velocity distribution at the correct temperature.",
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
        "nsamples": (
            InputValue,
            {
                "dtype": int,
                "default": 0,
                "help": "Number of samples to accumulate for each planetary step.",
            },
        ),
        "stride": (
            InputValue,
            {
                "dtype": int,
                "default": 1,
                "help": "How often the planetary calculation should actually be triggered.",
            },
        ),
        "nbeads": (
            InputValue,
            {
                "dtype": int,
                "default": -1,
                "help": "Number of beads for centroid-constrained dynamics (default same as main trajectory)",
            },
        ),
        "screen": (
            InputValue,
            {
                "dtype": float,
                "default": 0.0,
                "dimension": "length",
                "help": "Screening parameter for path-integral frequency matrix.",
            },
        ),
    }

    dynamic = {}

    default_help = (
        "Holds all the information for the planetary model frequency matrix calculator."
    )
    default_label = "PLANETARY"

    def store(self, plan):
        """Takes a planetary instance and stores a minimal representation of it.

        Args:
            dyn: An integrator object.
        """

        if plan == {}:
            return

        self.mode.store(plan.mode)
        self.timestep.store(plan.ccdyn.dt)
        self.thermostat.store(plan.ccdyn.thermostat)
        self.nmts.store(plan.ccdyn.nmts)
        self.nsamples.store(plan.nsamples)
        self.stride.store(plan.stride)
        self.screen.store(plan.screen)

    def fetch(self):
        """Creates an ensemble object.

        Returns:
            An ensemble object of the appropriate mode and with the appropriate
            objects given the attributes of the InputEnsemble object.
        """

        rv = super(InputPlanetary, self).fetch()
        rv["mode"] = self.mode.fetch()
        return rv

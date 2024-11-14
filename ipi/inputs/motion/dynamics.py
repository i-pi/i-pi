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
from ipi.inputs.barostats import InputBaro
from ipi.inputs.thermostats import InputThermo


__all__ = ["InputDynamics"]


class InputDynamics(InputDictionary):
    """Dynamics input class.

    Handles generating the appropriate ensemble class from the xml input file,
    and generating the xml checkpoint tags and data from an instance of the
    object.

    Attributes:
        mode: An optional string giving the mode (ensemble) to be simulated.
            Defaults to 'nve'.
        splitting: An optional string giving the Louiville splitting used for
            sampling the target ensemble. Defaults to 'obabo'.

    Fields:
        thermostat: The thermostat to be used for constant temperature dynamics.
        barostat: The barostat to be used for constant pressure or stress
            dynamics.
        timestep: An optional float giving the size of the timestep in atomic
            units. Defaults to 1.0.
        nmts: Number of iterations for each MTS level
    """

    attribs = {
        "mode": (
            InputAttribute,
            {
                "dtype": str,
                "default": "nve",
                "help": """The ensemble that will be sampled during the simulation.
                nve: constant-energy-volume; nvt: constant-temperature-volume;
                npt: constant-temperature-pressure(isotropic); nst: constant-temperature-stress(anisotropic);
                sc: Suzuki-Chin high-order NVT; scnpt: Suzuki-Chin high-order NpT;
                nvt-cc: constrained-centroid NVT;
                 """,
                "options": [
                    "nve",
                    "nvt",
                    "npt",
                    "nst",
                    "sc",
                    "scnpt",
                    "nvt-cc",
                ],
            },
        ),
        "splitting": (
            InputAttribute,
            {
                "dtype": str,
                "default": "obabo",
                "help": "The Louiville splitting used for sampling the target ensemble. ",
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
    }

    dynamic = {}

    default_help = "Holds all the information for the MD integrator, such as timestep, the thermostats and barostats that control it."
    default_label = "DYNAMICS"

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

    def fetch(self):
        """Creates an ensemble object.

        Returns:
            An ensemble object of the appropriate mode and with the appropriate
            objects given the attributes of the InputEnsemble object.
        """

        rv = super(InputDynamics, self).fetch()
        rv["mode"] = self.mode.fetch()
        rv["splitting"] = self.splitting.fetch()
        return rv

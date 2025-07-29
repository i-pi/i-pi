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
                    "nve-f",
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
        options = dyn.options
        self.friction.store(options["friction"])
        self.frictionSD.store(options["frictionSD"])
        if options["friction"]:
            self.fric_spec_dens.store(options["fric_spec_dens"])
            self.fric_spec_dens_ener.store(options["fric_spec_dens_ener"])

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

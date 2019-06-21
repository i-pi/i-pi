"""Creates objects that deal with the different ensembles."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

import numpy as np
import ipi.engine.thermostats
import ipi.engine.barostats
from ipi.utils.inputvalue import InputDictionary, InputAttribute, InputValue, InputArray, input_default
from ipi.inputs.barostats import InputBaro
from ipi.inputs.thermostats import InputThermo


__all__ = ['InputConstrainedDynamics']


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
        "mode": (InputAttribute, {"dtype": str,
                                  "default": 'nve',
                                  "help": "The ensemble that will be sampled during the simulation. ",
                                  "options": ['nve', 'nvt']}),
        "splitting": (InputAttribute, {"dtype": str,
                                       "default": 'baoab',
                                       "help": "The integrator used for sampling the target ensemble. ",
                      "options": ['rattle','geodesic','gobabo', 'baoab']})
    }

    fields = {
        "thermostat": (InputThermo, {"default": input_default(factory=ipi.engine.thermostats.Thermostat),
                                     "help": "The thermostat for the atoms, keeps the atom velocity distribution at the correct temperature."}),
        "barostat": (InputBaro, {"default": input_default(factory=ipi.engine.barostats.Barostat),
                                 "help": InputBaro.default_help}),
        "timestep": (InputValue, {"dtype": float,
                                  "default": 1.0,
                                  "help": "The time step.",
                                  "dimension": "time"}),
        "nmts": (InputArray, {"dtype": int,
                              "default": np.zeros(0, int),
                              "help": "Number of iterations for each MTS level (including the outer loop, that should in most cases have just one iteration)."}),
        "nsteps_o": (InputValue, {"dtype": int,
                                "default": 1,
                                "help": "The number of sub steps used in the evolution of the thermostat (used in function step_Oc). Relevant only for GLE thermostats" }),
        "nsteps_geo": (InputValue, {"dtype": int,
                                           "default": 1,
                                           "help": "The number of sub steps used in the evolution of the geodesic flow (used in function step_Ag)." }),
        "constrained_indices": (InputArray, {"dtype": int,
                              "default": np.zeros(0, int),
                              "help": "List of the form [i,j] containing the list of atoms indices that are to be constrained."}),
        "constrained_distances": (InputArray, {"dtype": float,
                              "default": np.zeros(0, float),
                              "help": "List of the form [rij] containing the bond distances."})
    }

    dynamic = {}

    default_help = "Holds all the information for the MD integrator, such as timestep, the thermostats and barostats that control it."
    default_label = "CONSTRAINED_DYNAMICS"

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
        self.constraint_list.store(dyn.constraint_list)
        self.splitting.store(dyn.splitting)

    def fetch(self):
        """Creates an ensemble object.

        Returns:
            An ensemble object of the appropriate mode and with the appropriate
            objects given the attributes of the InputEnsemble object.
        """

        rv = super(InputConstrainedDynamics, self).fetch()
        rv["mode"] = self.mode.fetch()
        rv["splitting"] = self.splitting.fetch()
        return rv

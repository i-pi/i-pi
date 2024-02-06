"""Creates objects that deal with the different ensembles."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np

from ipi.inputs.forces import InputForces
from ipi.engine.ensembles import *
from ipi.utils.inputvalue import *
from ipi.utils.units import *


__all__ = ["InputEnsemble"]


class InputEnsemble(Input):
    """Ensemble input class.

    Handles generating the appropriate ensemble class from the xml input file,
    and generating the xml checkpoint tags and data from an instance of the
    object.

    Attributes:
       mode: An optional string giving the mode of ensemble to be simulated.
          Defaults to 'unknown'.

    Fields:
       temperature: An optional float giving the temperature in atomic units.
          Defaults to 1.0.
       pressure: An optional float giving the external pressure in atomic units.
          Defaults to 1.0.
       eens: An optional float giving the ensemble contribution to the conserved
          quantity.
       stress: An optional array containing the terms of the stress tensor as
          [pxx, pxy, pxz, pyx, pyy .. pzy, pzz].
    """

    fields = {
        "temperature": (
            InputValue,
            {
                "dtype": float,
                "default": -1.0,
                "help": "The temperature of the system.",
                "dimension": "temperature",
            },
        ),
        "pressure": (
            InputValue,
            {
                "dtype": float,
                "default": -12345,  # hard-coded to signal unset pressure
                "help": "The external pressure.",
                "dimension": "pressure",
            },
        ),
        "stress": (
            InputArray,
            {
                "dtype": float,
                "default": -12345.0
                * np.identity(3, float),  # hard-coded to signal unset stress
                "help": "The external stress.",
                "dimension": "pressure",
            },
        ),
        "eens": (
            InputValue,
            {
                "dtype": float,
                "default": 0.0,
                "help": "The ensemble contribution to the conserved quantity.",
                "dimension": "energy",
            },
        ),
        "bias": (InputForces, {"help": InputForces.default_help, "default": []}),
        "bias_weights": (
            InputArray,
            {
                "dtype": float,
                "default": np.zeros(0),
                "help": "Bias weights.",
                "dimension": "undefined",
            },
        ),
        "hamiltonian_weights": (
            InputArray,
            {
                "dtype": float,
                "default": np.zeros(0),
                "help": "Hamiltonian weights.",
                "dimension": "undefined",
            },
        ),
        "time": (
            InputValue,
            {
                "dtype": float,
                "default": 0.0,
                "dimension": "time",
                "help": "The internal time for this system",
            },
        ),
    }
    dynamic = {}

    default_help = "Holds all the information that is ensemble specific, such as the temperature and the external pressure."
    default_label = "ENSEMBLE"

    def store(self, ens):
        """Takes an ensemble instance and stores a minimal representation of it.

        Args:
           ens: An ensemble object.
        """

        super(InputEnsemble, self).store(ens)
        self.temperature.store(ens.temp)
        self.pressure.store(ens.pext)
        self.stress.store(ens.stressext)
        self.eens.store(ens.eens)
        self.bias.store(ens.bcomp)
        self.bias_weights.store(ens.bweights)
        self.hamiltonian_weights.store(ens.hweights)
        self.time.store(ens.time)

    def fetch(self):
        """Creates an ensemble object.

        Returns:
           An ensemble object of the appropriate mode and with the appropriate
           objects given the attributes of the InputEnsemble object.
        """

        super(InputEnsemble, self).fetch()

        ens = Ensemble(
            eens=self.eens.fetch(),
            temp=self.temperature.fetch(),
            pext=self.pressure.fetch(),
            stressext=self.stress.fetch(),
            bcomponents=self.bias.fetch(),
            bweights=self.bias_weights.fetch(),
            hweights=self.hamiltonian_weights.fetch(),
            time=self.time.fetch(),
        )

        return ens

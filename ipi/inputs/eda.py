"""Creates objects that deal with constant pressure and stress simulations."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np

from ipi.engine.eda import ElectricField, BEC
from ipi.engine.barostats import *
from ipi.utils.inputvalue import *
from ipi.inputs.thermostats import *
from ipi.engine.cell import Cell
from ipi.inputs.cell import *
from copy import copy


__all__ = ["InputElectricField", "InputBEC"]


class InputElectricField(Input):

    """Electric field input class."""

    attribs = {}

    fields = {
        "amp": (
            InputArray,
            {
                "dtype": float,
                "default": np.zeros(3),
                "help": "The amplitude of the external electric field (in cartesian coordinates)",
                "dimension": "electric-field",
            },
        ),
        "freq": (
            InputValue,
            {
                "dtype": float,
                "default": 0.0,
                "help": "The pulsation of the external electric field",
                "dimension": "frequency",
            },
        ),
        "phase": (
            InputValue,
            {
                "dtype": float,
                "default": 0.0,
                "help": "The phase of the external electric field (in rad)",
            },
        ),
        "peak": (
            InputValue,
            {
                "dtype": float,
                "default": 0.0,
                "help": "The time when the external electric field gets its maximum value",
                "dimension": "time",
            },
        ),
        "sigma": (
            InputValue,
            {
                "dtype": float,
                "default": np.inf,
                "help": "The standard deviations (time) of the gaussian envelope function of the external electric field",
                "dimension": "time",
            },
        ),
    }

    default_help = "Simulates an external time dependent electric field"
    default_label = "Efield"

    def __init__(self, *argv, **kargw):
        super().__init__(*argv, **kargw)
        pass

    def store(self, Efield):
        for k in self.fields.keys():
            getattr(self, k).store(getattr(Efield, k))
        return

    def fetch(self):
        return ElectricField(
            amp=self.amp.fetch(),
            freq=self.freq.fetch(),
            phase=self.phase.fetch(),
            peak=self.peak.fetch(),
            sigma=self.sigma.fetch(),
        )


class InputBEC(InputTensor):

    """Born Effective Charges field input class."""

    attribs = copy(InputTensor.attribs)

    default_help = "Simulates an external time dependent electric field"
    default_label = "Efield"

    # def __init__(self,*argv, **kargw):
    #     super().__init__(*argv, **kargw)
    #     pass

    def store(self, BEC):
        for k in self.fields.keys():
            getattr(self, k).store(getattr(BEC, k))
        return

    def fetch(self):
        super().fetch()
        return BEC(cbec=self.mode == "otf", bec=self.value.reshape((-1, 3)))

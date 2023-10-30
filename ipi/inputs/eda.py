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
from ipi.utils.depend import dd, dpipe
from ipi.utils.depend import dobject, depend_array, depend_value


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

    # def __init__(self, *argv, **kargw):
    #     super().__init__(*argv, **kargw)
    #     pass

    def store(self, Efield: ElectricField):
        self.amp.store(Efield.amp)
        self.freq.store(Efield.freq)
        self.phase.store(Efield.phase)
        self.peak.store(Efield.peak)
        self.sigma.store(Efield.sigma)
        # for k in self.fields.keys():
        #     getattr(self, k).store(getattr(Efield, k))
        return

    def fetch(self):
        return ElectricField(
            amp=self.amp.fetch(),
            freq=self.freq.fetch(),
            phase=self.phase.fetch(),
            peak=self.peak.fetch(),
            sigma=self.sigma.fetch(),
        )


class InputBEC(Input):
    attribs = copy(InputArray.attribs)

    attribs["mode"] = (
        InputAttribute,
        {
            "dtype": str,
            "default": "none",
            "options": ["driver", "manual", "file", "none"],
            "help": "If 'mode' is 'DFPT', then the array is computed on the fly. "
            + InputArray.attribs["mode"][1]["help"],
        },
    )

    fields = {
        "bec": (
            InputArray,
            {
                "dtype" : float,
                "default": np.zeros(0),
            },
        )
    }

    def parse(self, xml=None, text=""):
        """Reads the data for an array from an xml file.

        Args:
           xml: An xml_node object containing the all the data for the parent
              tag.
           text: The data held between the start and end tags.
        """
        Input.parse(self, xml=xml, text=text)
        mode = self.mode.fetch()
        if mode in ["manual", "file"]:  # ['manual','file']
            self.bec.parse(xml, text)
        elif mode == "none":
            self.bec.value = np.full((0, 0), np.nan)
        elif mode == "driver":
            self.bec.value = np.full((0, 0), np.nan)
        else:
            raise ValueError("error in InputBEC.parse")

    def fetch(self):
        super().fetch()
        return BEC(
            cbec=self.mode.fetch() == "driver", bec=self.bec.fetch().reshape((-1, 3))
        )

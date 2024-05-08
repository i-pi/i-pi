"""Creates objects to deal with dynamics driven by external fields."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

from ipi.engine.motion.driven_dynamics import ElectricField, BEC
from ipi.utils.inputvalue import (
    input_default,
)
from ipi.inputs.motion.dynamics import InputDynamics
import numpy as np

from ipi.engine.barostats import *
from ipi.utils.inputvalue import *
from ipi.inputs.thermostats import *
from ipi.inputs.cell import *
from copy import copy


__all__ = ["InputDrivenDynamics", "InputElectricField", "InputBEC"]


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
                "dimension": "number",
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

    def store(self, Efield: ElectricField):
        self.amp.store(Efield.amp)
        self.freq.store(Efield.freq)
        self.phase.store(Efield.phase)
        self.peak.store(Efield.peak)
        self.sigma.store(Efield.sigma)
        return

    def fetch(self):
        return ElectricField(
            amp=self.amp.fetch(),
            freq=self.freq.fetch(),
            phase=self.phase.fetch(),
            peak=self.peak.fetch(),
            sigma=self.sigma.fetch(),
        )


class InputBEC(InputArray):
    """BEC input class.

    Class to set the Born Effective Charges.
    The rationale oif the functioning of this class is the same as in InputCell.
    """

    attribs = copy(InputArray.attribs)

    attribs["mode"] = (
        InputAttribute,
        {
            "dtype": str,
            "default": "none",
            "options": ["driver", "manual", "file", "none"],
            "help": "If 'mode' is 'driver', then the array is computed on the fly. "
            + InputArray.attribs["mode"][1]["help"],
        },
    )

    default_help = "Deals with the Born Effective Charges tensors"
    default_label = "BEC"

    def __init__(self, help=None, dimension=None, units=None, default=None, dtype=None):
        """Initializes InputBEC.

        Just calls the parent initialization function with appropriate arguments.
        """
        super().__init__(help=help, dimension="number", default=default, dtype=float)

    def store(self, bec):
        super(InputBEC, self).store(bec.bec)

    def parse(self, xml=None, text=""):
        """Reads the data for an array from an xml file.

        Args:
           xml: An xml_node object containing the all the data for the parent
              tag.
           text: The data held between the start and end tags.
        """
        Input.parse(self, xml=xml, text=text)
        mode = self.mode.fetch()
        if mode in ["manual", "file"]:
            super().parse(xml, text)
        elif mode == "none":
            self.value = np.full((0, 3), np.nan)
        elif mode == "driver":
            self.value = np.full((0, 3), np.nan)
        else:
            raise ValueError("error in InputBEC.parse")

    def fetch(self):
        bec = super(InputBEC, self).fetch()
        return BEC(cbec=self.mode.fetch() == "driver", bec=bec.reshape((-1, 3)))


class InputDrivenDynamics(InputDynamics):
    """Driven Dynamics input class.

    Attributes: same of InputDynamics

    Fields: same of InputDynamics

    """

    fields = {
        "efield": (
            InputElectricField,
            {
                "default": input_default(factory=ElectricField),
                "help": "The external electric field parameters:"
                + "plane-wave parameters (intensity/amplitude, angular frequency, and phase) and "
                + "gaussian envelope function parameters (peak time/mean of the gaussian, and pulse duration/standard deviation of the gaussian)",
            },
        ),
        "bec": (
            InputBEC,
            {
                "dtype": float,
                "default": input_default(factory=BEC),
                "dimension": "number",
                "help": "The Born Effective Charges tensors (cartesian coordinates)",
            },
        ),
    }

    fields.update(InputDynamics.fields)

    attribs = {
        "mode": (
            InputAttribute,
            {
                "dtype": str,
                "default": "eda-nve",
                "help": """The ensemble that will be sampled during the simulation.
                eda-nve: nve with an external electric field;
                 """,
                "options": [
                    "eda-nve",
                ],
            },
        ),
        "splitting": InputDynamics.attribs["splitting"],
    }

    dynamic = {}

    default_help = "Holds all the information for a driven dynamics."
    default_label = "DRIVEN_DYNAMICS"

    def __init__(self, *argc, **argv):
        super().__init__(*argc, **argv)

    def store(self, dyn):
        if dyn == {}:
            return
        super().store(dyn)
        self.efield.store(dyn.Electric_Field)
        self.bec.store(dyn.Born_Charges)

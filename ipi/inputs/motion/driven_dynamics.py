"""Creates objects that deal with the different ensembles."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

from ipi.engine.motion.eda import ElectricField, BEC
from ipi.utils.inputvalue import (
    input_default,
)
from ipi.inputs.motion.eda import InputElectricField, InputBEC
from ipi.inputs.motion.dynamics import InputDynamics


__all__ = ["InputDrivenDynamics"]


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

    attribs = {}
    attribs.update(InputDynamics.attribs)

    dynamic = {}

    default_help = "Holds all the information for a driven dynamics."
    default_label = "DRIVEN_DYNAMICS"

    def __init__(self, *argc, **argv):
        super().__init__(*argc, **argv)

    def store(self, dyn):
        if dyn == {}:
            return
        super().store(dyn)
        self.efield.store(dyn.efield)
        self.bec.store(dyn.bec)

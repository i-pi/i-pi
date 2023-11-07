import numpy as np
from ipi.engine.motion import *
from ipi.utils.inputvalue import *
from ipi.inputs.thermostats import *
from ipi.inputs.initializer import *
from ipi.utils.units import *

__all__ = ["InputBECTensorsCalculator"]


class InputBECTensorsCalculator(InputDictionary):

    """
    Contains options related with finite difference computation of Born Effective Charges.
    """

    attribs = {}
    fields = {
        "pos_shift": (
            InputValue,
            {
                "dtype": float,
                "default": 0.01,
                "dimension": "length",
                "help": "The finite displacement in position used to compute derivative of the polarization.",
            },
        ),
        "prefix": (
            InputValue,
            {"dtype": str, "default": "BEC", "help": "Prefix of the output files."},
        ),
        "asr": (
            InputValue,
            {
                "dtype": str,
                "default": "none",
                "options": ["none", "lin"],  # , "poly", "lin", "crystal"],
                "help": "Removes the zero frequency vibrational modes depending on the symmerty of the system.",
            },
        ),
        "bec": (
            InputArray,
            {
                "dtype": float,
                "default": np.zeros(0, float),
                "dimension": "number",
                "help": "Portion of the BECs known up to now.",
            },
        ),
        "atoms": (
            InputArray,
            {
                "dtype": str,
                "default": np.asarray(["all"]),
                "help": "Atoms whose BEC tensor has to be computed. It can be 'all', a chemical species ('Li', 'Mg') or an atom index. List of the previous cases are accepted.",
            },
        ),
    }

    dynamic = {}

    default_help = "Fill in."
    default_label = "BEC"

    def store(self, phonons):
        if phonons == {}:
            return
        self.pos_shift.store(phonons.deltax)
        self.prefix.store(phonons.prefix)
        self.asr.store(phonons.asr)
        self.bec.store(phonons.bec)
        self.atoms.store(phonons.atoms)

    def fetch(self):
        rv = super(InputBECTensorsCalculator, self).fetch()
        return rv

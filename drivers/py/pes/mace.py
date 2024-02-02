""" An interface for the [MACE](https://github.com/ACEsuit/mace) calculator """

import sys
from .ase import ASEDriver

try:
    from mace.calculators import MACECalculator
except:
    MACECalculator = None

__DRIVER_NAME__ = "mace"
__DRIVER_CLASS__ = "MACE_driver"

ERROR_MSG = """
MACE driver requires specification of a .json model,
and a template file that describes the chemical makeup of the structure.

Example: python driver.py -m mace -u -o model.json,template.xyz
"""


class MACE_driver(ASEDriver):
    def __init__(self, args=None, verbose=False):
        if MACECalculator is None:
            raise ImportError("Couldn't load mace bindings")

        super().__init__(args, verbose, ERROR_MSG)

    def check_arguments(self):
        """Check the arguments requuired to run the driver

        This loads the potential and atoms template in MACE
        """

        super().check_arguments()

        if len(self.args) < 2:
            sys.exit(self.error_msg)

        self.ase_calculator = MACECalculator(model_paths=self.args[1], device="cpu")

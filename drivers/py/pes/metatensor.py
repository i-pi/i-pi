"""
Interface with metatensor
(https://lab-cosmo.github.io/metatensor/latest/atomistic/index.html), that can
be used to perform calculations based on different types of machine learning
potentials
"""

import sys

from ipi.utils.messages import warning

from .ase import ASEDriver

try:
    from metatensor.torch.atomistic.ase_calculator import MetatensorCalculator
except ImportError:
    warning("Could not find or import metatensor.torch")
    MetatensorCalculator = None

__DRIVER_NAME__ = "metatensor"
__DRIVER_CLASS__ = "MetatensorDriver"

ERROR_MSG = """
The metatensor driver requires specification of a .pt TorchScript model and an
ASE-readable template file that describes the chemical makeup of the structure.

Example: python driver.py -m metatensor -u -o model.pt,template.xyz
"""


class MetatensorDriver(ASEDriver):
    def __init__(self, args=None, verbose=False, error_msg=ERROR_MSG):
        super().__init__(args, verbose, error_msg)

    def check_arguments(self):
        """Check the arguments required to run the driver

        This loads the potential and atoms template in metatensor
        """

        if MetatensorCalculator is None:
            raise ImportError("could not import metatensor.torch, is it installed?")
        super().check_arguments()

        if len(self.args) < 2:
            sys.exit(self.error_msg)
        self.model_path = self.args[1]

        self.ase_calculator = MetatensorCalculator(self.model_path)

        # Show the model metadata to the users
        print(self.ase_calculator.metadata())

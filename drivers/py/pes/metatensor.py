"""Interface with the [metatensor](https://lab-cosmo.github.io/metatensor/latest/atomistic/index.html)
calculator, that can be used to perform calculations based on different types of machine learning potentials"""

import sys
from .ase import ASEDriver

from ipi.utils.messages import warning

try:
    from metatensor.torch.atomistic.ase_calculator import MetatensorCalculator
except ImportError:
    warning("Could not find or import the metatensor module")
    MetatensorCalculator = None

__DRIVER_NAME__ = "metatensor"
__DRIVER_CLASS__ = "MetatensorDriver"

ERROR_MSG = """
The metatensor driver requires specification of a .pt torchscript model 
and an ASE-readable template file that describes the chemical makeup of the structure. 

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
            raise ImportError("Couldn't load metatensor bindings")
        super().check_arguments()

        if len(self.args) < 2:
            sys.exit(self.error_msg)
        self.model_path = self.args[1]

        self.ase_calculator = MetatensorCalculator(self.model_path)

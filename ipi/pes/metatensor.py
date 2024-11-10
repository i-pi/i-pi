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
    import metatensor.torch
    from metatensor.torch.atomistic.ase_calculator import MetatensorCalculator
except ImportError as e:
    warning(f"Could not find or import metatensor.torch: {e}")
    MetatensorCalculator = None

__DRIVER_NAME__ = "metatensor"
__DRIVER_CLASS__ = "MetatensorDriver"

ERROR_MSG = """
The metatensor driver requires specification of a .pt TorchScript model and an
ASE-readable template file that describes the chemical makeup of the structure.

Example: python driver.py -m metatensor -u -o template.xyz,model.pt,device=cpu,\
extensions=path/to/extensions,check_consistency=False
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

        metatensor_major, metatensor_minor, *_ = metatensor.torch.__version__.split(".")
        metatensor_major = int(metatensor_major)
        metatensor_minor = int(metatensor_minor)

        if metatensor_major != 0 or metatensor_minor != 5:
            raise ImportError(
                "this code is only compatible with metatensor-torch v0.5.x, "
                f"found version v{metatensor.torch.__version__} "
                f"at '{metatensor.torch.__file__}'"
            )

        super().check_arguments()

        if len(self.args) < 2:
            sys.exit(self.error_msg)
        self.model_path = self.args[1]

        device = None
        extensions_directory = None
        check_consistency = False
        for arg in self.args[2:]:
            if arg.startswith("device="):
                device = arg[7:]
            elif arg.startswith("extensions="):
                extensions_directory = arg[11:]
            elif arg.startswith("check_consistency="):
                arg = arg[18:]
                if arg == "True":
                    check_consistency = True
                elif arg == "False":
                    check_consistency = False
                else:
                    raise ValueError(
                        "invalid value for check_consistency, expected True or False, "
                        f"got {arg}"
                    )

            else:
                sys.exit(self.error_msg)

        self.ase_calculator = MetatensorCalculator(
            self.model_path,
            device=device,
            extensions_directory=extensions_directory,
            check_consistency=check_consistency,
        )

        # Show the model metadata to the users
        print(self.ase_calculator.metadata())

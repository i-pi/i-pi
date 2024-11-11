"""
Interface with metatensor
(https://lab-cosmo.github.io/metatensor/latest/atomistic/index.html), that can
be used to perform calculations based on different types of machine learning
potentials
"""

from ipi.utils.messages import warning

from .ase import ASEDriver

MetatensorCalculator = None
mtt = None

__DRIVER_NAME__ = "metatensor"
__DRIVER_CLASS__ = "MetatensorDriver"


class MetatensorDriver(ASEDriver):
    _error_msg = """
The metatensor driver requires specification of a .pt TorchScript model and an
ASE-readable template file that describes the chemical makeup of the structure.

Example: python driver.py -m metatensor -u -o template.xyz,model.pt,device=cpu,\
extensions=path/to/extensions,check_consistency=False
"""

    def __init__(
        self,
        template,
        model,
        device="cpu",
        extensions="",
        check_consistency=False,
        *args,
        **kwargs,
    ):
        global MetatensorCalculator, mtt
        if MetatensorCalculator is None:
            try:
                import metatensor.torch as mtt
                from metatensor.torch.atomistic.ase_calculator import (
                    MetatensorCalculator,
                )
            except ImportError as e:
                warning(f"Could not find or import metatensor.torch: {e}")

        self.model = model
        self.device = device
        self.extensions = extensions
        self.check_consistency = check_consistency
        super().__init__(template, *args, **kwargs)

    def check_parameters(self):
        """Check the arguments required to run the driver

        This loads the potential and atoms template in metatensor
        """

        if MetatensorCalculator is None:
            raise ImportError("could not import metatensor.torch, is it installed?")

        metatensor_major, metatensor_minor, *_ = mtt.__version__.split(".")
        metatensor_major = int(metatensor_major)
        metatensor_minor = int(metatensor_minor)

        if metatensor_major != 0 or metatensor_minor < 5:
            raise ImportError(
                "this code is only compatible with metatensor-torch >= v0.5, "
                f"found version v{mtt.__version__} "
                f"at '{mtt.__file__}'"
            )

        super().check_parameters()

        self.model_path = self.model

        device = self.device
        extensions_directory = self.extensions
        check_consistency = self.check_consistency

        self.ase_calculator = MetatensorCalculator(
            self.model_path,
            device=device,
            extensions_directory=extensions_directory,
            check_consistency=check_consistency,
        )

        # Show the model metadata to the users
        print(self.ase_calculator.metadata())

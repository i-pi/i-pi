""" An interface for the [MACE](https://github.com/ACEsuit/mace) calculator """

from .ase import ASEDriver

MACECalculator = None

__DRIVER_NAME__ = "mace"
__DRIVER_CLASS__ = "MACE_driver"


class MACE_driver(ASEDriver):
    _error_msg = """
MACE driver requires specification of a .json model,
and a template file that describes the chemical makeup of the structure.

Example: python driver.py -m mace -u -o model.json,template.xyz
"""

    def __init__(self, template, model, device="cpu", *args, **kwargs):

        global MACECalculator

        try:
            from mace.calculators import MACECalculator
        except:
            raise ImportError("Couldn't load mace bindings")

        self.model = model
        self.device = device
        super().__init__(template, *args, **kwargs)

    def check_parameters(self):
        """Check the arguments requuired to run the driver

        This loads the potential and atoms template in MACE
        """

        super().check_parameters()

        self.ase_calculator = MACECalculator(model_paths=self.model, device=self.device)

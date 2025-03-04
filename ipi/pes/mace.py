"""An interface for the [MACE](https://github.com/ACEsuit/mace) calculator"""

from .ase import ASEDriver

from ipi.utils.messages import verbosity, warning

MACECalculator = None

__DRIVER_NAME__ = "mace"
__DRIVER_CLASS__ = "MACE_driver"


class MACE_driver(ASEDriver):
    """
    Driver for the MACE MLIPs.
    The driver requires specification of a .json model,
    and a template file that describes the chemical makeup
    of the structure.

    Command-line:
    i-pi-py_driver.py -m mace -u -o template=template.xyz,model=model.json

    Parameters:
    :param template: string, filename of an ASE-readable structure file
        to initialize atomic number and types
    :param model: string, filename of the json-formatted model file
    """

    def __init__(self, template, model, device="cpu", *args, **kwargs):
        warning(
            "THIS PES HAS NOT BEEN TESTED FOLLOWING CONVERSION TO THE NEW PES API.",
            verbosity.low,
        )
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

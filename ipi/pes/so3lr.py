"""An interface for the [SO3LR](https://github.com/general-molecular-simulations/so3lr) calculator"""

from .ase import ASEDriver
from ipi.utils.messages import verbosity, warning

So3lrCalculator = None

__DRIVER_NAME__ = "so3lr"
__DRIVER_CLASS__ = "SO3LR_driver"


class SO3LR_driver(ASEDriver):
    """
    SO3LR driver requires a template file that describes the chemical makeup of the structure.

    Optionally, lr_cutoff and dispersion_energy_lr_cutoff_damping can be specified.

    Example: python driver.py -m so3lr -u -o template.xyz,lr_cutoff=12,dispersion_energy_lr_cutoff_damping=2
    """

    def __init__(self, *args, **kwargs):
        warning(
            "THIS PES HAS NOT BEEN TESTED FOLLOWING CONVERSION TO THE NEW PES API.",
            verbosity.low,
        )
        global So3lrCalculator
        try:
            from so3lr import So3lrCalculator
        except:
            raise ImportError("Couldn't load so3lr bindings")

        super().__init__(*args, **kwargs)

    def check_parameters(self):
        super().check_parameters()

        self.ase_calculator = So3lrCalculator(**self.kwargs)

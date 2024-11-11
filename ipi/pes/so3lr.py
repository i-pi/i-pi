""" An interface for the [SO3LR](https://github.com/general-molecular-simulations/so3lr) calculator """

from .ase import ASEDriver

So3lrCalculator = None

__DRIVER_NAME__ = "so3lr"
__DRIVER_CLASS__ = "SO3LR_driver"

ERROR_MSG = """
SO3LR driver requires a template file that describes the chemical makeup of the structure.

Optionally, lr_cutoff and dispersion_energy_lr_cutoff_damping can be specified.

Example: python driver.py -m so3lr -u -o template.xyz,lr_cutoff=12,dispersion_energy_lr_cutoff_damping=2
"""


class SO3LR_driver(ASEDriver):
    def __init__(self, args=None, verbose=False):

        global So3lrCalculator
        try:
            from so3lr import So3lrCalculator
        except:
            raise ImportError("Couldn't load so3lr bindings")

        super().__init__(args, verbose, ERROR_MSG)

    def check_arguments(self):
        super().check_arguments()

        args = self.args

        kwargs = {}
        if len(args) >= 2:
            _ = args[0]  # template we don't need

            for arg in args[1:]:
                key, value = arg.split("=")
                kwargs[key] = float(value)

        self.ase_calculator = So3lrCalculator(**kwargs)

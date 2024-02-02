# this is how the driver will be referred to in the input files
__DRIVER_NAME__ = "dummy"
__DRIVER_CLASS__ = "Dummy_driver"


class Dummy_driver(object):
    """Base class providing the structure of a PES for the python driver."""

    def __init__(
        self, args="", verbose=False, error_msg="Invalid arguments for the PES"
    ):
        """Initialized dummy drivers"""
        self.error_msg = error_msg
        self.args = args.split(",")
        self.verbose = verbose
        self.check_arguments()

    def check_arguments(self):
        """Dummy function that checks the arguments required to run the driver"""
        pass

    def __call__(self, cell, pos):
        """Does nothing, but returns properties that can be used by the driver loop."""
        pot = 0.0
        force = pos * 0.0  # makes a zero force with same shape as pos
        vir = cell * 0.0  # makes a zero virial with same shape as cell
        extras = "nada"
        return pot, force, vir, extras

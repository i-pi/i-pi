# this is how the driver will be referred to in the input files
__DRIVER_NAME__ = "dummy"
__DRIVER_CLASS__ = "Dummy_driver"

import json


class Dummy_driver(object):
    """Base class providing the structure of a PES for the python driver.

    Command line:
        i-pi-py_driver -m dummy [...]
    Init arguments:
        :param verbose: bool to determine whether the PES should output verbose info.
    """

    def __init__(self, verbose=False, *args, **kwargs):
        """Initialized dummy drivers"""
        self.verbose = verbose
        self.args = args
        self.kwargs = kwargs

        self.check_parameters()

    def check_parameters(self):
        """Dummy function that checks the arguments required to run the driver"""
        pass

    def __call__(self, cell, pos):
        """Does nothing, but returns properties that can be used by the driver loop."""
        pot = 0.0
        force = pos * 0.0  # makes a zero force with same shape as pos
        vir = cell * 0.0  # makes a zero virial with same shape as cell
        extras = json.dumps(
            {"dipole": [0.0, 0.0, 0.0]}
        )  # have json formatting to potentially work with some test examples. meaningless value
        return pot, force, vir, extras

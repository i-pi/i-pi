"""Interface with [elphmod](https://github.com/janberges/elphmod) MD driver."""

import sys
from .dummy import Dummy_driver

__DRIVER_NAME__ = "elphmod"
__DRIVER_CLASS__ = "ModelIIIDriver"


class ModelIIIDriver(Dummy_driver):
    """Wrapper around elphmod MD driver.
    A pickled driver instance is required (see elphmod.md.Driver.save).

    Example: python3 driver.py -u -m elphmod -o driver.pickle
    """

    def check_parameters(self):
        """Check arguments and load driver instance."""

        import elphmod

        if len(self.args) != 1:
            sys.exit(self.__doc__)

        self.driver = elphmod.md.Driver.load(self.args[0])

    def __call__(self, cell, pos):
        """Calculate energy and forces for given structure."""

        return self.driver(cell, pos)

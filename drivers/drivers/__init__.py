""" Small functions/classes providing access to driver PES to be called from driver.py """

from .harmonic import Harmonic_driver
from .spline import Spline_driver
from .dummy import Dummy_driver

__all__ = ["__drivers__", "Dummy_driver", "Harmonic_driver", "Spline_driver"]

# dictionary linking strings
__drivers__ = {
    "dummy": Dummy_driver,
    "harmonic": Harmonic_driver,
    "spline": Spline_driver,
}

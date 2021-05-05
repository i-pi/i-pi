""" Small functions/classes providing access to driver PES to be called from driver.py """

from .dummy import Dummy_driver
from .harmonic import Harmonic_driver
from .spline import Spline_driver
from .doublewell import DoubleWell_driver, DoubleWell_with_friction_driver
from .rascal import Rascal_driver

__all__ = [
    "__drivers__",
    "Dummy_driver",
    "Harmonic_driver",
    "Rascal_driver",
    "Spline_driver",
    "DoubleWell_driver",
    "DoubleWell_with_friction_driver",
]

# dictionary linking strings
__drivers__ = {
    "dummy": Dummy_driver,
    "harmonic": Harmonic_driver,
    "rascal": Rascal_driver,
    "DoubleWell": DoubleWell_driver,
    "DW_friction": DoubleWell_with_friction_driver,
    "spline": Spline_driver,
}

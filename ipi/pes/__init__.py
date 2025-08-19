"""Small functions/classes providing access to driver PES to be called from driver.py"""

import importlib
from .dummy import Dummy_driver

__drivers__ = {
    "ase": "ase",
    "bath": "bath",
    "double_double_well": "doubledouble_well",
    "DW": "doublewell",
    "DW_bath": "doublewell_with_bath",
    "DW_friction": "doublewell_with_friction",
    "driverdipole": "driverdipole",
    "dummy": "dummy",
    "elphmod": "elphmod",
    "harmonic": "harmonic",
    "mace": "mace",
    "metatensor": "metatensor",
    "metatomic": "metatomic",
    "MorseHarmonic": "morse",
    "pet": "pet",
    "psiflow": "psiflow",
    "rascal": "rascal",
    "so3lr": "so3lr",
    "Spherical_LJ": "spherical_LJ",
    "spline": "spline",
    "xtb": "xtb",
}


# sort alphabetically
__drivers__ = dict(sorted(__drivers__.items()))


def load_driver(module_name: str) -> Dummy_driver:
    """
    Dynamically load a driver class from a module.

    Imports the given module, looks up its `__DRIVER_CLASS__` attribute to
    determine the class name, and returns the corresponding class object.

    Args:
        module_name: Name of the module to import.

    Returns:
        The driver class defined in the module.

    Raises:
        AttributeError: If the module does not define `__DRIVER_CLASS__`.
    """

    module = importlib.import_module(module_name, __package__)

    driver_class_name = getattr(module, "__DRIVER_CLASS__", None)
    if driver_class_name is None:
        raise AttributeError(
            f"Module '{module_name}' does not define '__DRIVER_CLASS__'."
        )

    driver_class = getattr(module, driver_class_name)
    return driver_class

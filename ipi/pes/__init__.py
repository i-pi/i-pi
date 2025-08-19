"""Small functions/classes providing access to driver PES to be called from driver.py"""

import os, sys
import importlib
from importlib import util
from pathlib import Path
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


def load_driver(mode: str) -> Dummy_driver:
    """
    Load a driver class for the given mode.

    Looks for a local <mode>.py file or a module in `ipi.pes`.
    The module must define `__DRIVER_CLASS__`.

    Parameters
    ----------
    mode : str
        Driver name or local Python file (without `.py`).

    Returns
    -------
    Dummy_driver
        The driver class.

    Raises
    ------
    ValueError, AttributeError
    """

    if os.path.isfile(f"{mode}.py"):
        # import client from <args.mode>.py (in the current working directory)
        file_path = Path(f"{mode}.py")
        spec = util.spec_from_file_location(file_path.stem, file_path)
        module = util.module_from_spec(spec)
        sys.modules[file_path.stem] = module
        spec.loader.exec_module(module)
    else:
        # import client from ipi/pes/<module_name>.py
        if mode not in __drivers__:
            choices = ", ".join(__drivers__.keys())
            raise ValueError(f"Invalid mode '{mode}'. Available modes: {choices}")
        module_name = f"ipi.pes.{__drivers__[mode]}"
        module = importlib.import_module(module_name, __package__)

    driver_class_name = getattr(module, "__DRIVER_CLASS__", None)
    if driver_class_name is None:
        raise AttributeError(
            f"Module '{module_name}' does not define '__DRIVER_CLASS__'."
        )

    driver_class = getattr(module, driver_class_name)
    return driver_class

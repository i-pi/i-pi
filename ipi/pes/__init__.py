"""Small functions/classes providing access to driver PES to be called from driver.py"""

import sys
import importlib
from importlib import util
from pathlib import Path
from .dummy import Dummy_driver

__drivers__ = {
    "ase": "ase",
    "bath": "bath",
    "custom": None,
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


def load_driver(mode: str, file_path: str) -> Dummy_driver:
    """
    Load a driver class by mode.

    If `mode="custom"`, imports the driver from `file_path`.
    Otherwise, loads a registered driver from `ipi.pes`.
    The module must define `__DRIVER_CLASS__`.

    Parameters
    ----------
    mode : str
        "custom" or a registered driver name.
    file_path : str
        Path to a `.py` file (used only if mode="custom").

    Returns
    -------
    Dummy_driver
        The driver class.

    Raises
    ------
    ValueError
        If mode is invalid.
    AssertionError
        If the driver is unimplemented.
    AttributeError
        If `__DRIVER_CLASS__` is missing.
    """

    if mode == "custom":
        # import client from <file_path>
        file_path = Path(file_path)
        spec = util.spec_from_file_location(file_path.stem, file_path)
        module = util.module_from_spec(spec)
        sys.modules[file_path.stem] = module
        spec.loader.exec_module(module)
    else:
        # import client from ipi/pes/<module_name>.py
        if mode not in __drivers__:
            choices = ", ".join(__drivers__.keys())
            raise ValueError(f"Invalid mode '{mode}'. Available modes: {choices}")
        file_path = __drivers__[mode]
        assert file_path is not None, f"Driver '{mode}' is not implemented."
        module_name = f"ipi.pes.{file_path}"
        module = importlib.import_module(module_name, __package__)

    driver_class_name = getattr(module, "__DRIVER_CLASS__", None)
    if driver_class_name is None:
        raise AttributeError(
            f"Module '{module_name}' does not define '__DRIVER_CLASS__'."
        )

    driver_class = getattr(module, driver_class_name)
    return driver_class

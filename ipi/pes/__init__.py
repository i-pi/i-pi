"""Small functions/classes providing access to driver PES to be called from driver.py"""

import os
import sys
import ast
import importlib
from importlib import util
from pathlib import Path
from .dummy import Dummy_driver
from ipi.utils.messages import warning
import pkgutil


def scan_pes_file(pyfile: str):
    """
    Return (__DRIVER_CLASS__, __DRIVER_NAME__) if both are string literals, else (None, None).

    Parse a Python source file without 'running' the file and extract values of __DRIVER_CLASS__ and __DRIVER_NAME__
    if they are assigned as string literals at the top level of the module.
    """
    try:
        with open(pyfile, "r", encoding="utf-8") as f:
            tree = ast.parse(f.read(), filename=pyfile)
    except OSError as e:
        # raise an error if the file cannot be opened/read
        # this will prevent from adding badly formatted files to this folder
        raise ValueError(f"Could not open or read file '{pyfile}'\n{e}")

    driver_class = driver_name = None
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id in {
                    "__DRIVER_CLASS__",
                    "__DRIVER_NAME__",
                }:
                    try:
                        value = ast.literal_eval(node.value)
                    except Exception:
                        value = None
                    if target.id == "__DRIVER_CLASS__" and isinstance(value, str):
                        driver_class = value
                    elif target.id == "__DRIVER_NAME__" and isinstance(value, str):
                        driver_name = value
    return driver_class, driver_name


def load_pes(mode: str, pes_path: str) -> Dummy_driver:
    """
    Load a driver PES by mode.

    If `mode="custom"`, imports the driver from `pes_path`.
    Otherwise, loads a registered driver from `ipi.pes`.
    The module must define `__DRIVER_CLASS__`.

    Parameters
    ----------
    mode : str
        "custom" or a registered driver name.
    pes_path : str
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
        # import client from <pes_path>
        pes_path = Path(pes_path)
        spec = util.spec_from_file_location(pes_path.stem, pes_path)
        module = util.module_from_spec(spec)
        sys.modules[pes_path.stem] = module
        spec.loader.exec_module(module)
    else:
        # import client from ipi/pes/<module_name>.py
        if mode not in __drivers__:
            choices = ", ".join(__drivers__.keys())
            raise ValueError(f"Invalid mode '{mode}'. Available modes: {choices}")
        pes_path = __drivers__[mode]
        assert pes_path is not None, f"Driver '{mode}' is not implemented."
        module_name = f"ipi.pes.{pes_path}"
        try:
            module = importlib.import_module(module_name, __package__)
        except:
            print(f"Could not import module '{module_name}'")
            raise

    driver_class_name = getattr(module, "__DRIVER_CLASS__", None)
    if driver_class_name is None:
        raise AttributeError(
            f"Module '{module_name}' does not define '__DRIVER_CLASS__'."
        )

    driver_class = getattr(module, driver_class_name)
    return driver_class


__all__ = []

# Dictionary to store driver name to class mapping
__drivers__ = {}

# Iterate through all modules in the current package folder to detect
# available PES drivers
for loader, module_name, is_pkg in pkgutil.iter_modules(__path__):
    spec = importlib.util.find_spec(f"{__package__}.{module_name}")
    if not spec or not spec.origin or not spec.origin.endswith(".py"):
        continue

    if os.path.basename(spec.origin) == "tools.py":
        continue  # skip private modules

    driver_class, driver_name = scan_pes_file(spec.origin)
    if not (driver_class and driver_name):
        # let's raise a warning for future developers and be sure that we are not adding badly formatted files
        warning(
            f"Module '{module_name}' does not define both __DRIVER_CLASS__ and __DRIVER_NAME__ as string literals. Skipping."
        )
        continue  # not a permissible driver

    # If both class and name are defined, update __drivers__,
    # but don't load the module just yet
    if driver_class and driver_name:
        __drivers__[driver_name] = module_name

# sort alphabetically
__drivers__ = dict(sorted(__drivers__.items()))
__all__.append("__drivers__")

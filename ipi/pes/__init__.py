"""Small functions/classes providing access to driver PES to be called from driver.py"""

import importlib
import traceback
from .dummy import Dummy_driver

__drivers__ = [
    "ase",
    "bath",
    "double_double_well",
    "DW",
    "DW_bath",
    "DW_friction",
    "driverdipole",
    "dummy",
    "elphmod",
    "harmonic",
    "mace",
    "metatensor",
    "metatomic",
    "MorseHarmonic",
    "pet",
    "psiflow",
    "rascal",
    "so3lr",
    "Spherical_LJ",
    "spline",
    "xtb",
]

# sort alphabetically
__drivers__ = sorted(__drivers__)


def load_driver(mode: str):
    """
    Loads a PES driver module and retrieves its driver class.

    Parameters:
        mode (str): The name of the PES mode/module to load.

    Returns:
        driver_class: The loaded driver class object.

    Raises:
        ImportError: If the module or driver class cannot be loaded.
    """
    module_name = f"ipi.pes.{mode}"

    try:
        module = importlib.import_module(module_name, __package__)

        driver_class_name = getattr(module, "__DRIVER_CLASS__", None)
        if driver_class_name is None:
            raise AttributeError(
                f"Module '{module_name}' does not define '__DRIVER_CLASS__'."
            )

        driver_class = getattr(module, driver_class_name)
        return driver_class

    except ModuleNotFoundError as e:
        print(f"\n[ERROR] PES module '{module_name}' could not be found.")
        traceback.print_exc()
        raise ImportError(
            f"Could not import PES module '{module_name}'. "
            f"Please check the mode argument: '{mode}'."
        ) from e

    except AttributeError as e:
        print(f"\n[ERROR] Driver class not found in module '{module_name}'.")
        traceback.print_exc()
        raise ImportError(
            f"'{module_name}' does not define the expected class: {e}"
        ) from e

    except Exception as e:
        print(
            f"\n[ERROR] Unexpected error while importing or accessing driver class from '{module_name}'."
        )
        traceback.print_exc()
        raise ImportError(
            f"An unexpected error occurred while loading PES module '{module_name}': {e}"
        ) from e

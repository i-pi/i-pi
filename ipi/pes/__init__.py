"""Small functions/classes providing access to driver PES to be called from driver.py"""

import pkgutil
import importlib
import traceback

__all__ = []

# Dictionary to store driver name to class mapping
__drivers__ = {}

# Iterate through all modules in the current package folder
for loader, module_name, is_pkg in pkgutil.iter_modules(__path__):
    # Import the module
    try:
        module = importlib.import_module("." + module_name, __package__)
    except Exception:
        print(f"!! Could not import PES module {module_name} !!")
        traceback.print_exc()

    # Get the driver class and name from the module
    driver_class = getattr(module, "__DRIVER_CLASS__", None)
    driver_name = getattr(module, "__DRIVER_NAME__", None)

    # If both class and name are defined, update __all__ and __drivers__
    if driver_class and driver_name:
        __all__.append(driver_class)
        __drivers__[driver_name] = getattr(module, driver_class)
        globals()[driver_class] = getattr(module, driver_class)  # add class to globals
    else:
        if driver_class != "driver_tools":
            raise ImportError(
                f"PES module `{module_name}` does not define __DRIVER_CLASS__ and __DRIVER_NAME__"
            )

__all__.append("__drivers__")

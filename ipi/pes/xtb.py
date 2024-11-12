import numpy as np
import json

from ipi.utils.units import unit_to_internal, unit_to_user
from ipi.utils.messages import verbosity, warning
from .dummy import Dummy_driver

tb = None

__DRIVER_NAME__ = "xtb"
__DRIVER_CLASS__ = "TBLiteDriver"


class TBLiteDriver(Dummy_driver):
    """
    Interface with tblite to provide GFN1-xTB and GFN2-xTB calculators.

    Example::

        python driver.py -m xtb -u -o config.json
    """

    def __init__(self, json_input, *args, **kwargs):
        warning(
            "THIS PES HAS NOT BEEN TESTED FOLLOWING CONVERSION TO THE NEW PES API.",
            verbosity.low,
        )
        config = json.load(open(json_input))

        global tb
        try:
            import tblite.interface as tb
        except ImportError:
            raise ModuleNotFoundError(
                "Could not find tblite for xtb driver. Please install tblite-python with mamba"
            )

        try:
            self.method = config["method"]
            self.numbers = np.asarray(config["numbers"])
        except KeyError as e:
            raise ValueError(f"Required key {str(e)} not found.")
        self.charge = config.get("charge")
        self.uhf = config.get("uhf")
        self.periodic = config.get("periodic")
        self.verbosity = 1 if self.verbose else 0

        super().__init__(*args, **kwargs)

    def __call__(self, cell, pos):
        """
        Get energies, forces, and stresses from the tblite library
        This routine assumes that the client will take positions
        in bohr, and return energies in hartree, and forces
        in hartree/abohr.
        """

        pos = unit_to_user("length", "atomic_unit", pos)
        cell = unit_to_user("length", "atomic_unit", cell.T)

        calc = tb.Calculator(
            self.method,
            self.numbers,
            np.asarray(pos),
            self.charge,
            self.uhf,  # unpaired electrons
            np.asarray(cell) if self.periodic else None,
            np.repeat(self.periodic, 3) if self.periodic else None,
        )
        calc.set("verbosity", self.verbosity)

        # Do the actual calculation
        results = calc.singlepoint()
        pot = results.get("energy")
        grad = results.get("gradient")
        vir = results.get("virial")

        # converts to internal quantities
        pot_ipi = np.asarray(unit_to_internal("energy", "atomic_unit", pot), np.float64)
        force_ipi = np.asarray(
            unit_to_internal("force", "atomic_unit", -grad), np.float64
        )
        vir_ipi = np.array(
            unit_to_internal("energy", "atomic_unit", -vir.T), np.float64
        )
        extras = ""

        return pot_ipi, force_ipi, vir_ipi, extras

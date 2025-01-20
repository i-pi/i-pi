import numpy as np

from ipi.utils.units import unit_to_internal, unit_to_user
from ipi.utils.io import read_file_name
from .dummy import Dummy_driver

tb = None

__DRIVER_NAME__ = "xtb"
__DRIVER_CLASS__ = "TBLiteDriver"


class TBLiteDriver(Dummy_driver):
    """
    Interface with tblite to provide GFN1-xTB and GFN2-xTB calculators.

    Example::

        i-pi-py_driver -m xtb -u -o template=input.xyz,method=GFN2-xTB
    """

    def __init__(
        self, template, method, charge=None, uhf=None, periodic=None, *args, **kwargs
    ):

        super().__init__(*args, **kwargs)

        global tb
        try:
            import tblite.interface as tb
        except ImportError:
            raise ModuleNotFoundError(
                "Could not find tblite for xtb driver. Please install tblite-python with mamba"
            )

        input_data = read_file_name(template)
        atoms = input_data["atoms"]
        symbols = atoms.names[:]
        numbers = np.asarray(tb.symbols_to_numbers(symbols))
        self.periodic = periodic
        self.verbosity = 1 if self.verbose else 0

        pos = unit_to_user("length", "atomic_unit", atoms.q[:])
        cell = unit_to_user("length", "atomic_unit", input_data["cell"].h.T)

        self.calc = tb.Calculator(
            method,
            numbers,
            np.asarray(pos),
            charge,
            uhf,  # unpaired electrons
            np.asarray(cell) if self.periodic else None,
            np.repeat(self.periodic, 3) if self.periodic else None,
        )
        self.calc.set("verbosity", self.verbosity)

    def __call__(self, cell, pos):
        """
        Get energies, forces, and stresses from the tblite library
        This routine assumes that the client will take positions
        in bohr, and return energies in hartree, and forces
        in hartree/abohr.
        """

        pos = unit_to_user("length", "atomic_unit", pos)
        cell = unit_to_user("length", "atomic_unit", cell.T)

        self.calc.update(np.asarray(pos), np.asarray(cell) if self.periodic else None)

        # Do the actual calculation
        results = self.calc.singlepoint()
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

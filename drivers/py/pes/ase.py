"""Interface with ASE calculators"""

import sys
import numpy as np
from .dummy import Dummy_driver

from ipi.utils.units import unit_to_internal, unit_to_user
from ipi.utils.messages import warning

try:
    from ase.io import read
except ImportError:
    warning("Could not find or import the ASE module")
    MetatensorCalculator = None

__DRIVER_NAME__ = "ase"
__DRIVER_CLASS__ = "ASEDriver"

ERROR_MSG = """
This ASE driver requires specification of and ASE calculator
and an ASE-readable template file that describes the chemical makeup of the structure.

Example: python driver.py -m ase -u -o template.xyz,model_parameters
"""


class ASEDriver(Dummy_driver):
    """Abstract base class using an arbitrary ASE calculator as i-pi driver"""

    def __init__(self, args=None, verbose=False, error_msg=ERROR_MSG):
        super().__init__(args, verbose, error_msg=error_msg)

    def check_arguments(self):
        """Check the arguments required to run the driver

        This loads the potential and atoms template in metatensor
        """

        if len(self.args) >= 1:
            self.template = self.args[0]
        else:
            sys.exit(self.error_msg)

        self.template_ase = read(self.template)

        self.ase_calculator = None

    def __call__(self, cell, pos):
        """Get energies, forces, and stresses from the ASE calculator
        This routine assumes that the client will take positions
        in angstrom, and return energies in electronvolt, and forces
        in ev/ang.
        """

        # ASE calculators assume angstrom and eV units
        pos = unit_to_user("length", "angstrom", pos)
        # ASE expects cell-vectors-as-rows
        cell = unit_to_user("length", "angstrom", cell.T)
        # applies the cell and positions to the template
        structure = self.template_ase.copy()
        structure.positions[:] = pos
        structure.cell[:] = cell
        structure.calc = self.ase_calculator

        # Do the actual calculation
        properties = structure.get_properties(["energy", "forces", "stress"])
        pot = properties["energy"]
        force = properties["forces"]
        stress = properties["stress"]
        if len(stress) == 6:
            # converts from voight notation
            stress = np.array(stress[[0, 5, 4, 5, 1, 3, 4, 3, 2]])

        # converts to internal quantities
        pot_ipi = np.asarray(
            unit_to_internal("energy", "electronvolt", pot), np.float64
        )
        force_ipi = np.asarray(unit_to_internal("force", "ev/ang", force), np.float64)
        vir_calc = -stress * structure.get_volume()
        vir_ipi = np.array(
            unit_to_internal("energy", "electronvolt", vir_calc.T), dtype=np.float64
        )
        extras = ""

        return pot_ipi, force_ipi, vir_ipi, extras

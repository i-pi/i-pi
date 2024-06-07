"""Interface with [QMC]"""

import sys
import numpy as np
from .dummy import Dummy_driver

from ipi.utils.units import unit_to_internal, unit_to_user
from ipi.utils.messages import warning

try:
    from QMC import SingleStructCalculator as QMCCalc

    try:
        from ase.io import read
    except ImportError:
        warning("The QMC driver has an ASE dependency")
        raise
except ImportError:
    warning("Could not find or import the QMC module")
    QMCCalc = None

__DRIVER_NAME__ = "QMC"
__DRIVER_CLASS__ = "QMC_driver"


class QMC_driver(Dummy_driver):
    def __init__(self, args=None, verbose=False):
        self.error_msg = """
TODO
"""

        super().__init__(args, verbose)

        if QMCCalc is None:
            raise ImportError("Couldn't load QMC bindings")

    def check_arguments(self):
        """ """
        args = self.args

        pass

        self.template_ase = read(self.template)
        self.template_ase.arrays["forces"] = np.zeros_like(self.template_ase.positions)
        self.QMC_calc = QMCCalc(
            self.model_path,
            **kwargs,
        )

    def __call__(self, cell, pos):
        """Get energies, forces, and stresses from the QMC model
        This routine assumes that the client will take positions
        in angstrom, and return energies in electronvolt, and forces
        in ev/ang.
        """

        pos_QMC = unit_to_user("length", "angstrom", pos)
        # QMC expects ASE-format, cell-vectors-as-rows
        cell_QMC = unit_to_user("length", "angstrom", cell.T)
        # applies the cell and positions to the template
        QMC_structure = self.template_ase.copy()
        QMC_structure.positions = pos_QMC
        QMC_structure.cell = cell_QMC

        # Do the actual calculation
        pot, force = self.QMC_calc.forward(QMC_structure)
        pot_ipi = pot
        force_ipi = np.asarray(unit_to_internal("force", "ev/ang", force), np.float64)
        # QMC does not yet compute stress
        vir_QMC = 0 * np.eye(3)
        vir_ipi = unit_to_internal("energy", "electronvolt", vir_QMC.T)
        extras = ""
        return pot_ipi, force_ipi, vir_ipi, extras

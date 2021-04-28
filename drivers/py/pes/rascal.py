"""Interface with librascal to run machine learning potentials"""

import sys
from .dummy import Dummy_driver

from ipi.utils.mathtools import det_ut3x3
from ipi.utils.units import unit_to_internal, unit_to_user

try:
    from rascal.models.genericmd import GenericMDCalculator as RascalCalc
except:
    RascalCalc = None


class Rascal_driver(Dummy_driver):
    def __init__(self, args=None):

        self.error_msg = """Rascal driver requires specification of a .json model file fitted with librascal, 
                            and a template file that describes the chemical makeup of the structure. 
                            Example: python driver.py -m rascal -u -o model.json,template.xyz"""

        super().__init__(args)

        if RascalCalc is None:
            raise ImportError("Couldn't load librascal bindings")

    def check_arguments(self):
        """Check the arguments required to run the driver

        This loads the potential and atoms template in librascal
        """
        try:
            arglist = self.args.split(",")
        except ValueError:
            sys.exit(self.error_msg)

        if len(arglist) == 2:
            self.model = arglist[0]
            self.template = arglist[1]
        else:
            sys.exit(self.error_msg)

        self.rascal_calc = RascalCalc(self.model, True, self.template)

    def __call__(self, cell, pos):
        """Get energies, forces, and stresses from the librascal model"""
        pos_rascal = unit_to_user("length", "angstrom", pos)
        cell_rascal = unit_to_user("length", "angstrom", cell)
        # Do the actual calculation
        pot, force, stress = self.rascal_calc.calculate(pos_rascal, cell_rascal)
        pot_ipi = unit_to_internal("energy", "electronvolt", pot)
        force_ipi = unit_to_internal("force", "ev/ang", force)
        # The rascal stress is normalized by the cell volume (in rascal units)
        vir_rascal = -1 * stress * det_ut3x3(cell_rascal)
        vir_ipi = unit_to_internal("energy", "electronvolt", vir_rascal)
        extras = ""
        return pot_ipi, force_ipi, vir_ipi, extras

"""Interface with librascal to run machine learning potentials"""

import sys, os
import numpy as np
from .dummy import Dummy_driver

from ipi.utils.units import unit_to_internal, unit_to_user

try:
    sys.path.append(os.getcwd() + "/pet/src")
    from single_struct_calculator import SingleStructCalculator as PETCalc

    try:
        from ase.io import read
    except:
        raise ImportError("The PET driver has an ASE dependency")
except:
    PETCalc = None

__DRIVER_NAME__ = "pet"
__DRIVER_CLASS__ = "PET_driver"


class PET_driver(Dummy_driver):
    def __init__(self, args=None):
        self.error_msg = """
The PET driver requires specification of a .json model file fitted with librascal, 
and a template file that describes the chemical makeup of the structure. 

Example: python driver.py -m pet -u -o model.json,template.xyz
"""

        super().__init__(args)

        if PETCalc is None:
            raise ImportError("Couldn't load PET bindings")

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

        self.template_ase = read(self.template)
        self.template_ase.arrays["forces"] = np.zeros_like(self.template_ase.positions)
        self.pet_calc = PETCalc(self.model)

    def __call__(self, cell, pos):
        """Get energies, forces, and stresses from the PET model"""
        pos_pet = unit_to_user("length", "angstrom", pos)
        # librascal expects ASE-format, cell-vectors-as-rows
        cell_pet = unit_to_user("length", "angstrom", cell.T)
        # applies the cell and positions to the template
        pet_structure = self.template_ase.copy()
        pet_structure.positions = pos_pet
        pet_structure.cell = cell_pet

        # Do the actual calculation
        pot, force = self.pet_calc.forward(pet_structure)
        pot_ipi = np.asarray(
            unit_to_internal("energy", "electronvolt", pot), np.float64
        )
        force_ipi = np.asarray(unit_to_internal("force", "ev/ang", force), np.float64)
        # PET does not yet compute stress
        vir_pet = 0 * np.eye(3)
        vir_ipi = unit_to_internal("energy", "electronvolt", vir_pet.T)
        extras = ""
        return pot_ipi, force_ipi, vir_ipi, extras

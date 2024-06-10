"""Interface with [PET](https://github.com/serfg/pet) models to run machine learning potentials"""

import sys
import numpy as np
from .dummy import Dummy_driver

from ipi.utils.units import unit_to_internal, unit_to_user
from ipi.utils.messages import warning

try:
    from pet import SingleStructCalculator as PETCalc

    try:
        from ase.io import read
    except ImportError:
        warning("The PET driver has an ASE dependency")
        raise
except ImportError:
    warning("Could not find or import the PET module")
    PETCalc = None

__DRIVER_NAME__ = "pet"
__DRIVER_CLASS__ = "PET_driver"


class PET_driver(Dummy_driver):
    def __init__(self, args=None, verbose=False):
        self.error_msg = """
The PET driver requires (a) a path to the results/experiment_name folder emitted by pet_train
                        (b) a path to an ase.io.read-able file with a prototype structure

Other arguments to the pet.SingleStructCalculator class can be optionally
supplied in key=value form after the required arguments.

Example: python driver.py -m pet -u -o "path/to/results/name,template.xyz,device=cuda"
"""

        super().__init__(args, verbose)

        if PETCalc is None:
            raise ImportError("Couldn't load PET bindings")

    def check_arguments(self):
        """Check the arguments required to run the driver

        This loads the potential and atoms template in PET
        """
        args = self.args

        if len(args) >= 2:
            self.model_path = args[0]
            self.template = args[1]
            kwargs = {}
            if len(args) > 2:
                for arg in args[2:]:
                    key, value = arg.split("=")
                    lookup = {
                        "None": None,
                        "False": False,
                        "True": True,
                    }

                    if value in lookup:
                        value = lookup[value]

                    kwargs[key] = value

        else:
            sys.exit(self.error_msg)

        self.template_ase = read(self.template)
        self.template_ase.arrays["forces"] = np.zeros_like(self.template_ase.positions)
        self.pet_calc = PETCalc(
            self.model_path,
            **kwargs,
        )

    def __call__(self, cell, pos):
        """Get energies, forces, and stresses from the PET model
        This routine assumes that the client will take positions
        in angstrom, and return energies in electronvolt, and forces
        in ev/ang.
        """

        pos_pet = unit_to_user("length", "angstrom", pos)
        # PET expects ASE-format, cell-vectors-as-rows
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

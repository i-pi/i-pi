"""Interface with [PET](https://github.com/serfg/pet) models to run machine learning potentials"""

import numpy as np
from .dummy import Dummy_driver

from ipi.utils.units import unit_to_internal, unit_to_user
from ipi.utils.messages import warning

PETCalc = None
read = None

__DRIVER_NAME__ = "pet"
__DRIVER_CLASS__ = "PET_driver"


class PET_driver(Dummy_driver):
    _error_msg = """
The PET driver requires (a) a path to the results/experiment_name folder emitted by pet_train
                        (b) a path to an ase.io.read-able file with a prototype structure

Other arguments to the pet.SingleStructCalculator class can be optionally
supplied in key=value form after the required arguments.

Example: python driver.py -m pet -u -o "path/to/results/name,template.xyz,device=cuda"
"""

    def __init__(self, model_path, template, *args, **kwargs):
        global PETCalc, read

        try:
            from pet import SingleStructCalculator as PETCalc

            try:
                from ase.io import read
            except ImportError:
                warning("The PET driver has an ASE dependency")
                raise
        except ImportError:
            raise ImportError("Could not find or import the PET module")

        self.model_path = model_path
        self.template = template
        super().__init__(*args, **kwargs)

    def check_parameters(self):
        """Check the arguments required to run the driver

        This loads the potential and atoms template in PET
        """

        self.template_ase = read(self.template)
        self.template_ase.arrays["forces"] = np.zeros_like(self.template_ase.positions)
        self.pet_calc = PETCalc(
            self.model_path,
            **self.kwargs,
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

import sys
import numpy as np
from ipi.utils.units import unit_to_internal, unit_to_user
from .dummy import Dummy_driver

from ase.io import read

try:
    from mace.calculators import MACECalculator
except:
    MACECalculator = None

__DRIVER_NAME__ = "mace"
__DRIVER_CLASS__ = "MACE_driver"


class MACE_driver(Dummy_driver):
    def __init__(self, args=None, verbose=False):
        self.error_msg = """MACE driver requires specification of a .json model,
                            and a template file that describes the chemical makeup of the structure.
                            Example: python driver.py -m mace -u -o model.json,template.xyz"""

        super().__init__(args, verbose)

    def check_arguments(self):
        """Check the arguments requuired to run the driver

        This loads the potential and atoms template in MACE
        """
        try:
            arglist = self.args.split(",")
        except ValueError:
            sys.exit(self.error_msg)

        self.model_atoms = read(arglist[0])
        self.driver_example_atoms = arglist[0]
        self.driver_model_path = arglist[1]

    def __call__(self, cell, pos):
        """Get energies, forces, and stresses from the MACE model
        This routine assumes that the client will take positions
        in angstrom, and return energies in electronvolt, and forces
        in ev/ang.
        """
        pos_calc = unit_to_user("length", "angstrom", pos)
        cell_calc = unit_to_user("length", "angstrom", cell.T)

        atoms = read(self.driver_example_atoms)
        atoms.set_pbc([True, True, True])
        atoms.set_cell(cell_calc, scale_atoms=True)
        atoms.set_positions(pos_calc)

        mace_calculator = MACECalculator(
            model_path=self.driver_model_path, device="cpu"
        )
        atoms.set_calculator(mace_calculator)

        pot = atoms.get_potential_energy()
        force = atoms.get_forces()
        stress = atoms.get_stress(voigt=False)

        pot_ipi = float(unit_to_internal("energy", "electronvolt", pot))
        force_ipi = np.array(
            unit_to_internal("force", "ev/ang", force.reshape(-1, 3)), dtype=np.float64
        )

        vir_calc = -stress * self.model_atoms.get_volume()
        vir_ipi = np.array(
            unit_to_internal("energy", "electronvolt", vir_calc.T), dtype=np.float64
        )
        extras = ""

        return pot_ipi, force_ipi, vir_ipi, extras

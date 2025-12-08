"""Interface with ASE calculators"""

import json
import numpy as np
from .dummy import Dummy_driver

from ipi.utils.units import unit_to_internal, unit_to_user

from ase.io import read
from ase import Atoms


__DRIVER_NAME__ = "ase"
__DRIVER_CLASS__ = "ASEDriver"


class ASEDriver(Dummy_driver):
    """
    Base class using an arbitrary ASE calculator as i-pi driver.
    Should not be called directly as it does not set a calculator.

    Parameters:
        :param verbose: bool, whether to print verbose output
        :param template: string, ASE-readable filename where to get the structure data from
        :param has_energy: bool, whether the model can compute energy
        :param has_forces: bool, whether the model can compute forces
        :param has_stress: bool, whether the model can compute stress
    """

    def __init__(
        self,
        template,
        has_energy=True,
        has_forces=True,
        has_stress=True,
        *args,
        **kwargs
    ):

        self.template = template
        self.capabilities = []
        if has_energy:
            self.capabilities.append("energy")
        if has_forces:
            self.capabilities.append("forces")
        if has_stress:
            self.capabilities.append("stress")
        super().__init__(*args, **kwargs)

    def check_parameters(self):
        """Check the arguments required to run the driver

        This loads the potential and atoms template for ASE
        """

        self.template_ase = read(self.template)
        self.ase_calculator = None

    def convert_units(self, cell: np.ndarray, pos: np.ndarray):
        # ASE calculators assume angstrom and eV units
        pos = unit_to_user("length", "angstrom", pos)
        # ASE expects cell-vectors-as-rows
        cell = unit_to_user("length", "angstrom", cell.T)
        # applies the cell and positions to the template
        return cell, pos

    def compute_structure(self, cell: np.ndarray, pos: np.ndarray):
        """Get energies, forces, and stresses from the ASE calculator
        This routine assumes that the client will take positions
        in angstrom, and return energies in electronvolt, and forces
        in ev/ang.
        """
        cell, pos = self.convert_units(cell, pos)

        structure = self.template_ase.copy()
        structure.positions[:] = pos
        structure.cell[:] = cell
        structure.calc = self.ase_calculator

        # Do the actual calculation
        properties = structure.get_properties(self.capabilities)

        return self.post_process(properties, structure)

    def post_process(self, properties, structure: Atoms):

        pot = properties["energy"] if "energy" in self.capabilities else 0.0
        force = (
            properties["forces"]
            if "forces" in self.capabilities
            else np.zeros_like(structure.positions)
        )
        stress = properties["stress"] if "stress" in self.capabilities else np.zeros(9)
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

        # Extra information apart from the "mandatory" ones (energy, forces, and stress).
        # The model could return the dipole, polarizability, Born charges for example.
        # These information will be stored in the 'extras' variable (a python dict)
        # with the same keys as returned by the model and values converted to float or list.
        extras = {}
        for key in properties:
            if key not in ["energy", "forces", "stress"]:
                value = properties[key]
                if isinstance(value, np.ndarray):
                    extras[key] = value.tolist()
                else:
                    extras[key] = float(value)
        if extras == {}:
            extras = ""
        else:
            extras = json.dumps(extras)

        return pot_ipi, force_ipi, vir_ipi, extras

"""Interface with ASE calculators"""

import numpy as np
from .dummy import Dummy_driver

from ipi.utils.units import unit_to_internal, unit_to_user
from ipi.utils.messages import warning


read = None

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
        global read
        try:
            from ase.io import read
        except ImportError:
            warning("Could not find or import the ASE module")

        global all_changes
        from ase.calculators.calculator import all_changes

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

        This loads the potential and atoms template in metatensor
        """

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
        properties = structure.get_properties(self.capabilities)

        pot = properties["energy"] if "energy" in self.capabilities else 0.0
        force = (
            properties["forces"]
            if "forces" in self.capabilities
            else np.zeros_like(pos)
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
        extras = ""

        return pot_ipi, force_ipi, vir_ipi, extras

"""Spherical Lennard-Jones potential driver.

This module implements a spherical Lennard-Jones (LJ) potential centered around a given point
and applies it to a specified set of atomic species. It supports both NumPy and PyTorch backends,
although NumPy is used by default due to similar performance.

The driver is designed to be used with i-PI and ASE frameworks, functioning both as an ASE
calculator and an i-PI driver.
"""

import json
import numpy as np
from ipi.utils.units import unit_to_internal, unit_to_user
from ipi.utils.messages import warning

# Attempt relative import first (for i-PI structure), fall back to absolute
try:
    from .ase import ASEDriver
except:
    from ase import ASEDriver

from ase.calculators.calculator import Calculator
from ase import Atoms

# Compatibility for different ASE versions
try:
    from ase.calculators.calculator import all_changes
except ImportError:
    from ase.calculators import all_changes

__DRIVER_NAME__ = "Spherical_LJ"
__DRIVER_CLASS__ = "Spherical_LJ_driver"

# ---------------------- #
# Attention/comments:
# - the docstring of the class explains how to use the driven and which parameters to provide.
# - a similar confined potential has been used in https://doi.org/10.1103/PhysRevResearch.1.033145
# - pay attention that the center of mass of your systems should roughly be at the center of the sphere (see "center" parameter),
#       however, if any atoms are outside the sphere, the code will raise a warning (see 'warning_message' variable).


# ---------------------- #
class Spherical_LJ_driver(ASEDriver, Calculator):
    """
    Spherical Lennard-Jones potential driver.
    Driver implementing a Spherical Lennard-Jones (LJ) potential.
    Parameters must be passed via a dictionary or a JSON file with the following keys:
    {
        "center": [0,0,0],              // center of the sphere
        "radius": 10,                   // radius of the spherical potential
        "sigma": 2.569,                 // sigma parameter of the LJ potential
        "epsilon": 2.754,               // epsilon parameter of the LJ potential

        "center_unit": "angstrom",      // unit of the center (supported: angstrom, atomic_unit[bohr] + derived quantities)
        "radius_unit": "angstrom",      // unit of the radius (supported: angstrom, atomic_unit[bohr] + derived quantities)
        "sigma_unit": "angstrom",       // unit of sigma      (supported: angstrom, atomic_unit[bohr] + derived quantities)
        "epsilon_unit": "kilocal/mol",  // unit of epsilon    (supported: electronvolt, atomic_unit[Hartree], j/mol, cal/mol, kelvin + derived quantities)

        "symbols": ["O"],               // list of chemical symbols to whom the potential applies
        "first_power": 9,               // first power of the LJ potential
        "second_power": 3               // second power of the LJ potential
    }
    Attention: no default parameters are provided, so you must specify all of them.
    You can just copy and paste the above example and modify the values as needed.

    This class is both a i-PI driver and an ASE calculator:
        - 'Spherical_LJ_driver.calculate' will compute the potential energy and forces (in eV and eV/ang respectively).
        - 'ASEDriver.__call__' will take care of converting the results into atomic units.
    """

    dimensions = {
        "center": "length",
        "radius": "length",
        "sigma": "length",
        "epsilon": "energy",
    }

    units = {"length": "angstrom", "energy": "electronvolt"}

    def __init__(
        self,
        template,
        instructions,
        has_energy=True,
        has_forces=True,
        has_stress=False,
        *args,
        **kwargs,
    ):

        assert (
            not has_stress
        ), "Spherical_LJ_driver does not support stress calculation."

        Calculator.__init__(self, *args, **kwargs)
        ASEDriver.__init__(
            self, template, has_energy, has_forces, has_stress, *args, **kwargs
        )

        # Read instructions from file or use provided dictionary
        if isinstance(instructions, str):
            with open(instructions, "r") as f:
                instructions = json.load(f)
        elif isinstance(instructions, dict):
            pass
        else:
            raise ValueError("`instructions` can be `str` or `dict` only.")

        # convert parameters to the required units
        to_delete = list()
        for k, dimension in self.dimensions.items():
            variable = f"{k}_unit"
            if variable in instructions:
                to_delete.append(variable)
                if instructions[variable] is not None:
                    factor = convert(
                        1,
                        dimension,
                        _from=instructions[variable],
                        _to=self.units[dimension],
                    )
                    instructions[k] = process_input(instructions[k]) * factor
        for k in to_delete:
            del instructions[k]

        self.instructions = instructions

    def check_parameters(self):
        """Check the arguments required to run the driver"""
        super().check_parameters()
        self.ase_calculator = self

    def calculate(self, atoms: Atoms, properties=None, system_changes=all_changes):
        """
        Core method that calculates energy and forces for given atoms using
        the spherical Lennard-Jones potential.
        """
        super().calculate(atoms, properties, system_changes)

        # atomic positions [ang]
        pos = atoms.get_positions()

        # atoms to be considered
        indices = [
            n
            for n, atom in enumerate(atoms)
            if atom.symbol in self.instructions["symbols"]
        ]

        # ---------------------- #
        potential, forces = compute_energy_and_forces(
            pos[indices], self.instructions, False
        )

        # ---------------------- #
        all_forces = np.zeros((atoms.get_global_number_of_atoms(), 3))
        all_forces[indices] = forces

        self.results = {
            "free_energy": potential,
            "energy": potential,
            "forces": all_forces,
            "stress": np.zeros(6),  # No stress calculation
        }
        return


# ---------------------- #
# Here starts a list of functions that are used to compute
# the potential and forces withing the Spherical_LJ_driver class.


# ---------------------- #
def compute_energy_and_forces(pos, instructions):
    """
    Computes potential and analytical forces using NumPy.

    Returns:
        potential (float), forces (ndarray)
    """
    center = np.asarray(instructions["center"])
    radius = float(instructions["radius"])
    sigma = float(instructions["sigma"])
    epsilon = float(instructions["epsilon"])
    first_power = int(instructions["first_power"])
    second_power = int(instructions["second_power"])

    delta = pos - center
    r = radius - np.linalg.norm(delta, axis=1)

    # Check if any atom is outside the spherical potential
    if np.any(r <= 0):
        warning(
            "Some atoms are outside the spherical potential. This can lead to numerical instability."
        )

    # Calculate the potential
    sr = sigma / r
    potential = epsilon * (2.0 / 15.0 * sr**first_power - sr**second_power)
    potential = np.sum(potential)

    # Calculate the forces
    r_vec = np.linalg.norm(delta, axis=1)  # ||delta||

    # Correct derivative of the potential including the 2/15 factor
    dU_dr = -epsilon * (
        -first_power * (2.0 / 15.0) * (sigma**first_power) / (r ** (first_power + 1))
        + second_power * (sigma**second_power) / (r ** (second_power + 1))
    )

    unit_vectors = delta / r_vec[:, np.newaxis]  # normalize delta
    forces = -dU_dr[:, np.newaxis] * unit_vectors

    return potential, forces


# ---------------------- #
def convert(
    what: float,
    family: str = None,
    _from: str = "atomic_unit",
    _to: str = "atomic_unit",
) -> float:
    """
    Converts a physical quantity between units of the same type (length, energy, etc.)
    Example:
    value = convert(7.6,'length','angstrom','atomic_unit')
    arr = convert([1,3,4],'energy','atomic_unit','millielectronvolt')
    """
    # from ipi.utils.units import unit_to_internal, unit_to_user
    if family is not None:
        factor = unit_to_internal(family, _from, 1)
        factor *= unit_to_user(family, _to, 1)
        return what * factor
    else:
        return what


# ---------------------- #
def process_input(value):
    """
    Standardizes user input into numerical format (float or np.array).
    """
    if isinstance(value, float):
        return value
    elif isinstance(value, int):
        return value
    elif isinstance(value, list):
        return np.array(value)
    else:
        raise TypeError("Input must be a float or a list.")

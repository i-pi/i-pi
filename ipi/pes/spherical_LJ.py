"""Spherical Lennard-Jones potential"""

import json
import numpy as np
from ipi.utils.units import unit_to_internal, unit_to_user
from ipi.utils.messages import warning

try:
    from .ase import ASEDriver
except:
    from ase import ASEDriver

from ase.calculators.calculator import Calculator
from ase import Atoms

try:
    from ase.calculators.calculator import all_changes
except ImportError:
    from ase.calculators import all_changes

__DRIVER_NAME__ = "Spherical_LJ"
__DRIVER_CLASS__ = "Spherical_LJ_driver"

USE_TORCH = True
DEBUG = True
warning_message = "Some atoms are outside the spherical potential. This can lead to numerical instability."

if USE_TORCH:
    import torch

    torch.set_default_dtype(torch.float64)

def convert(
    what: float,
    family: str = None,
    _from: str = "atomic_unit",
    _to: str = "atomic_unit",
) -> float:
    """Convert a quantity from one unit of a specific family to another.
    Example:
    value = convert(7.6,'length','angstrom','atomic_unit')
    arr = convert([1,3,4],'energy','atomic_unit','millielectronvolt')"""
    # from ipi.utils.units import unit_to_internal, unit_to_user
    if family is not None:
        factor = unit_to_internal(family, _from, 1)
        factor *= unit_to_user(family, _to, 1)
        return what * factor
    else:
        return what


def process_input(value):
    """Process the input value to ensure it is in the correct format."""
    if isinstance(value, float):
        return value
    elif isinstance(value, int):
        return value
    elif isinstance(value, list):
        return np.array(value)
    else:
        raise TypeError("Input must be a float or a list.")


class Spherical_LJ_driver(ASEDriver, Calculator):
    """
    Spherical Lennard-Jones potential driver.

    You need to provide a JSON file or a dict structured as follows:
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
        super().calculate(atoms, properties, system_changes)

        # with torch.autograd.set_detect_anomaly(False):

        pos = atoms.get_positions()

        # Attention, do not use any inplace operations with torch, as it can lead to issues with gradients
        indices = [
            n
            for n, atom in enumerate(atoms)
            if atom.symbol in self.instructions["symbols"]
        ]
        pos_to_use_np = pos[indices]

        #----------------------#    
        if USE_TORCH or DEBUG:
            potential, forces = calculate_torch(
                pos_to_use_np, self.instructions,True
            )
            
            # compute the potential and forces with numpy for debugging
            if DEBUG:
                potential_np, forces_np = calculate_numpy(
                pos_to_use_np, self.instructions,False
            )
            # Check that the potential and forces calculated with torch and numpy are the same
            assert np.allclose(
                potential, potential_np
            ), "Potential mismatch between torch and numpy calculations."
            assert np.allclose(
                forces, forces_np
            ), "Forces mismatch between torch and numpy calculations."
            
        #----------------------#
        else:
            potential, forces = calculate_numpy(
                pos_to_use_np, self.instructions,False
            )

        #----------------------#
        all_forces = np.zeros((atoms.get_global_number_of_atoms(), 3))
        all_forces[indices] = forces

        self.results = {
            "free_energy": potential,
            "energy": potential,
            "forces": all_forces,
            "stress": np.zeros(6),  # No stress calculation
        }
        return
    
def calculate_torch(pos_to_use_np,instructions,use_torch):
    
    center = np.asarray(instructions["center"])
    radius = float(instructions["radius"])
    sigma = float(instructions["sigma"])
    epsilon = float(instructions["epsilon"])
    first_power = int(instructions["first_power"])
    second_power = int(instructions["second_power"])
        
    pos_leaf = torch.tensor(pos_to_use_np, requires_grad=True)
    pos_leaf.retain_grad()
    center_torch = torch.tensor(center)

    # The shifting happens inside LJ_potential / LJ_forces
    potential, _, _ = LJ_potential(
        pos_leaf,
        center_torch,
        radius,
        sigma,
        epsilon,
        first_power,
        second_power,
        use_torch,
    )
    # forces = LJ_forces(
    #     pos_leaf,
    #     potential,
    #     center_torch,
    #     None,
    #     None,
    #     None,
    #     None,
    #     None,
    #     use_torch,
    # )
    
    potential.backward()
    if pos_leaf.grad is None:
        raise RuntimeError(
            "No gradient found. Ensure pos_leaf is a leaf tensor with requires_grad=True."
        )
    forces = -pos_leaf.grad

    potential = float(potential.detach().cpu().numpy())
    forces = forces.detach().cpu().numpy()
    
    return potential, forces

def calculate_numpy(pos_to_use_np, instructions,use_torch):
    
    center = np.asarray(instructions["center"])
    radius = float(instructions["radius"])
    sigma = float(instructions["sigma"])
    epsilon = float(instructions["epsilon"])
    first_power = int(instructions["first_power"])
    second_power = int(instructions["second_power"])
    
    pos_to_use = pos_to_use_np
    potential, r, delta  = LJ_potential(
        pos_to_use,
        center,
        radius,
        sigma,
        epsilon,
        first_power,
        second_power,
        use_torch,
    )
    # forces = LJ_forces(
    #     pos_to_use,
    #     None,
    #     center,
    #     radius,
    #     sigma,
    #     epsilon,
    #     first_power,
    #     second_power,
    #     use_torch,
    # )
    
    # delta = pos_to_use - center
    r_vec = np.linalg.norm(delta, axis=1)       # ||delta||
    # r = radius - r_vec                          # same as in potential
    if np.any(r <= 0):
        warning(warning_message)

    # Correct derivative of the potential including the 2/15 factor
    dU_dr = epsilon * (
        -first_power * (2.0 / 15.0) * (sigma ** first_power) / (r ** (first_power + 1)) +
        second_power * (sigma ** second_power) / (r ** (second_power + 1))
    )

    unit_vectors = delta / r_vec[:, np.newaxis]  # normalize delta
    forces = dU_dr[:, np.newaxis] * unit_vectors
    
    return potential, forces

def to_numpy(x):
    """Convert input to numpy array if it is a torch tensor."""
    if isinstance(x, torch.Tensor):
        return x.detach().cpu().numpy()
    return x

def potential_expression(sr,epsilon,first_power,second_power):
    return epsilon * (2.0 / 15.0 * sr**first_power - sr**second_power)

def LJ_potential(pos, center, radius, sigma, epsilon, first_power, second_power, use_torch):

    delta = pos - center
    if use_torch:
        r = radius - torch.norm(delta, dim=1)
    else:
        r = radius - np.linalg.norm(delta, axis=1)
        
    if np.any(to_numpy(r) <= 0):
        warning(warning_message)

    sr = sigma / r
    potential = potential_expression(sr,epsilon,first_power,second_power)
    
    if use_torch:
        return torch.sum(potential), r, delta
    else:
        return np.sum(potential), r, delta

def LJ_forces(
    pos_leaf,
    potential,
    center,
    radius,
    sigma,
    epsilon,
    first_power,
    second_power,
    use_torch,
):
    if use_torch:
        potential.backward()
        if pos_leaf.grad is None:
            raise RuntimeError(
                "No gradient found. Ensure pos_leaf is a leaf tensor with requires_grad=True."
            )
        return -pos_leaf.grad
    else:
        delta = pos_leaf - center
        r_vec = np.linalg.norm(delta, axis=1)       # ||delta||
        r = radius - r_vec                          # same as in potential
        if np.any(r <= 0):
            warning(warning_message)

        # Correct derivative of the potential including the 2/15 factor
        dU_dr = epsilon * (
            -first_power * (2.0 / 15.0) * (sigma ** first_power) / (r ** (first_power + 1)) +
            second_power * (sigma ** second_power) / (r ** (second_power + 1))
        )

        unit_vectors = delta / r_vec[:, np.newaxis]  # normalize delta
        return dU_dr[:, np.newaxis] * unit_vectors

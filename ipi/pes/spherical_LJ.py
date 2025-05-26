"""Spherical Lennard-Jones potential"""

import json

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

import numpy as np
from ipi.utils.units import unit_to_internal, unit_to_user

__DRIVER_NAME__ = "Spherical_LJ"
__DRIVER_CLASS__ = "Spherical_LJ_driver"

USE_TORCH = True

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

    dimensions = {"center": "length", "radius": "length", "sigma": "length", "epsilon": "energy"}

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
        """Check the arguments required to run the driver
        """
        super().check_parameters()
        self.ase_calculator = self

    def calculate(self, atoms:Atoms, properties=None, system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        
        with torch.autograd.set_detect_anomaly(USE_TORCH):
        
            pos = atoms.get_positions()
            # if USE_TORCH:
            #     center = torch.tensor(self.instructions["center"])
            #     radius = torch.tensor(self.instructions["radius"])
            #     sigma = torch.tensor(self.instructions["sigma"])
            #     epsilon = torch.tensor(self.instructions["epsilon"])
            #     first_power = torch.tensor(self.instructions["first_power"],dtype=torch.int32)
            #     second_power = torch.tensor(self.instructions["second_power"],dtype=torch.int32)
                
            # else:
            center = np.asarray(self.instructions["center"])
            radius = float(self.instructions["radius"])
            sigma = float(self.instructions["sigma"])
            epsilon = float(self.instructions["epsilon"])
            first_power = int(self.instructions["first_power"])
            second_power = int(self.instructions["second_power"])
                
            delta = pos - center
            indices = [n  for n,atom in enumerate(atoms) if atom.symbol in self.instructions["symbols"] ]
            pos_to_use = delta[indices]
            
            if USE_TORCH:
                pos_to_use = torch.from_numpy(pos_to_use).clone().detach().requires_grad_()
                pos_to_use.retain_grad()
        
            # potential, forces = self.potential_and_forces(pos_to_use) 
            
            # here inside torch will be used
            potential = LJ_potential(pos_to_use,radius,sigma,epsilon,first_power,second_power,USE_TORCH)        
            forces = LJ_forces(pos_to_use,potential,radius,sigma,epsilon,first_power,second_power,USE_TORCH)
            
            if USE_TORCH:
                potential =  float(potential.detach().cpu().numpy())
                forces = forces.detach().cpu().numpy()
            
            # ----------------------# 
            
            all_forces = np.zeros((atoms.get_global_number_of_atoms(), 3))
            all_forces[indices] = forces
            
            self.results = {
                "free_energy": potential,
                "energy": potential,
                "forces": all_forces,
                "stress": np.zeros(6),  # No stress calculation
            }
        return
    
    # def potential_and_forces(self, pos:np.ndarray):
        
    #     if USE_TORCH:
            
    #         pos = torch.tensor(pos,requires_grad=True)
    #         radius = torch.tensor(self.instructions["radius"])
    #         sigma = torch.tensor(self.instructions["sigma"])
    #         epsilon = torch.tensor(self.instructions["epsilon"])
    #         first_power = torch.tensor(self.instructions["first_power"],dtype=torch.int32)
    #         second_power = torch.tensor(self.instructions["second_power"],dtype=torch.int32)
    #     else:
    #         pos = np.array(pos, dtype=np.float64)
    #         radius = float(self.instructions["radius"])
    #         sigma = float(self.instructions["sigma"])
    #         epsilon = float(self.instructions["epsilon"])
    #         first_power = int(self.instructions["first_power"])
    #         second_power = int(self.instructions["second_power"])
            
        
    #     potential = LJ_potential(pos,radius,sigma,epsilon,first_power,second_power,USE_TORCH)        
    #     forces = LJ_forces(pos,potential,radius,sigma,epsilon,first_power,second_power,USE_TORCH)
        
    #     if USE_TORCH:
    #         return float(potential.detach().cpu().numpy()), forces.detach().cpu().numpy()
    #     else:
    #         return potential, forces
    
def LJ_potential(pos: np.ndarray,radius:float,sigma:float,epsilon:float, first_power:int, second_power:int,use_torch:bool=USE_TORCH) -> np.ndarray:
    assert pos.ndim == 2 and pos.shape[1] == 3, "Positions must be a 2D array with shape (n_atoms, 3)"
    if use_torch:
        r = torch.norm(pos, dim=1)
    else:
        r:np.ndarray = np.linalg.norm(pos, axis=1)
    # Attention: don't use inplace operations with torch, as it can lead to issues with gradients
    # for example, 'r -= radius' will not work
    r = r - radius # distance from the wall of the sphere
    assert r.ndim == 1, "Input must be a 1D array of distances."
    potential = 4 * epsilon * ((sigma / r) ** first_power - (sigma / r) ** second_power)
    if use_torch:
        return torch.sum(potential)
    else:
        return np.sum(potential)

def LJ_forces(pos: np.ndarray,potential:float, radius: float, sigma: float, epsilon: float, first_power: int, second_power: int,use_torch:bool=USE_TORCH) -> np.ndarray:
    assert pos.ndim == 2 and pos.shape[1] == 3, "Positions must be a 2D array with shape (n_atoms, 3)"
    if use_torch:
        assert pos.requires_grad, "Positions must require gradient for force calculation."
        potential.backward()  # Compute gradients
        potential.retain_grad()
        forces = -pos.grad  # Forces are the negative gradient of the potential
        return forces
    
    else:        
        r = np.linalg.norm(pos, axis=1)  # shape (n_atoms,)
        d = r - radius  # distance from the wall

        # Compute scalar part of the derivative
        coeff = 4 * epsilon
        dU_dr = coeff * (
            -first_power * (sigma ** first_power) / (d ** (first_power + 1)) +
            second_power * (sigma ** second_power) / (d ** (second_power + 1))
        )

        # Unit vectors: pos / r
        unit_vectors = pos / r[:, np.newaxis]

        forces = -dU_dr[:, np.newaxis] * unit_vectors  # shape (n_atoms, 3)
        
        return forces
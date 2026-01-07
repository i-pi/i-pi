"""=Lennard-Jones potential driver.

This module implements a spherical Lennard-Jones (LJ) potential centered around a given point
and applies it to a specified set of atomic species. It supports both NumPy and PyTorch backends,
although NumPy is used by default due to similar performance.

The driver is designed to be used with i-PI and ASE frameworks, functioning both as an ASE
calculator and an i-PI driver.
"""

import json
import numpy as np
from ipi.pes.tools import convert, process_input
from ipi.utils.messages import warning
from ipi.pes.dummy import Dummy_driver

# ---------------------- #
# Attention/comments:
# - the docstring of the class explains how to use the driven and which parameters to provide.
# - a similar confined potential has been used in https://doi.org/10.1103/PhysRevResearch.1.033145
# - pay attention that the center of mass of your systems should roughly be at the center of the sphere (see "center" parameter),
#       however, if any atoms are outside the sphere, the code will raise a warning (see 'warning_message' variable).


# ---------------------- #
class LJdriver(Dummy_driver):
    """
    Spherical Lennard-Jones potential driver.
    Driver implementing a Spherical Lennard-Jones (LJ) potential.
    Parameters must be passed via a dictionary or a JSON file with the following keys:
    {
        "z_plane": 10,                  // z-coordinate of the plane
        "sigma": 2.569,                 // sigma parameter of the LJ potential
        "epsilon": 2.754,               // epsilon parameter of the LJ potential

        "z_plane_unit": "angstrom",     // unit of the z-coordinate (supported: angstrom, atomic_unit[bohr] + derived quantities)
        "sigma_unit": "angstrom",       // unit of sigma      (supported: angstrom, atomic_unit[bohr] + derived quantities)
        "epsilon_unit": "kilocal/mol",  // unit of epsilon    (supported: electronvolt, atomic_unit[Hartree], j/mol, cal/mol, kelvin + derived quantities)

        "symbols": ["O"],               // list of chemical symbols to whom the potential applies
        "first_power": 9,               // first power of the LJ potential
        "second_power": 3               // second power of the LJ potential
    }
    Attention: no default parameters are provided, so you must specify all of them.
    You can just copy and paste the above example and modify the values as needed.
    """

    dimensions = {
        "sigma": "length",
        "epsilon": "energy",
    }

    units = {"length": "atomic_unit", "energy": "atomic_unit"}

    def __init__(
        self,
        template: str,
        instructions: dict,
        symbols: list = None,
        has_energy=True,
        has_forces=True,
        has_stress=False,
        *args,
        **kwargs,
    ):
        assert has_energy, "LJWall_driver requires energy calculation."
        assert has_forces, "LJWall_driver requires forces calculation."
        assert (
            not has_stress
        ), "LJWall_driver does not support stress calculation."

        super().__init__(self, *args, **kwargs)

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

        # Initialize symbols
        self.atoms = None
        if symbols is None:
            # Let's try to use ASE ...
            try:
                from ase.io import read

                atoms = read(template)
                self.symbols = atoms.get_chemical_symbols()
                assert atoms.positions.shape == (
                    len(self.symbols),
                    3,
                ), "Positions shape mismatch."
            except ImportError:
                warning("Could not find or import the ASE module")
        else:
            # ... but the user can also provide symbols directly
            self.symbols = symbols

    def compute(self, cell:np.ndarray, pos:np.ndarray):
        """
        Core method that calculates energy and forces for given atoms using
        the spherical Lennard-Jones potential.
        """
        
        if not isinstance(cell, list):
            return self.compute([cell], [pos])
        
        assert len(cell) == len(pos), "Cell and position lists must have the same length."
        assert cell.ndim == 3 and cell.shape[1:] == (3, 3), "Cell must be a list of 3x3 matrices."
        assert pos.ndim == 3 and pos.shape[1:] == (len(self.symbols), 3), "Position must be a list of Nx3 matrices."

        potential, forces, vir, extras = super().compute(cell, pos)

        # atoms to be considered
        indices = [
            n
            for n, symbol in enumerate(self.symbols)
            if symbol in self.instructions["symbols"]
        ]

        potential, some_forces = self.compute_energy_and_forces(
            pos[indices]
        )

        forces[indices, :] = some_forces

        return potential, forces, vir, extras
    
    def compute_energy_and_forces(self,pos:np.ndarray)->tuple:
        raise ValueError("This method should have been overridden.")
    
    def lennard_jones(self,r:np.ndarray):
        """
        Computes potential and analytical forces using NumPy.

        Returns:
            potential (float), forces (ndarray)
        """
        
        assert r.ndim == 2, "Input positions must be a 2D array of distances."
        
        sigma = float(self.instructions["sigma"])
        epsilon = float(self.instructions["epsilon"])
        first_power = int(self.instructions["first_power"])
        second_power = int(self.instructions["second_power"])

        # Check if any atom is outside the spherical potential
        if np.any(r <= 0):
            warning(
                "Some atoms are outside the spherical potential. This can lead to numerical instability."
            )

        # Calculate the potential
        sr = sigma / r
        potential = epsilon * ( sr**first_power - sr**second_power)
        potential = np.sum(potential)

        # Correct derivative of the potential including the 2/15 factor
        forces = epsilon * (
            -first_power * (sigma**first_power) / (r ** (first_power + 1))
            + second_power * (sigma**second_power) / (r ** (second_power + 1))
        )
        
        assert forces.shape == r.shape, "Forces shape mismatch."

        return potential, forces

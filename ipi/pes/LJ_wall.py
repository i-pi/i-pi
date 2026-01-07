"""Spherical Lennard-Jones potential driver.

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
from ipi.pes.LennardJones import LJdriver

# Attempt relative import first (for i-PI structure), fall back to absolute
try:
    from .dummy import Dummy_driver
except:
    from dummy import Dummy_driver


__DRIVER_NAME__ = "LJWall"
__DRIVER_CLASS__ = "LJWall_driver"

# ---------------------- #
# Attention/comments:
# - the docstring of the class explains how to use the driven and which parameters to provide.
# - a similar confined potential has been used in https://doi.org/10.1103/PhysRevResearch.1.033145
# - pay attention that the center of mass of your systems should roughly be at the center of the sphere (see "center" parameter),
#       however, if any atoms are outside the sphere, the code will raise a warning (see 'warning_message' variable).


# ---------------------- #
class LJWall_driver(LJdriver):
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
        **LJdriver.dimensions,
        **{
            "z_plane": "length",
        },
    }
    
    def compute_energy_and_forces(self,pos:np.ndarray)->tuple: 
        

    



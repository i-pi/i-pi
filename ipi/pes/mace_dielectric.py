""" An interface for the [MACE](https://github.com/ACEsuit/mace) calculator to get dielectric properties"""

import os
import json
import numpy as np
from .mace import MACE_driver
from ipi.utils.messages import warning
from ipi.utils.units import unit_to_internal

MACECalculator = None

__DRIVER_NAME__ = "mace_dielectric"
__DRIVER_CLASS__ = "MACE_dielectric_driver"

METHODS_AVAILABLE = [
    "fixed-D",  # https://doi.org/10.1103/PhysRevB.93.144201
    "fixed-E",  # https://doi.org/10.48550/arXiv.2502.02413
    "none",  # just for debugging purposes
]

# ToDo:
# - check the units of the electric displacement
# - write fixed_D
# - add flag to chose whether to :
#   -- compute the energy and forces, dipole and BEC, and then compute the sum of the forces
#   -- compute the energy, add E x mu, and compute the total forces directly (easier and faster)


class MACE_dielectric_driver(MACE_driver):

    method: str
    instructions: str

    def __init__(self, method: str, instructions: str, *args, **kwargs):
        """Driver for the MACE MLIPs to compute dielectric properties."""

        super().__init__(*args, **kwargs)

        method = str(method)
        instructions = str(instructions)

        if method not in METHODS_AVAILABLE:
            raise ValueError(f"method must be one of {METHODS_AVAILABLE}")
        if not os.path.exists(instructions):
            raise ValueError(f"File '{instructions}' does not exist.")
        if not str(instructions).endswith(".json"):
            raise ValueError(f"File '{instructions}' must be a json file.")

        self.method = method

        with open(instructions, "r") as f:
            tmp_data = json.load(f)

        # ------------------------------ #
        # fixed electric field (E)
        if self.method == "fixed-E":

            # 'instructions' must contain the following
            # {
            #     "E": 0.4,          // mandatory
            #     "E-unit": "V/ang", // optional (default: atomic_unit)
            # }

            if "E" not in tmp_data:
                raise ValueError(
                    f"File '{instructions}' must contain 'E' for method 'fixed-E'."
                )
            if "E-unit" not in tmp_data:
                warning(
                    f"File '{instructions}' does not contain 'E-unit': we will be assumed that 'E' is in 'atomic_unit'."
                )
                tmp_data["E-unit"] = "atomic_unit"

            E = unit_to_internal("electric-field", tmp_data["E-unit"], tmp_data["E"])
            self.instructions = {
                "E": E,
            }

        # ------------------------------ #
        # fixed electric displacement (D)
        elif self.method == "fixed-D":

            # 'instructions' must contain the following
            # {
            #     "D": 4,           // mandatory
            #     "D-unit": "eang", // optional (default: atomic_unit)
            #     "V": 123,         // mandatory
            #     "V-unit": "ang3"  // optional (default: atomic_unit)
            # }

            if "D" not in self.instructions:
                raise ValueError(
                    f"File '{instructions}' must contain 'D' (electric displacement) for method 'fixed-D'."
                )
            if "D-unit" not in self.instructions:
                warning(
                    f"File '{instructions}' does not contain 'D-unit' (electric displacement unit): we will assumed that 'D' is in 'atomic_unit'."
                )
                tmp_data["D-unit"] = "atomic_unit"
            if "V" not in self.instructions:
                raise ValueError(
                    f"File '{instructions}' must contain 'V' (volume) for method 'fixed-D'."
                )
            if "V-unit" not in self.instructions:
                warning(
                    f"File '{instructions}' does not contain 'volume' (volume unit): we will assumed that 'volume' is in 'atomic_unit'."
                )
                tmp_data["V-unit"] = "atomic_unit"

            D = unit_to_internal("electric-field", tmp_data["D-unit"], tmp_data["D"])
            V = unit_to_internal("volume", tmp_data["V-unit"], tmp_data["V"])
            self.instructions = {
                "D": D,
                "V": V,
            }

    def __call__(self, cell, pos):

        pot_ipi, force_ipi, vir_ipi, extras_ipi = super().__call__(cell, pos)

        # just returns what provided by the MACE model
        if self.method == "none":
            return pot_ipi, force_ipi, vir_ipi, extras_ipi

        # add the contribution to the fixed electric field (E) to the potential and forces
        elif self.method == "fixed-E":
            pot, force, vir, extras = fixed_E(extras_ipi, self.instructions["E"])

        # add the contribution to the fixed electric displacement (D) to the potential and forces
        elif self.method == "fixed-D":
            pot, force, vir, extras = fixed_D(
                extras_ipi, self.instructions["D"], self.instructions["V"]
            )

        else:  # should never reach here
            raise ValueError(
                f"This is an implementation error. 'method' = {self.method} but it should be one of {METHODS_AVAILABLE}."
            )

        pot_ipi += pot
        force_ipi += force
        vir_ipi += vir
        extras_ipi = {**extras_ipi, **extras}

        return pot_ipi, force_ipi, vir_ipi, extras_ipi


# ------------------------------ #
def fixed_E(extras: dict, Efield: np.ndarray):
    pot = -np.asarray(extras["dipole"]) @ Efield
    force = extras["BEC"] @ Efield
    vir = np.zeros((3, 3))
    return pot, force, vir, extras


# ------------------------------ #
def fixed_D(extras: dict, displacement: np.ndarray, volume: float):
    polarization = extras["dipole"] / volume
    Efield = displacement - 4 * np.pi * polarization
    return fixed_E(extras, Efield)

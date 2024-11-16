from typing import List

import numpy as np

from ipi.utils.units import unit_to_internal, unit_to_user
from ipi.utils.io import read_file_name
from .dummy import Dummy_driver

tb = None

__DRIVER_NAME__ = "xtb"
__DRIVER_CLASS__ = "TBLiteDriver"


class TBLiteDriver(Dummy_driver):
    """
    Interface with tblite to provide GFN1-xTB and GFN2-xTB calculators.

    Example::

        i-pi-py_driver -m xtb -u -o template=input.xyz,method=GFN2-xTB
    """

    def __init__(
        self, template, method, charge=None, uhf=None, periodic=None, *args, **kwargs
    ):

        super().__init__(*args, **kwargs)

        global tb
        try:
            import tblite.interface as tb
        except ImportError:
            raise ModuleNotFoundError(
                "Could not find tblite for xtb driver. Please install tblite-python with mamba"
            )

        input_data = read_file_name(template)
        atoms = input_data["atoms"]
        symbols = atoms.names[:]
        numbers = np.asarray(symbols_to_numbers(symbols))
        self.periodic = periodic
        self.verbosity = 1 if self.verbose else 0

        pos = unit_to_user("length", "atomic_unit", atoms.q[:])
        cell = unit_to_user("length", "atomic_unit", input_data["cell"].h.T)

        self.calc = tb.Calculator(
            method,
            numbers,
            np.asarray(pos),
            charge,
            uhf,  # unpaired electrons
            np.asarray(cell) if self.periodic else None,
            np.repeat(self.periodic, 3) if self.periodic else None,
        )
        self.calc.set("verbosity", self.verbosity)

    def __call__(self, cell, pos):
        """
        Get energies, forces, and stresses from the tblite library
        This routine assumes that the client will take positions
        in bohr, and return energies in hartree, and forces
        in hartree/abohr.
        """

        pos = unit_to_user("length", "atomic_unit", pos)
        cell = unit_to_user("length", "atomic_unit", cell.T)

        self.calc.update(np.asarray(pos), np.asarray(cell) if self.periodic else None)

        # Do the actual calculation
        results = self.calc.singlepoint()
        pot = results.get("energy")
        grad = results.get("gradient")
        vir = results.get("virial")

        # converts to internal quantities
        pot_ipi = np.asarray(unit_to_internal("energy", "atomic_unit", pot), np.float64)
        force_ipi = np.asarray(
            unit_to_internal("force", "atomic_unit", -grad), np.float64
        )
        vir_ipi = np.array(
            unit_to_internal("energy", "atomic_unit", -vir.T), np.float64
        )
        extras = ""

        return pot_ipi, force_ipi, vir_ipi, extras


def symbols_to_numbers(symbols: List[str]) -> List[int]:
    return [symbol_to_number(symbol) for symbol in symbols]


def symbol_to_number(symbol: str) -> int:
    return {
        "H": 1,
        "He": 2,
        "Li": 3,
        "Be": 4,
        "B": 5,
        "C": 6,
        "N": 7,
        "O": 8,
        "F": 9,
        "Ne": 10,
        "Na": 11,
        "Mg": 12,
        "Al": 13,
        "Si": 14,
        "P": 15,
        "S": 16,
        "Cl": 17,
        "Ar": 18,
        "K": 19,
        "Ca": 20,
        "Sc": 21,
        "Ti": 22,
        "V": 23,
        "Cr": 24,
        "Mn": 25,
        "Fe": 26,
        "Co": 27,
        "Ni": 28,
        "Cu": 29,
        "Zn": 30,
        "Ga": 31,
        "Ge": 32,
        "As": 33,
        "Se": 34,
        "Br": 35,
        "Kr": 36,
        "Rb": 37,
        "Sr": 38,
        "Y": 39,
        "Zr": 40,
        "Nb": 41,
        "Mo": 42,
        "Tc": 43,
        "Ru": 44,
        "Rh": 45,
        "Pd": 46,
        "Ag": 47,
        "Cd": 48,
        "In": 49,
        "Sn": 50,
        "Sb": 51,
        "Te": 52,
        "I": 53,
        "Xe": 54,
        "Cs": 55,
        "Ba": 56,
        "La": 57,
        "Ce": 58,
        "Pr": 59,
        "Nd": 60,
        "Pm": 61,
        "Sm": 62,
        "Eu": 63,
        "Gd": 64,
        "Tb": 65,
        "Dy": 66,
        "Ho": 67,
        "Er": 68,
        "Tm": 69,
        "Yb": 70,
        "Lu": 71,
        "Hf": 72,
        "Ta": 73,
        "W": 74,
        "Re": 75,
        "Os": 76,
        "Ir": 77,
        "Pt": 78,
        "Au": 79,
        "Hg": 80,
        "Tl": 81,
        "Pb": 82,
        "Bi": 83,
        "Po": 84,
        "At": 85,
        "Rn": 86,
        "Fr": 87,
        "Ra": 88,
        "Ac": 89,
        "Th": 90,
        "Pa": 91,
        "U": 92,
        "Np": 93,
        "Pu": 94,
        "Am": 95,
        "Cm": 96,
        "Bk": 97,
        "Cf": 98,
        "Es": 99,
        "Fm": 100,
        "Md": 101,
        "No": 102,
        "Lr": 103,
        "Rf": 104,
        "Db": 105,
        "Sg": 106,
        "Bh": 107,
        "Hs": 108,
        "Mt": 109,
        "Ds": 110,
        "Rg": 111,
        "Cn": 112,
        "Nh": 113,
        "Fl": 114,
        "Mc": 115,
        "Lv": 116,
        "Ts": 117,
        "Og": 118,
    }[symbol]

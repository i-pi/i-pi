"""Contains atomic masses, fundamental constants, and unit conversions
to/from atomic units.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import re

from ipi.utils.messages import verbosity, info

try:
    import ase
    from ase.units import Bohr, Hartree, eV, Angstrom, Pascal, kB, Debye, _amu
except ImportError:
    ase = None

__all__ = ["Constants", "Elements", "unit_to_internal", "unit_to_user"]


class Constants(object):
    """Class whose members are fundamental constants.

    Attributes:
        kb: Boltzmann constant.
        hbar: Reduced Planck's constant.
        amu: Atomic mass unit.
    """

    kb = 1.0
    hbar = 1.0
    amu = 1822.8885
    e = 1.0  # elementary charge, electron charge = - |e|


class Elements(dict):
    """Class which contains the mass of different elements.

    Attributes:
        mass_list: A dictionary containing the masses of different elements.
            Has the form {"label": Mass in a.m.u.}. Note that the generic "X"
            label is assumed to be an electron.
    """

    mass_list = {
        "L": 1.0000 / Constants.amu,  # for cavity photon
        "X": 1.0000 / Constants.amu,
        "H": 1.00794,
        "Z": 1.382943,  # an interpolated H-D atom, based on y=1/sqrt(m) scaling
        "D": 2.0141,
        "H2": 2.0160,
        "T": 3.016049,
        "He": 4.002602,
        "Li": 6.9410,
        "Be": 9.012182,
        "B": 10.811,
        "C": 12.0107,
        "C13": 13.00335,
        "C14": 14.003241,
        "N": 14.00674,
        "O": 15.9994,
        "O17": 16.9991315,
        "O18": 17.9991610,
        "F": 18.998403,
        "Ne": 20.1797,
        "Na": 22.989770,
        "Mg": 24.3050,
        "Al": 26.981538,
        "Si": 28.0855,
        "P": 30.973761,
        "S": 32.066,
        "Cl": 35.4527,
        "Cl37": 36.96590259,
        "Ar": 39.9480,
        "K": 39.0983,
        "Ca": 40.078,
        "Sc": 44.955910,
        "Ti": 47.867,
        "V": 50.9415,
        "Cr": 51.9961,
        "Mn": 54.938049,
        "Fe": 55.845,
        "Co": 58.933200,
        "Ni": 58.6934,
        "Cu": 63.546,
        "Zn": 65.39,
        "Ga": 69.723,
        "Ge": 72.61,
        "As": 74.92160,
        "Se": 78.96,
        "Br": 79.904,
        "Br81": 80.9163,
        "Kr": 83.80,
        "Rb": 85.4678,
        "Sr": 87.62,
        "Y": 88.90585,
        "Zr": 91.224,
        "Nb": 92.90638,
        "Mo": 95.94,
        "Tc": 98,
        "Ru": 101.07,
        "Rh": 102.90550,
        "Pd": 106.42,
        "Ag": 107.8682,
        "Cd": 112.411,
        "In": 114.818,
        "Sn": 118.710,
        "Sb": 121.760,
        "Te": 127.60,
        "I": 126.90447,
        "Xe": 131.29,
        "Cs": 132.90545,
        "Ba": 137.327,
        "La": 138.9055,
        "Ce": 140.166,
        "Pr": 140.90765,
        "Nd": 144.24,
        "Pm": 145,
        "Sm": 150.36,
        "Eu": 151.964,
        "Gd": 157.25,
        "Tb": 158.92534,
        "Dy": 162.50,
        "Ho": 164.93032,
        "Er": 167.26,
        "Tm": 168.93241,
        "Yb": 173.04,
        "Lu": 174.967,
        "Hf": 178.49,
        "Ta": 180.9479,
        "W": 183.84,
        "Re": 186.207,
        "Os": 190.23,
        "Ir": 192.217,
        "Pt": 195.078,
        "Au": 196.96655,
        "Hg": 200.59,
        "Tl": 204.3833,
        "Pb": 207.2,
        "Bi": 208.98038,
        "Po": 209,
        "At": 210,
        "Rn": 222,
        "Fr": 223,
        "Ra": 226,
        "Ac": 227,
        "Th": 232.0381,
        "Pa": 231.03588,
        "U": 238.0289,
        "Np": 237,
        "Pu": 244,
        "Am": 243,
        "Cm": 247,
        "Bk": 247,
        "Cf": 251,
        "Es": 252,
        "Fm": 257,
        "Md": 258,
        "No": 259,
        "Lr": 262,
        "Rf": 267,
        "Db": 268,
        "Sg": 269,
        "Bh": 270,
        "Hs": 269,
        "Mt": 278,
        "Ds": 281,
        "Rg": 280,
        "Cn": 285,
        "Uut": 286,
        "Fl": 289,
        "Uup": 288,
        "Lv": 293,
        "Uus": 294,
        "Uuo": 294,
    }

    @classmethod
    def mass(cls, label):
        """Function to access the mass_list attribute.

        Note that this does not require an instance of the Elements class to be
        created, as this is a class method. Therefore using Elements.mass(label)
        will give the mass of the element with the atomic symbol given by label.

        Args:
            label: The atomic symbol of the atom whose mass is required.

        Returns:
            A float giving the mass of the atom with atomic symbol label.
        """

        try:
            return cls.mass_list[label] * Constants.amu
        except KeyError:
            info("Unknown element given, you must specify the mass", verbosity.low)
            return -1.0


# these are the conversion FROM the unit stated to internal (atomic) units
UnitMap = {
    "undefined": {"": 1.00, "automatic": 1.00, "atomic_unit": 1.00},
    "energy": {
        "": 1.00,
        "automatic": 1.00,
        "atomic_unit": 1.00,
        "electronvolt": 0.036749326,
        "j/mol": 0.00000038087989,
        "cal/mol": 0.0000015946679,
        "kelvin": 3.1668152e-06,
    },
    "temperature": {
        "": 1.00,
        "automatic": 1.00,
        "atomic_unit": 1.00,
        "kelvin": 3.1668152e-06,
    },
    "time": {"": 1.00, "automatic": 1.00, "atomic_unit": 1.00, "second": 4.1341373e16},
    "frequency": {  # NB Internally, ANGULAR frequencies are used.
        "": 1.00,
        "automatic": 1.00,
        "atomic_unit": 1.00,
        "inversecm": 4.5563353e-06,
        "hertz*rad": 2.4188843e-17,
        "hertz": 1.5198298e-16,
        "hz": 1.5198298e-16,
        "thz": 1.5198298e-4,
    },
    "electric-field": {  # Hartree/Bohr radius
        # 1Hartree = 27.2113862459 eV
        #    1Bohr = 0.5291772109  A
        "": 1.00,
        "automatic": 1.00,
        "atomic_unit": 1.00,
        "ev/a": 0.019446903811441516,  # 0.5291772109/27.2113862459
    },
    "electric-dipole": {  # electron charge * Bohr
        "": 1.00,
        "automatic": 1.00,
        "atomic_unit": 1.00,
        "eang": 1.8897261,  # electron charge * angstrom
    },
    "ms-momentum": {  # TODO fill up units here (mass-scaled momentum)
        "": 1.00,
        "automatic": 1.00,
        "atomic_unit": 1.00,
    },
    "length": {  # TODO move angles to separate entry;
        "": 1.00,
        "automatic": 1.00,
        "atomic_unit": 1.00,
        "angstrom": 1.8897261,
        "meter": 1.8897261e10,
        "radian": 1.00,
        "degree": 0.017453292519943295,
    },
    "volume": {
        "": 1.00,
        "automatic": 1.00,
        "atomic_unit": 1.00,
        "angstrom3": 6.748334231,
    },
    "velocity": {"": 1.00, "automatic": 1.00, "atomic_unit": 1.00, "m/s": 4.5710289e-7},
    "momentum": {"": 1.00, "automatic": 1.00, "atomic_unit": 1.00},
    "mass": {
        "": 1.00,
        "automatic": 1.00,
        "atomic_unit": 1.00,
        "dalton": 1.00 * Constants.amu,
        "electronmass": 1.00,
    },
    "pressure": {
        "": 1.00,
        "automatic": 1.00,
        "atomic_unit": 1.00,
        "bar": 3.398827377e-9,
        "atmosphere": 3.44386184e-9,
        "pascal": 3.398827377e-14,
        "ev/ang3": 0.0054456877,
    },
    "density": {"": 1.00, "automatic": 1.00, "atomic_unit": 1.00, "g/cm3": 162.67263},
    "force": {
        "": 1.00,
        "automatic": 1.00,
        "atomic_unit": 1.00,
        "newton": 12137805,
        "ev/ang": 0.019446904,
    },
    "hessian": {
        "": 1.00,
        "automatic": 1.00,
        "atomic_unit": 1.00,
        "ev/ang^2": 0.010290858,
    },
}

# Conditionally includes "ase" units for each quantity
if ase is not None:

    # Conversion factors from atomic units to ASE units
    conversion_factors = {
        "undefined": {"ase": 1},
        "energy": {"ase": 1 / (Hartree / eV)},
        "temperature": {"ase": kB},
        "time": {
            "ase": 1
            / ((Bohr / Angstrom) * (1 / (Hartree / eV)) ** 0.5 * (_amu / 1) ** 0.5)
        },
        "frequency": {
            "ase": (Bohr / Angstrom) * (1 / (Hartree / eV)) ** 0.5 * (_amu / 1) ** 0.5
        },
        "electric-field": {"ase": 1 / (Hartree / (Bohr * eV / Angstrom))},
        "electric-dipole": {"ase": 1 / (Bohr / Debye)},
        "length": {"ase": 1 / (Bohr / Angstrom)},
        "volume": {"ase": 1 / ((Bohr**3) / (Angstrom**3))},
        "velocity": {
            "ase": 1 / ((Hartree / eV) ** 0.5 * (Bohr / Angstrom) / (_amu / 1) ** 0.5)
        },
        "momentum": {"ase": 1 / ((Bohr / Angstrom) * (Hartree / eV) ** 0.5)},
        "mass": {"ase": 1},  # Atomic unit of mass to amu
        "pressure": {"ase": 1 / (Hartree / (Bohr**3 * Pascal))},
        "density": {
            "ase": 1 / (_amu / (Angstrom**3))
        },  # Atomic unit of density to g/cm^3
        "force": {"ase": 1 / (Hartree / Bohr)},  # Atomic unit of force to N
        "hessian": {
            "ase": 1 / (Hartree / (Bohr**2))
        },  # Atomic unit of Hessian to eV/Å^2
    }

    # Update UnitMap with ASE conversion factors
    for quantity, factors in conversion_factors.items():
        if quantity in UnitMap:
            UnitMap[quantity].update(factors)


# a list of magnitude prefixes
UnitPrefix = {
    "": 1.0,
    "yotta": 1e24,
    "zetta": 1e21,
    "exa": 1e18,
    "peta": 1e15,
    "tera": 1e12,
    "giga": 1e9,
    "mega": 1e6,
    "kilo": 1e3,
    "milli": 1e-3,
    "micro": 1e-6,
    "nano": 1e-9,
    "pico": 1e-12,
    "femto": 1e-15,
    "atto": 1e-18,
    "zepto": 1e-21,
    "yocto": 1e-24,
}


# builds a RE to match prefix and split out the base unit
UnitPrefixRE = ""
for key in UnitPrefix:
    UnitPrefixRE = UnitPrefixRE + key + "|"
UnitPrefixRE = " *(" + UnitPrefixRE[1:] + ")(.*) *"
UnitPrefixRE = re.compile(UnitPrefixRE)


#
# Atomic units are used EVERYWHERE internally. In order to quickly    #
# interface with any "outside" unit, we set up a simple conversion    #
# library.                                                            #
#


def unit_to_internal(family, unit, number):
    """Converts a number of given dimensions and units into internal units.

    Args:
        family: The dimensionality of the number.
        unit: The units 'number' is originally in.
        number: The value of the parameter in the units 'unit'.

    Returns:
        The number in internal units.

    Raises:
        ValueError: Raised if the user specified units aren't given in the
            UnitMap dictionary.
        IndexError: Raised if the programmer specified dimensionality for the
            parameter isn't in UnitMap. Shouldn't happen, for obvious reasons.
        TypeError: Raised if the prefix is correct, but the base unit is not, in
            the user specified unit string.
    """

    if not (family == "number" or family in UnitMap):
        raise IndexError(family + " is an undefined units kind.")
    if family == "number":
        return number

    if unit == "":
        prefix = ""
        base = ""
    else:
        m = UnitPrefixRE.match(unit)
        if m is None:
            raise ValueError(
                "Unit " + unit + " is not structured with a prefix+base syntax."
            )
        prefix = m.group(1)
        base = m.group(2)

    if prefix not in UnitPrefix:
        raise TypeError(prefix + " is not a valid unit prefix.")
    if not base.lower() in UnitMap[family]:
        print(base.lower(), UnitMap[family])
        raise TypeError(base + " is an undefined unit for kind " + family + ".")

    return number * UnitMap[family][base.lower()] * UnitPrefix[prefix]


def unit_to_user(family, unit, number):
    """Converts a number of given dimensions from internal to user units.

    Args:
        family: The dimensionality of the number.
        unit: The units 'number' should be changed to.
        number: The value of the parameter in internal units.

    Returns:
        The number in the user specified units
    """

    return number / unit_to_internal(family, unit, 1.0)

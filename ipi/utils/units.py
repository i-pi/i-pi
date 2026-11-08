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
    from ase.units import Bohr, Hartree, eV, Angstrom, kB
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

    Element masses are the IUPAC standard atomic weights 2013 (Meija et al.,
    Pure Appl. Chem. 88, 265 (2016), Table 1, uncertainties dropped), the same
    values as ASE's ase.data.atomic_masses. Where the standard weight is an
    interval, the conventional weight of Table 3 is used; elements without
    stable isotopes use the mass of their most stable isotope (Table 4).
    "D", "T" and "Z" are isotope masses, not standard weights.
    """

    mass_list = {
        "L": 1.0000 / Constants.amu,  # for cavity photon
        "X": 1.0000 / Constants.amu,
        "H": 1.008,
        "Z": 1.382943,  # an interpolated H-D atom, based on y=1/sqrt(m) scaling
        "D": 2.0141,
        "H2": 2.0160,
        "T": 3.016049,
        "He": 4.002602,
        "Li": 6.94,
        "Be": 9.0121831,
        "B": 10.81,
        "C": 12.011,
        "C13": 13.00335,
        "C14": 14.003241,
        "N": 14.007,
        "O": 15.999,
        "O17": 16.9991315,
        "O18": 17.9991610,
        "F": 18.998403163,
        "Ne": 20.1797,
        "Na": 22.98976928,
        "Mg": 24.305,
        "Al": 26.9815385,
        "Si": 28.085,
        "P": 30.973761998,
        "S": 32.06,
        "Cl": 35.45,
        "Cl37": 36.96590259,
        "Ar": 39.948,
        "K": 39.0983,
        "Ca": 40.078,
        "Sc": 44.955908,
        "Ti": 47.867,
        "V": 50.9415,
        "Cr": 51.9961,
        "Mn": 54.938044,
        "Fe": 55.845,
        "Co": 58.933194,
        "Ni": 58.6934,
        "Cu": 63.546,
        "Zn": 65.38,
        "Ga": 69.723,
        "Ge": 72.63,
        "As": 74.921595,
        "Se": 78.971,
        "Br": 79.904,
        "Br81": 80.9163,
        "Kr": 83.798,
        "Rb": 85.4678,
        "Sr": 87.62,
        "Y": 88.90584,
        "Zr": 91.224,
        "Nb": 92.90637,
        "Mo": 95.95,
        "Tc": 97.90721,
        "Ru": 101.07,
        "Rh": 102.9055,
        "Pd": 106.42,
        "Ag": 107.8682,
        "Cd": 112.414,
        "In": 114.818,
        "Sn": 118.71,
        "Sb": 121.76,
        "Te": 127.6,
        "I": 126.90447,
        "Xe": 131.293,
        "Cs": 132.90545196,
        "Ba": 137.327,
        "La": 138.90547,
        "Ce": 140.116,
        "Pr": 140.90766,
        "Nd": 144.242,
        "Pm": 144.91276,
        "Sm": 150.36,
        "Eu": 151.964,
        "Gd": 157.25,
        "Tb": 158.92535,
        "Dy": 162.5,
        "Ho": 164.93033,
        "Er": 167.259,
        "Tm": 168.93422,
        "Yb": 173.054,
        "Lu": 174.9668,
        "Hf": 178.49,
        "Ta": 180.94788,
        "W": 183.84,
        "Re": 186.207,
        "Os": 190.23,
        "Ir": 192.217,
        "Pt": 195.084,
        "Au": 196.966569,
        "Hg": 200.592,
        "Tl": 204.38,
        "Pb": 207.2,
        "Bi": 208.9804,
        "Po": 208.98243,
        "At": 209.98715,
        "Rn": 222.01758,
        "Fr": 223.01974,
        "Ra": 226.02541,
        "Ac": 227.02775,
        "Th": 232.0377,
        "Pa": 231.03588,
        "U": 238.02891,
        "Np": 237.04817,
        "Pu": 244.06421,
        "Am": 243.06138,
        "Cm": 247.07035,
        "Bk": 247.07031,
        "Cf": 251.07959,
        "Es": 252.083,
        "Fm": 257.09511,
        "Md": 258.09843,
        "No": 259.101,
        "Lr": 262.11,
        "Rf": 267.122,
        "Db": 268.126,
        "Sg": 271.134,
        "Bh": 270.133,
        "Hs": 269.1338,
        "Mt": 278.156,
        "Ds": 281.165,
        "Rg": 281.166,
        "Cn": 285.177,
        "Uut": 286,
        "Fl": 289.19,
        "Uup": 288,
        "Lv": 293.204,
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
        "v/ang": 0.019446903811441516,  # 0.5291772109/27.2113862459
    },
    "electric-dipole": {  # electron charge * Bohr
        "": 1.00,
        "automatic": 1.00,
        "atomic_unit": 1.00,
        "eang": 1.8897261,  # electron charge * angstrom
        "debye": 0.3934303,  # https://en.wikipedia.org/wiki/Debye#cite_note-10
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
    "action": {
        "": 1.00,
        "automatic": 1.00,
        "atomic_unit": 1.00,
        "joule*second": 9.4825216e33,
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
    "inverse-pressure": {
        "": 1.00,
        "automatic": 1.00,
        "atomic_unit": 1.00,
        "bar^-1": 1.0 / 3.398827377e-9,
        "1/bar": 1.0 / 3.398827377e-9,
        "atmosphere^-1": 1.0 / 3.44386184e-9,
        "1/atmosphere": 1.0 / 3.44386184e-9,
        "pascal^-1": 1.0 / 3.398827377e-14,
        "1/pascal": 1.0 / 3.398827377e-14,
        "ang3/ev": 1.0 / 0.0054456877,
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

    UnitMap["undefined"]["ase"] = 1.0
    UnitMap["energy"]["ase"] = eV / Hartree
    UnitMap["length"]["ase"] = Angstrom / Bohr
    UnitMap["temperature"]["ase"] = kB * eV / Hartree
    UnitMap["time"]["ase"] = (Angstrom / Bohr) * (Constants.amu / (eV / Hartree)) ** 0.5
    UnitMap["velocity"]["ase"] = UnitMap["length"]["ase"] / UnitMap["time"]["ase"]
    UnitMap["force"]["ase"] = UnitMap["energy"]["ase"] / UnitMap["length"]["ase"]
    UnitMap["volume"]["ase"] = UnitMap["length"]["ase"] ** 3
    UnitMap["frequency"]["ase"] = 1 / UnitMap["time"]["ase"]
    UnitMap["action"]["ase"] = UnitMap["energy"]["ase"] * UnitMap["time"]["ase"]
    UnitMap["mass"]["ase"] = (
        UnitMap["energy"]["ase"] / UnitMap["velocity"]["ase"] ** 2
    )  # to check
    UnitMap["pressure"]["ase"] = (
        UnitMap["energy"]["ase"] / UnitMap["volume"]["ase"]
    )  # to check
    UnitMap["inverse-pressure"]["ase"] = 1.0 / UnitMap["pressure"]["ase"]
    UnitMap["momentum"]["ase"] = (
        UnitMap["mass"]["ase"] * UnitMap["velocity"]["ase"]
    )  # to check
    UnitMap["density"]["ase"] = (
        UnitMap["mass"]["ase"] / UnitMap["volume"]["ase"]
    )  # to check
    UnitMap["hessian"]["ase"] = (
        UnitMap["energy"]["ase"] / UnitMap["length"]["ase"] ** 2
    )  # to check
    # since ASe handles charges in atomic_units, dipole and electric field should have the same units as lenth and force, respectively
    UnitMap["electric-dipole"]["ase"] = UnitMap["length"]["ase"]
    UnitMap["electric-field"]["ase"] = UnitMap["force"]["ase"]

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


def unit_to_internal(family, unit, number=1.0):
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


def unit_to_user(family, unit, number=1.0):
    """Converts a number of given dimensions from internal to user units.

    Args:
        family: The dimensionality of the number.
        unit: The units 'number' should be changed to.
        number: The value of the parameter in internal units.

    Returns:
        The number in the user specified units
    """

    return number / unit_to_internal(family, unit)

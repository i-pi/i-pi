"""Small functions/classes providing access to driver PES to be called from driver.py"""

from .dummy import Dummy_driver

__drivers__ = [
    "ase",
    "bath",
    "double_double_well",
    "DW",
    "DW_bath",
    "DW_friction",
    "driverdipole",
    "dummy",
    "extmace",
    "elphmod",
    "harmonic",
    "mace",
    "metatensor",
    "metatomic",
    "MorseHarmonic",
    "pet",
    "psiflow",
    "rascal",
    "so3lr",
    "Spherical_LJ",
    "spline",
    "xtb",
]

# sort alphabetically
__drivers__ = sorted(__drivers__)

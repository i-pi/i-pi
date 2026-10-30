"""Tests unit conversion."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


from ipi.utils import units


def test_case_insensitive():
    angstrom = units.unit_to_internal("length", "angstrom", 1.0)
    Angstrom = units.unit_to_internal("length", "Angstrom", 1.0)
    if angstrom != Angstrom:
        raise ValueError("angstrom != Angstrom")


def test_inverse_pressure_roundtrip():
    """Common compressibility units convert to and from atomic units."""

    value = 4.5e-5
    internal = units.unit_to_internal("inverse-pressure", "bar^-1", value)
    converted = units.unit_to_user("inverse-pressure", "bar^-1", internal)
    assert abs(converted - value) < 1e-15

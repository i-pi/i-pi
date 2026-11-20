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

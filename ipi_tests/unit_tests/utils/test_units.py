"""Tests unit conversion."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


from ipi.utils import units
import numpy as np


def test_case_insensitive():
    angstrom = units.unit_to_internal("length", "angstrom", 1.0)
    Angstrom = units.unit_to_internal("length", "Angstrom", 1.0)
    if angstrom != Angstrom:
        raise ValueError("angstrom != Angstrom")


def test_cyclic_frequency_units_are_not_angular():
    cycles_per_au = units.unit_to_internal("frequency-cyclic", "GHz", 1.0)
    angular_per_au = units.unit_to_internal("frequency", "GHz", 1.0)
    np.testing.assert_allclose(angular_per_au, 2.0 * np.pi * cycles_per_au)

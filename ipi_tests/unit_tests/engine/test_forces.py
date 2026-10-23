"""Tests for the processing of forcefield results."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

import numpy as np

from ipi.engine import forces
from ipi.engine.forces import ForceBead
from ipi.utils.depend import depend_value


class MockForceField:
    """Forcefield with just a name."""

    def __init__(self, name):
        self.name = name


def bead_with_virial(vir, ffname):
    """Returns a ForceBead holding a given forcefield result."""

    bead = ForceBead()
    bead.ff = MockForceField(ffname)
    bead._ufvx = depend_value(name="ufvx", value=[0.0, np.zeros(3), vir, {}])
    return bead


def test_virial_upper_triangular_copy():
    """ForceBead returns the upper triangle of the virial, leaving the result unchanged"""

    vir = np.arange(9.0).reshape(3, 3)
    bead = bead_with_virial(vir.copy(), "copy")
    assert np.array_equal(bead.get_vir(), np.triu(vir))
    assert np.array_equal(bead.ufvx[2], vir)


def test_virial_symmetry_warning(mocker):
    """A non-symmetric full virial gives one warning per forcefield"""

    warn = mocker.patch.object(forces, "warning")
    bead_with_virial(np.arange(9.0).reshape(3, 3), "asymmetric").get_vir()
    bead_with_virial(np.arange(9.0).reshape(3, 3), "asymmetric").get_vir()
    bead_with_virial(np.eye(3) + 1.0, "symmetric").get_vir()
    bead_with_virial(np.triu(np.ones((3, 3))), "upper").get_vir()
    assert warn.call_count == 1

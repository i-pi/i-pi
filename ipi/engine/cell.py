"""Classes which deal with the system box.

Used for implementing the minimum image convention.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np

from ipi.utils.array_backend import xp
from ipi.utils.depend import *
from ipi.utils.mathtools import *

__all__ = ["Cell", "GenericCell"]


class Cell:
    """Base class to represent the simulation cell in a periodic system.

    This class has the base attributes required for either flexible or
    isotropic cell dynamics. Uses an upper triangular lattice vector matrix to
    represent the cell.

    Depend objects:
       h: An array giving the lattice vector matrix.
       ih: An array giving the inverse of the lattice vector matrix.
       V: The volume of the cell.
    """

    def __init__(self, h=None):
        """Initialises base cell class.

        Args:
           h: Optional array giving the initial lattice vector matrix. The
              reference cell matrix is set equal to this. Must be an upper
              triangular 3*3 matrix. Defaults to a 3*3 zeroes matrix.
        """

        if h is None:
            h = np.zeros((3, 3), float)

        self._h = depend_array(name="h", value=h)
        self._ih = depend_array(
            name="ih",
            value=np.zeros((3, 3), float),
            func=self.get_ih,
            dependencies=[self._h],
        )
        self._V = depend_value(name="V", func=self.get_volume, dependencies=[self._h])

    def clone(self):
        return Cell(xp.asarray(dstrip(self.h), copy=True))

    def get_ih(self):
        """Inverts the lattice vector matrix."""

        return invert_ut3x3(self.h)

    def get_volume(self):
        """Calculates the volume of the system box."""

        return det_ut3x3(self.h)

    def apply_pbc(self, atom):
        """Minimum-image wrap of a single atom's position into the cell."""

        s = dstrip(self.ih) @ dstrip(atom.q)
        s = s - xp.round(s)
        return dstrip(self.h) @ s

    def array_pbc(self, pos):
        """Minimum-image wrap of a position array (n_atoms*3,) in place."""

        s = dstrip(pos).reshape(-1, 3)
        s = dstrip(self.ih) @ s.T
        s = s - xp.round(s)
        s = (dstrip(self.h) @ s).T
        pos[:] = s.reshape(-1)

    def minimum_distance(self, atom1, atom2):
        """Minimum-image vector between two atoms.

        Only rigorously accurate for a cubic cell but correct up to a cutoff
        shorter than the smallest width between parallel cell faces.
        """

        s = dstrip(self.ih) @ (dstrip(atom1.q) - dstrip(atom2.q))
        s = s - xp.round(s)
        return dstrip(self.h) @ s


dproperties(Cell, ["h", "ih", "V"])


class GenericCell(Cell):
    """A cell class that does not assume upper-triangular values.

    Depend objects:
       h: An array giving the lattice vector matrix.
       ih: An array giving the inverse of the lattice vector matrix.
       V: The volume of the cell.
    """

    def __init__(self, h=None):
        """Initialises base cell class.

        Args:
           h: Optional array giving the initial lattice vector matrix. The
              reference cell matrix is set equal to this.
        """

        if h is None:
            h = np.zeros((3, 3), float)

        self._h = depend_array(name="h", value=h)
        self._ih = depend_array(
            name="ih",
            value=np.zeros((3, 3), float),
            func=self.get_ih,
            dependencies=[self._h],
        )
        self._V = depend_value(name="V", func=self.get_volume, dependencies=[self._h])

    def clone(self):
        return GenericCell(xp.asarray(dstrip(self.h), copy=True))

    def get_ih(self):
        """Inverts the lattice vector matrix."""

        return xp.linalg.inv(dstrip(self.h))

    def get_volume(self):
        """Calculates the volume of the system box."""

        return xp.linalg.det(dstrip(self.h))


dproperties(GenericCell, ["h", "ih", "V"])

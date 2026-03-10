"""Creates objects that deal with the simulation box."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


from copy import copy

from ipi.engine.cell import *
from ipi.utils.inputvalue import *


__all__ = ["InputCell"]


class InputCell(InputArray):
    """Cell input class.

    Handles generating the appropriate cell class from the xml input file,
    and generating the xml checkpoint tags and data from an instance of the
    object.
    """

    attribs = copy(InputArray.attribs)

    default_help = """
Describes with the cell parameters. Takes as array which can be used to initialize the cell vector matrix.
N.B.: the cell parameters are stored with the lattice vectors in the columns, and the cell must be oriented
in such a way that the array is upper-triangular (i.e. with the first vector aligned along x and the second
vector in the xy plane).     
"""
    default_label = "CELL"

    def __init__(self, help=None, dimension=None, units=None, default=None, dtype=None):
        """Initializes InputCell.

        Just calls the parent initialization function with appropriate arguments.
        """

        super(InputCell, self).__init__(
            dtype=float, dimension="length", default=default, help=help
        )

    def store(self, cell):
        """Takes a Cell instance and stores of minimal representation of it.

        Args:
           cell: A cell object.
        """

        super(InputCell, self).store(cell.h)
        self.shape.store((3, 3))

    def fetch(self):
        """Creates a cell object.

        Returns:
           A cell object of the appropriate type and with the appropriate
           properties given the attributes of the InputCell object.
        """

        h = super(InputCell, self).fetch()
        h.shape = (3, 3)

        return Cell(h=h)

    def check(self):
        """Checks for optional parameters."""

        super(InputCell, self).check()

        h = self.value.reshape(9)

        if (h[[3, 6, 7]] ** 2).sum() > 1e-20:
            raise ValueError("The cell matrix must be upper triangular.")

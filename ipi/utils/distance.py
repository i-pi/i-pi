"""Python version of vector_separation() in drivers/distance.f90"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

import numpy as np

__all__ = ["vector_separation"]


def vector_separation(cell_h, cell_ih, qi, qj):
    """Calculates the vector separating two atoms.

    Note that minimum image convention is used, so only the image of
    atom j that is the shortest distance from atom i is considered.

    Also note that while this may not work if the simulation
    box is highly skewed from orthorhombic, as
    in this case it is possible to return a distance less than the
    nearest neighbour distance. However, this will not be of
    importance unless the cut-off radius is more than half the
    width of the shortest face-face distance of the simulation box,
    which should never be the case.

    Args:
       cell_h: The simulation box cell vector matrix.
       cell_ih: The inverse of the simulation box cell vector matrix.
       qi: The position vector of atom i.
       qj: The position vectors of one or many atoms j shaped as (N, 3).
    Returns:
       dij: The vectors separating atoms i and {j}.
       rij: The distances between atoms i and {j}.
    """

    sij = np.dot(cell_ih, (qi - qj).T)  # column vectors needed
    sij -= np.rint(sij)

    dij = np.dot(cell_h, sij).T  # back to i-pi shape
    rij = np.linalg.norm(dij, axis=1)

    return dij, rij

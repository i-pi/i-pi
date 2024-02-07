"""Classes which deal with classical atoms.

Used for holding information about the atoms, including their positions, masses
momenta and kinetic energy. Has separate classes for accessing the global
arrays of atoms and for individual atoms.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np

from ipi.utils.depend import *


__all__ = ["Atoms", "Atom"]


class Atom:
    """Represent an atom, with position, velocity, mass and related properties.

    This is actually only an interface to the Atoms class, i.e. only stores
    views of the large arrays which contain all the coordinates.

    Attributes:
       kin: The kinetic energy of the atom.
       kstress: The contribution of the atom to the kinetic stress tensor.

    Depend objects:
       p: The three components of the momentum of the atom.
       q: The three components of the position of the atom.
       m: The mass of the atom.
       name: The name of the atom.
       m3: An array of 3 elements with each element being the mass of the atom.
          Used when each degree of freedom needs to be divided by the mass.
    """

    def __init__(self, system, index):
        """Initialises Atom.

        Args:
           system: An Atoms object containing the required atom.
           index: An integer giving the index of the required atom in the atoms
              list. Note that indices start from 0.
        """

        self.p = system.p[3 * index : 3 * index + 3]
        self.q = system.q[3 * index : 3 * index + 3]
        self.m = system.m[index : index + 1]
        self.name = system.names[index : index + 1]
        self.m3 = system.m3[3 * index : 3 * index + 3]

    @property
    def kin(self):
        """Calculates the contribution of the atom to the kinetic energy."""

        return np.dot(self.p, self.p) / (2.0 * self.m)

    @property
    def kstress(self):
        """Calculates the contribution of the atom to the kinetic stress
        tensor.
        """

        p = dstrip(self.p)
        ks = np.zeros((3, 3), float)
        for i in range(3):
            for j in range(i, 3):
                ks[i, j] = p[i] * p[j]
        return ks / self.m


# dproperties(Atom, ["p", "q", "m", "m3", "name"])


class Atoms:
    """Storage for the atoms' positions, masses and velocities.

    Everything is stored as 3*n sized contiguous arrays,
    and a convenience-access is provided through a list of Atom objects.

    Attributes:
       natoms: The number of atoms.

    Depend objects:
       p: An array giving the components of the atom momenta.
       q: An array giving the components of the atom position.
       m: An array giving the atom masses.
       names: An array giving the atom names.
       m3: An array of 3*n elements where each element of m has been copied
          three times. Used when each degree of freedom needs to be divided
          by the mass.
       M: The total mass of all the atoms.
       kin: The total kinetic energy of the atoms. Depends on p and m3.
       kstress: The contribution of the atoms to the kinetic stress tensor.
          Depends on px, py, pz and m.
       qx: An array giving the x components of the positions.
       qy: An array giving the y components of the positions.
       qz: An array giving the z components of the positions.
       px: An array giving the x components of the momenta.
       py: An array giving the y components of the momenta.
       pz: An array giving the z components of the momenta.
    """

    def __init__(self, natoms, _prebind=None):
        """Initialises Atoms.

        Each replica and the centroid coordinate are all held as Atoms objects,
        and so slices of the global position and momentum arrays must be used in
        the initialisation so that they always agree with each other.

        Args:
           natoms: An integer giving the number of atoms.
           _prebind: An optional tuple of four elements; a depend_array of length
              3*natoms for the positions, another for the momenta, a depend_array
              of length natoms for the masses and another for the names.
        """

        self.natoms = natoms

        if _prebind is None:
            self._q = depend_array(name="q", value=np.zeros(3 * natoms, float))
            self._p = depend_array(name="p", value=np.zeros(3 * natoms, float))
            self._m = depend_array(name="m", value=np.zeros(natoms, float))
            self._names = depend_array(
                name="names", value=np.zeros(natoms, np.dtype("|U6"))
            )
        else:
            self._q = _prebind[0]
            self._p = _prebind[1]
            self._m = _prebind[2]
            self._names = _prebind[3]

        self._m3 = depend_array(
            name="m3",
            value=np.zeros(3 * natoms, float),
            func=self.mtom3,
            dependencies=[self._m],
        )

        self._M = depend_value(name="M", func=self.get_msum, dependencies=[self._m])
        self._kin = depend_value(
            name="kin", func=self.get_kin, dependencies=[self._p, self._m3]
        )
        self._kstress = depend_value(
            name="kstress", func=self.get_kstress, dependencies=[self._p, self._m]
        )

    def copy(self):
        """Creates a new Atoms object.

        Returns:
           An Atoms object with the same q, p, m and names arrays as the original.
        """

        newat = Atoms(self.natoms)
        newat.q[:] = self.q
        newat.p[:] = self.p
        newat.m[:] = self.m
        newat.names[:] = self.names
        return newat

    def __len__(self):
        """Length function.

        This is called whenever the standard function len(atoms) is used.

        Returns:
           The number of atoms.
        """

        return self.natoms

    def __iter__(self):
        """Iterator.

        This is called whenever one iterates over an Atoms object.

        Returns:
           Itertor over all atoms in this Atoms object.
        """

        for index in range(len(self)):
            yield Atom(self, index)

    def __getitem__(self, index):
        """Overwrites standard getting function.

        This is called whenever the standard function atoms[index] is used.
        Returns an Atom object with the appropriate position and momenta arrays.
        Note that they are dynamically generated each time an Atom needs to be
        accessed, as this reduces the number of depend objects that need to be
        held at any one time.

        Args:
           index: The index of the atom to be accessed.

        Returns:
           The atom given by the index.
        """

        return Atom(self, index)

    def __setitem__(self, index, value):
        """Overwrites standard setting function.

        This is called whenever the standard function atoms[index]=value is used.
        Changes the position and momenta of the appropriate slice of the global
        position and momentum arrays to those given by value.
        Note that they are dynamically generated each time an Atom needs to be
        accessed, as this reduces the number of depend objects that need to be
        held at any one time.

        Args:
           index: The atom to be changed.
           value: The Atom object that holds the new values.
        """

        pat = Atom(self, index)
        pat.p = value.p
        pat.q = value.q
        pat.m = value.m
        pat.name = value.name

    def get_msum(self):
        """Calculates the total mass."""

        return self.m.sum()

    def mtom3(self):
        """Returns a 3*n mass array.

        Returns:
           An array of 3*n elements where each element of m has been copied
           three times. Used when each degree of freedom needs to be divided
           by the mass.
        """

        m3 = np.zeros(3 * self.natoms, float)
        m3[0 : 3 * self.natoms : 3] = self.m
        m3[1 : 3 * self.natoms : 3] = m3[0 : 3 * self.natoms : 3]
        m3[2 : 3 * self.natoms : 3] = m3[0 : 3 * self.natoms : 3]
        return m3

    def get_kin(self):
        """Calculates the total kinetic energy of the system."""

        p = dstrip(self.p)
        return 0.5 * np.dot(p, p / dstrip(self.m3))

    def get_kstress(self):
        """Calculates the total contribution of the atoms to the kinetic stress
        tensor -- not volume-scaled
        """

        p = dstrip(self.p)
        m = dstrip(self.m)
        px = p[0::3]
        py = p[1::3]
        pz = p[2::3]

        ks = np.zeros((3, 3), float)
        ks[0, 0] = np.dot(px, px / m)
        ks[1, 1] = np.dot(py, py / m)
        ks[2, 2] = np.dot(pz, pz / m)
        ks[0, 1] = np.dot(px, py / m)
        ks[0, 2] = np.dot(px, pz / m)
        ks[1, 2] = np.dot(py, pz / m)
        return ks


dproperties(Atoms, ["p", "q", "m", "m3", "names", "M", "kin", "kstress"])

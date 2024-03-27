"""Classes which deal with all the PI beads.

Used for holding information about the beads, including their positions, masses
momenta and kinetic energy. Has different objects for the position and normal
mode representations, and has a special centroid atoms object for when the
centroid coordinate is required.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np

from ipi.utils.depend import *
from ipi.engine.atoms import Atoms


__all__ = ["Beads"]


class Beads:
    """Storage for the beads positions and velocities.

    Everything is stored as (nbeads,3*natoms) sized contiguous arrays,
    and a convenience-access to each replica of the system is provided through a
    list of Atoms objects. Contains arrays of both the normal mode representation
    and the position representation, and various sized arrays for the atom
    labels and masses. Also contains the potential and force between
    neighbouring replicas.

    Attributes:
       natoms: The number of atoms.
       nbeads: The number of beads.
       _blist: A list of Atoms objects for each replica of the system. Each
          replica is assumed to have the same mass and atom label.
       centroid: An atoms object giving the centroid coordinate of the beads.

    Depend objects:
       names: An array giving the atom names.
       m: An array giving the atom masses.
       m3: An array giving the mass associated with each degree of freedom.
       sm3: An array giving the square root of m3.
       q: An array giving all the bead positions.
       p: An array giving all the bead momenta.
       qc: An array giving the centroid positions. Depends on qnm.
       pc: An array giving the centroid momenta. Depends on pnm.
       vpath: The spring potential between the beads, divided by omegan**2.
          Depends on q.
       fpath: The spring force between the beads, divided by omegan**2.
          Depends on q.
       kins: A list of the kinetic energy of each replica.
       kin: The total kinetic energy of the system. Note that this is not the
          same as the estimate of the kinetic energy of the system, which is
          contained in the properties module.
       kstress: The total kinetic stress tensor for the system.
       rg: An array giving the radius of gyration of each atom.
    """

    def __init__(self, natoms, nbeads):
        """Initialises Beads.

        Args:
           natoms: Number of atoms.
           nbeads: Number of beads.
        """

        self.resize(natoms, nbeads)

    def resize(self, natoms, nbeads):
        """Creates all the data arrays needed in the simulation.

        Effectively initializes the whole Beads object, according to the
        specified number of atoms and beads. Is also used, as the name suggests,
        to resize the data to a new number of beads when this is necessary, for
        example in initialization from a simulation with a different number of
        beads.

        Also creates, or recreates, the dependency network, as this requires
        the data arrays to be created for it to work.

        Args:
           natoms: The number of atoms.
           nbeads: The number of beads.
        """

        self.natoms = natoms
        self.nbeads = nbeads

        self._names = depend_array(
            name="names", value=np.zeros(natoms, np.dtype("|U6"))
        )

        # atom masses, and mass-related arrays
        self._m = depend_array(
            name="m", value=np.zeros(natoms, float)
        )  # this is the prototype mass array (just one independent of bead n)
        self._m3 = depend_array(
            name="m3",
            value=np.zeros(
                (nbeads, 3 * natoms), float
            ),  # this is m conveniently replicated to be (nb,3*nat)
            func=self.mtom3,
            dependencies=[self._m],
        )
        self._sm3 = depend_array(
            name="sm3",
            value=np.zeros(
                (nbeads, 3 * natoms), float
            ),  # this is just the square root of m3
            func=self.m3tosm3,
            dependencies=[self._m3],
        )

        # positions and momenta. bead representation, base storage used everywhere
        self._q = depend_array(name="q", value=np.zeros((nbeads, 3 * natoms), float))
        self._p = depend_array(name="p", value=np.zeros((nbeads, 3 * natoms), float))

        # position and momentum of the centroid
        self._qc = depend_array(
            name="qc",
            value=np.zeros(3 * natoms, float),
            func=self.get_qc,
            dependencies=[self._q],
        )
        self._pc = depend_array(
            name="pc",
            value=np.zeros(3 * natoms, float),
            func=self.get_pc,
            dependencies=[self._p],
        )

        # path springs potential and force
        self._vpath = depend_value(
            name="vpath", func=self.get_vpath, dependencies=[self._q, self._m3]
        )
        self._fpath = depend_array(
            name="fpath",
            value=np.zeros((nbeads, 3 * natoms), float),
            func=self.get_fpath,
            dependencies=[self._q],
        )

        # create proxies to access the individual beads as Atoms objects
        # TODO: ACTUALLY THIS IS ONLY USED HERE METHINK, SO PERHAPS WE COULD REMOVE IT TO DECLUTTER THE CODE.
        self._blist = [
            Atoms(natoms, _prebind=(self.q[i, :], self.p[i, :], self.m, self.names))
            for i in range(nbeads)
        ]

        # kinetic energies of thhe beads, and total (classical) kinetic stress tensor
        self._kins = depend_array(
            name="kins",
            value=np.zeros(nbeads, float),
            func=self.kin_gather,
            dependencies=[b._kin for b in self._blist],
        )
        self._kin = depend_value(
            name="kin", func=self.get_kin, dependencies=[self._kins]
        )
        self._kstress = depend_array(
            name="kstress",
            value=np.zeros((3, 3), float),
            func=self.get_kstress,
            dependencies=[b._kstress for b in self._blist],
        )

    def copy(self, nbeads=-1):
        """Creates a new beads object with newP <= P beads from the original.

        Returns:
           A Beads object with the first newP q, p, m and names arrays as the original.
        """

        if nbeads > self.nbeads:
            raise ValueError("Cannot copy to an object with larger number of beads")
        elif nbeads == -1:
            nbeads = self.nbeads

        newbd = Beads(self.natoms, nbeads)
        newbd.q[:] = self.q[:nbeads]
        newbd.p[:] = self.p[:nbeads]
        newbd.m[:] = self.m
        newbd.names[:] = self.names
        return newbd

    def m3tosm3(self):
        """Takes the mass array and returns the square rooted mass array."""

        return np.sqrt(dstrip(self.m3))

    def mtom3(self):
        """Takes the mass array for each bead and returns one with an element
        for each degree of freedom.

        Returns:
           An array of size (nbeads,3*natoms), with each element corresponding
           to the mass associated with the appropriate degree of freedom in q.
        """

        m3 = np.zeros((self.nbeads, 3 * self.natoms), float)
        m3[:, 0 : 3 * self.natoms : 3] = self.m
        m3[:, 1 : 3 * self.natoms : 3] = m3[:, 0 : 3 * self.natoms : 3]
        m3[:, 2 : 3 * self.natoms : 3] = m3[:, 0 : 3 * self.natoms : 3]
        return m3

    def get_qc(self):
        """Gets the centroid coordinates."""

        return np.sum(dstrip(self.q), axis=0) / self.nbeads
        # return np.dot(np.ones(self.nbeads, float), dstrip(self.q)) / float(self.nbeads)

    def get_pc(self):
        """Gets the centroid momenta."""

        return np.sum(dstrip(self.p), axis=0) / self.nbeads
        # return np.dot(np.ones(self.nbeads, float), dstrip(self.p)) / float(self.nbeads)

    def kin_gather(self):
        """Gets the kinetic energy for all the replicas.

        Returns:
           A list of the kinetic energy for each system.
        """

        return np.array([b.kin for b in self._blist])

    def get_kin(self):
        """Gets the total kinetic energy of all the replicas.

        Note that this does not correspond to the total kinetic energy estimate
        for the system.

        Returns:
           The sum of the kinetic energy of each replica.
        """

        return self.kins.sum()

    def get_kstress(self):
        """Calculates the total kinetic stress tensor of all the replicas.

        Note that this does not correspond to the quantum kinetic stress tensor
        estimate for the system.

        Returns:
           The sum of the kinetic stress tensor of each replica.
        """

        ks = np.zeros((3, 3), float)
        for b in range(self.nbeads):
            ks += self[b].kstress
        return ks

    def get_vpath(self):
        """Calculates the spring potential between the replicas.

        Note that this is actually the harmonic potential without being
        multiplied by the factor omegan**2, which is only available in the
        ensemble as the temperature is required to calculate it.
        """

        epath = 0.0
        q = dstrip(self.q)
        m = dstrip(self.m3)[0]
        for b in range(self.nbeads):
            if b > 0:
                dq = q[b, :] - q[b - 1, :]
            else:
                dq = q[b, :] - q[self.nbeads - 1, :]
            epath += np.dot(dq, m * dq)
        print(
            "WARNING: RETURNS AN INCORRECT RESULT IF OPEN PATHS ARE BEING USED. CALL NM.VSPRING INSTEAD!!"
        )
        return epath * 0.5

    def get_fpath(self):
        """Calculates the spring force between the replicas.

        Note that this is actually the harmonic force without being
        multiplied by the factor omegan**2, which is only available in the
        ensemble as the temperature is required to calculate it.
        """

        nbeads = self.nbeads
        natoms = self.natoms
        f = np.zeros((nbeads, 3 * natoms), float)

        q = dstrip(self.q)
        m = dstrip(self.m3)[0]
        for b in range(nbeads):
            if b > 0:
                dq = q[b, :] - q[b - 1, :]
            else:
                dq = q[b, :] - q[self.nbeads - 1, :]
            dq *= m
            f[b] -= dq
            if b > 0:
                f[b - 1] += dq
            else:
                f[nbeads - 1] += dq
        return f

    # A set of functions to access individual beads as Atoms objects
    def __len__(self):
        """Length function.

        This is called whenever the standard function len(beads) is used.

        Returns:
           The number of beads.
        """

        return self.nbeads

    def __getitem__(self, index):
        """Overwrites standard getting function.

        This is called whenever the standard function beads[index] is used.
        Returns an Atoms object with the appropriate position and momenta arrays.

        Args:
           index: The index of the replica of the system to be accessed.

        Returns:
           The replica of the system given by the index.
        """

        return self._blist[index]

    def __setitem__(self, index, value):
        """Overwrites standard setting function.

        This is called whenever the standard function beads[index]=value is used.
        Changes the position and momenta of the appropriate slice of the global
        position and momentum arrays to those given by value.

        Args:
           index: The replica of the system to be changed.
           value: The Atoms object that holds the new values.
        """

        self._blist[index].p[:] = value.p
        self._blist[index].q[:] = value.q
        self._blist[index].m[:] = value.m
        self._blist[index].names[:] = value.names


dproperties(
    Beads,
    [
        "p",
        "q",
        "pc",
        "qc",
        "m",
        "m3",
        "names",
        "sm3",
        "kin",
        "kstress",
        "vpath",
        "fpath",
    ],
)

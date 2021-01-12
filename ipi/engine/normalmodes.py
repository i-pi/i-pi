"""Contains the classes that deal with the normal mode representation.

Deals with the normal mode transformation, including the complications
introduced by PA-CMD when the bead masses are rescaled. Also deals with
the change in the dynamics introduced by this mass-scaling, and has its
own functions to calculate the kinetic energy, and the exact propagator
in the normal mode representation under the ring polymer Hamiltonian.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np

from ipi.utils.depend import *
from ipi.utils import units
from ipi.utils import nmtransform
from ipi.utils.messages import verbosity, warning, info
from ipi.utils.exchange import *

__all__ = ["NormalModes"]


class NormalModes(dobject):

    """Handles the path normal modes.

    Normal-modes transformation, determination of path frequencies,
    dynamical mass matrix change, etc.

    Attributes:
       natoms: The number of atoms.
       nbeads: The number of beads.
       beads: The beads object for which the normal mode transformation should
          be done.
       ensemble: The ensemble object, specifying the temperature to hold the
          system to.
       motion: The motion object that will need normal-mode transformation and propagator
       transform: A nm_trans object that contains the functions that are
          required for the normal mode transformation.

    Depend objects:
       mode: A string specifying how the bead masses are chosen.
       transform_method: A string specifying how to do the normal mode
          transformation.
       nm_freqs: An array that specifies how the normal mode frequencies
          of the ring polymers are to be calculated, and thus how the
          bead masses should be chosen.
       qnm: The bead positions in the normal mode representation. Depends on
          beads.q.
       pnm: The bead momenta in the normal mode representation. Depends on
          beads.p.
       omegan: The effective vibrational frequency for the interaction
          between the replicas. Depends on the simulation temperature.
       omegan2: omegan**2.
       omegak: The normal mode frequencies for the free ring polymer.
          Depends on omegan.
       prop_pq: An array holding the exact normal mode propagator for the
          free ring polymer, using mass scaled coordinates.
          See J. Chem. Phys. 133, 124101 (2010). Depends on the bead masses
          and the timestep.
       nm_factor: An array of dynamical mass factors associated with each of
          the normal modes. Depends on nm_freqs and mode.
       dynm3: An array that gives the dynamical masses of individual atoms in the
          normal modes representation. Depends on nm_factor and beads.m3.
       dynomegak: The scaled vibrational frequencies. Depends on nm_factor and
          omegak.
       kins: A list of the kinetic energy for each normal mode, as
          calculated in the normal mode representation, using the
          dynamical mass factors. Depends on beads.sm3, beads.p and nm_factor.
       kin: The total kinetic energy, as calculated in the normal mode
          representation, using the dynamical mass factors.
       kstress: The kinetic stress tensor, as calculated in the normal mode
          representation, using the dynamical mass factors. Depends on
          beads.sm3, beads.p and nm_factor.
    """

    def __init__(
        self,
        mode="rpmd",
        transform_method="fft",
        propagator="exact",
        freqs=None,
        open_paths=None,
        bosons=None,
        dt=1.0,
        nmts=1,
    ):
        """Initializes NormalModes.

        Sets the options for the normal mode transform.

        Args:
           mode: A string specifying how to calculate the bead masses.
           transform_method: A string specifying how to do the normal mode
              transformation.
           freqs: A list of data used to calculate the dynamical mass factors.
        """

        if freqs is None:
            freqs = []
        if open_paths is None:
            open_paths = []
        if bosons is None:
            bosons = []
        self.open_paths = np.asarray(open_paths, int)
        self.bosons = np.asarray(bosons, int)
        dself = dd(self)
        dself.nmts = depend_value(name="nmts", value=nmts)
        dself.dt = depend_value(name="dt", value=dt)
        dself.mode = depend_value(name="mode", value=mode)
        dself.transform_method = depend_value(
            name="transform_method", value=transform_method
        )
        dself.propagator = depend_value(name="propagator", value=propagator)
        dself.nm_freqs = depend_array(name="nm_freqs", value=np.asarray(freqs, float))

    def copy(self, freqs=None):
        """Creates a new beads object from the original.

        Returns:
           A Beads object with the same q, p, m and names arrays as the original.
        """

        if freqs is None:
            freqs = self.nm_freqs.copy()

        newnm = NormalModes(
            self.mode,
            self.transform_method,
            self.propagator,
            freqs,
            self.open_paths,
            self.bosons,
            self.dt,
            self.nmts,
        )
        return newnm

    def bind(self, ensemble, motion, beads=None, forces=None):
        """Initializes the normal modes object and binds to beads and ensemble.

        Do all the work down here as we need a full-formed necklace and ensemble
        to know how this should be done.

        Args:
           beads: A beads object to be bound.
           ensemble: An ensemble object to be bound.
        """

        self.ensemble = ensemble
        self.motion = motion
        if beads is None:
            self.beads = motion.beads
        else:
            self.beads = beads

        self.forces = forces
        self.nbeads = beads.nbeads
        self.natoms = beads.natoms

        # if ( (len(self.bosons) > 0) and (len(self.bosons) < self.natoms) ):
        # raise(IOError("@NormalModes : Currently, only full bosonic/distinguishable simulations are allowed"))
        if len(self.bosons) > self.natoms:
            raise IOError

        dself = dd(self)

        # stores a reference to the bound beads and ensemble objects
        self.ensemble = ensemble
        dpipe(dd(motion).dt, dself.dt)

        # sets up what's necessary to perform nm transformation.
        if self.nbeads == 1:  # classical trajectory! don't waste time doing anything!
            self.transform = nmtransform.nm_noop(nbeads=self.nbeads)
        elif self.transform_method == "fft":
            self.transform = nmtransform.nm_fft(
                nbeads=self.nbeads, natoms=self.natoms, open_paths=self.open_paths
            )
        elif self.transform_method == "matrix":
            self.transform = nmtransform.nm_trans(
                nbeads=self.nbeads, open_paths=self.open_paths
            )

        # creates arrays to store normal modes representation of the path.
        # must do a lot of piping to create "ex post" a synchronization between the beads and the nm
        sync_q = synchronizer()
        sync_p = synchronizer()
        dself.qnm = depend_array(
            name="qnm",
            value=np.zeros((self.nbeads, 3 * self.natoms), float),
            func={"q": (lambda: self.transform.b2nm(dstrip(self.beads.q)))},
            synchro=sync_q,
        )
        dself.pnm = depend_array(
            name="pnm",
            value=np.zeros((self.nbeads, 3 * self.natoms), float),
            func={"p": (lambda: self.transform.b2nm(dstrip(self.beads.p)))},
            synchro=sync_p,
        )

        # must overwrite the functions
        dd(self.beads).q._func = {
            "qnm": (lambda: self.transform.nm2b(dstrip(self.qnm)))
        }
        dd(self.beads).p._func = {
            "pnm": (lambda: self.transform.nm2b(dstrip(self.pnm)))
        }
        dd(self.beads).q.add_synchro(sync_q)
        dd(self.beads).p.add_synchro(sync_p)

        # also within the "atomic" interface to beads
        for b in range(self.nbeads):
            dd(self.beads._blist[b]).q._func = {
                "qnm": (lambda: self.transform.nm2b(dstrip(self.qnm)))
            }
            dd(self.beads._blist[b]).p._func = {
                "pnm": (lambda: self.transform.nm2b(dstrip(self.pnm)))
            }
            dd(self.beads._blist[b]).q.add_synchro(sync_q)
            dd(self.beads._blist[b]).p.add_synchro(sync_p)

        # finally, we mark the beads as those containing the set positions
        dd(self.beads).q.update_man()
        dd(self.beads).p.update_man()

        # forces can be converted in nm representation, but here it makes no sense to set up a sync mechanism,
        # as they always get computed in the bead rep
        if self.forces is not None:
            dself.fnm = depend_array(
                name="fnm",
                value=np.zeros((self.nbeads, 3 * self.natoms), float),
                func=(lambda: self.transform.b2nm(dstrip(self.forces.f))),
                dependencies=[dd(self.forces).f],
            )
        else:  # have a fall-back plan when we don't want to initialize a force mechanism, e.g. for ring-polymer initialization
            dself.fnm = depend_array(
                name="fnm",
                value=np.zeros((self.nbeads, 3 * self.natoms), float),
                func=(
                    lambda: depraise(
                        ValueError(
                            "Cannot access NM forces when initializing the NM object without providing a force reference!"
                        )
                    )
                ),
                dependencies=[],
            )

        # create path-frequencies related properties
        dself.omegan = depend_value(
            name="omegan", func=self.get_omegan, dependencies=[dd(self.ensemble).temp]
        )
        dself.omegan2 = depend_value(
            name="omegan2", func=self.get_omegan2, dependencies=[dself.omegan]
        )
        dself.omegak = depend_array(
            name="omegak",
            value=np.zeros(self.beads.nbeads, float),
            func=self.get_omegak,
            dependencies=[dself.omegan],
        )
        dself.omegak2 = depend_array(
            name="omegak2",
            value=np.zeros(self.beads.nbeads, float),
            func=(lambda: self.omegak ** 2),
            dependencies=[dself.omegak],
        )

        # Add o_omegak to calculate the freq in the case of open path
        dself.o_omegak = depend_array(
            name="o_omegak",
            value=np.zeros(self.beads.nbeads, float),
            func=self.get_o_omegak,
            dependencies=[dself.omegan],
        )

        # sets up "dynamical" masses -- mass-scalings to give the correct RPMD/CMD dynamics
        dself.nm_factor = depend_array(
            name="nm_factor",
            value=np.zeros(self.nbeads, float),
            func=self.get_nmm,
            dependencies=[dself.nm_freqs, dself.mode],
        )
        # add o_nm_factor for the dynamical mass in the case of open paths
        dself.o_nm_factor = depend_array(
            name="nmm",
            value=np.zeros(self.nbeads, float),
            func=self.get_o_nmm,
            dependencies=[dself.nm_freqs, dself.mode],
        )
        dself.dynm3 = depend_array(
            name="dynm3",
            value=np.zeros((self.nbeads, 3 * self.natoms), float),
            func=self.get_dynm3,
            dependencies=[dself.nm_factor, dd(self.beads).m3],
        )
        dself.dynomegak = depend_array(
            name="dynomegak",
            value=np.zeros(self.nbeads, float),
            func=self.get_dynwk,
            dependencies=[dself.nm_factor, dself.omegak],
        )

        dself.dt = depend_value(name="dt", value=1.0)
        dpipe(dd(self.motion).dt, dself.dt)
        dself.prop_pq = depend_array(
            name="prop_pq",
            value=np.zeros((self.beads.nbeads, 2, 2)),
            func=self.get_prop_pq,
            dependencies=[dself.omegak, dself.nm_factor, dself.dt, dself.propagator],
        )
        dself.o_prop_pq = depend_array(
            name="o_prop_pq",
            value=np.zeros((self.beads.nbeads, 2, 2)),
            func=self.get_o_prop_pq,
            dependencies=[
                dself.o_omegak,
                dself.o_nm_factor,
                dself.dt,
                dself.propagator,
            ],
        )

        # if the mass matrix is not the RPMD one, the MD kinetic energy can't be
        # obtained in the bead representation because the masses are all mixed up
        dself.kins = depend_array(
            name="kins",
            value=np.zeros(self.nbeads, float),
            func=self.get_kins,
            dependencies=[dself.pnm, dd(self.beads).sm3, dself.nm_factor],
        )
        dself.kin = depend_value(
            name="kin", func=self.get_kin, dependencies=[dself.kins]
        )
        dself.kstress = depend_array(
            name="kstress",
            value=np.zeros((3, 3), float),
            func=self.get_kstress,
            dependencies=[dself.pnm, dd(self.beads).sm3, dself.nm_factor],
        )

        # Array that holds both vspring and fspring for bosons
        dself.vspring_and_fspring_B = depend_value(
            name="v_and_fs_B",
            value=[None, None],
            func=self.get_vspring_and_fspring_B,
            dependencies=[dself.beads.q, dself.beads.m3, dself.omegan2],
        )

        # spring energy, calculated in normal modes
        dself.vspring = depend_value(
            name="vspring",
            value=0.0,
            func=self.get_vspring,
            dependencies=[
                dself.qnm,
                dself.omegak,
                dself.o_omegak,
                dd(self.beads).m3,
                dself.vspring_and_fspring_B,
            ],
        )

        # spring forces on normal modes
        dself.fspringnm = depend_array(
            name="fspringnm",
            value=np.zeros((self.nbeads, 3 * self.natoms), float),
            func=self.get_fspringnm,
            dependencies=[dself.qnm, dself.omegak, dd(self.beads).m3],
        )

        # spring forces on beads, transformed from normal modes
        dself.fspring = depend_array(
            name="fs",
            value=np.zeros((self.nbeads, 3 * self.natoms), float),
            # func=(lambda: self.transform.nm2b(dstrip(self.fspringnm))),
            func=self.get_fspring,
            dependencies=[dself.fspringnm, dself.vspring_and_fspring_B],
        )

    def get_fspringnm(self):
        """Returns the spring force calculated in NM representation."""

        return -self.beads.m3 * self.omegak[:, np.newaxis] ** 2 * self.qnm

    def get_vspring(self):
        """Returns the spring energy calculated in NM representation for distinguishable particles.
        For bosons, get the first element of vspring_and_fspring_B[0]
        For a mixture of both, calculate separately and combine.
        """

        if self.nbeads == 1:
            return 0.0

        if len(self.bosons) == 0:
            sqnm = dstrip(self.qnm) * dstrip(self.beads.sm3)
            q2 = (sqnm ** 2).sum(axis=1)

            vspring = (self.omegak2 * q2).sum()

            for j in self.open_paths:
                vspring += (
                    self.beads.m[j]
                    * (self.o_omegak ** 2 - self.omegak ** 2)
                    * (
                        self.qnm[:, 3 * j] ** 2
                        + self.qnm[:, 3 * j + 1] ** 2
                        + self.qnm[:, 3 * j + 2] ** 2
                    )
                ).sum()

            return vspring * 0.5

        elif len(self.bosons) is self.natoms:
            return self.vspring_and_fspring_B[0]
        else:
            # Sum over only those particles who are distinguishable.
            vspring = 0.0

            notbosons = list(set(range(self.natoms)) - set(self.bosons))
            for j in notbosons:
                vspring += (
                    self.beads.m[j]
                    * self.omegak ** 2
                    * (
                        self.qnm[:, 3 * j] ** 2
                        + self.qnm[:, 3 * j + 1] ** 2
                        + self.qnm[:, 3 * j + 2] ** 2
                    )
                ).sum()

            return vspring * 0.5 + self.vspring_and_fspring_B[0]

    def get_omegan(self):
        """Returns the effective vibrational frequency for the interaction
        between replicas.
        """

        return (
            self.ensemble.temp * self.nbeads * units.Constants.kb / units.Constants.hbar
        )

    def get_omegan2(self):
        """Returns omegan**2."""

        return self.omegan ** 2

    def get_omegak(self):
        """Gets the normal mode frequencies.

        Returns:
           A list of the normal mode frequencies for the free ring polymer.
           The first element is the centroid frequency (0.0).
        """

        return self.omegan * nmtransform.nm_eva(self.nbeads)

    def get_o_omegak(self):
        """Gets the normal mode frequencies for a open path.

        Returns:
           A list of the normal mode frequencies for the free polymer.
           The first element is the centroid frequency (0.0).
        """

        return self.omegan * nmtransform.o_nm_eva(self.nbeads)

    def get_dynwk(self):
        """Gets the dynamical normal mode frequencies.

        Returns:
           A list of the scaled normal mode frequencies for the free ring polymer.
           The first element is the centroid frequency (0.0).
        """

        return self.omegak / np.sqrt(self.nm_factor)

    def get_prop_pq(self):
        """Gets the exact or Cayley-transformed normal mode propagator matrix.
        The latter allows for longer timestep (nve) and more efficient sampling (nvt).


        Note the special treatment for the centroid normal mode, which is
        propagated using the standard velocity Verlet algorithm as required.
        Note that both the normal mode positions and momenta are propagated
        using this matrix.

        Returns:
           An array of the form (nbeads, 2, 2). Each 2*2 array prop_pq[i,:,:]
           gives the exact propagator or Cayley propagator for the i-th normal mode of the
           ring polymer.
        """

        dt = self.dt
        pqk = np.zeros((self.nbeads, 2, 2), float)
        pqk[0] = np.array([[1, 0], [dt, 1]])

        # Note that the propagator uses mass-scaled momenta.
        if self.propagator == "cayley":
            for b in range(1, self.nbeads):
                sk = np.sqrt(self.nm_factor[b])
                square = (self.omegak[b] * dt / 2) ** 2
                pqk[b, 0, 0] = (1 - square) / (1 + square)
                pqk[b, 1, 1] = (1 - square) / (1 + square)
                pqk[b, 0, 1] = -(4 * square / dt * sk) / (1 + square)
                pqk[b, 1, 0] = dt / sk / (1 + square)
        else:  # exact propagator
            for b in range(1, self.nbeads):
                sk = np.sqrt(self.nm_factor[b])
                dtomegak = self.omegak[b] * dt / sk
                c = np.cos(dtomegak)
                s = np.sin(dtomegak)
                pqk[b, 0, 0] = c
                pqk[b, 1, 1] = c
                pqk[b, 0, 1] = -s * self.omegak[b] * sk
                pqk[b, 1, 0] = s / (self.omegak[b] * sk)
        return pqk

    def get_o_prop_pq(self):
        """Gets the normal mode propagator matrix for the open case.

        Note the special treatment for the centroid normal mode, which is
        propagated using the standard velocity Verlet algorithm as required.
        Note that both the normal mode positions and momenta are propagated
        using this matrix.

        Returns:
           An array of the form (nbeads, 2, 2). Each 2*2 array o_prop_pq[i,:,:]
           gives the exact propagator for the i-th normal mode of the
           ring polymer.
        """

        dt = self.dt
        pqk = np.zeros((self.nbeads, 2, 2), float)
        pqk[0] = np.array([[1, 0], [dt, 1]])

        # Note that the propagator uses mass-scaled momenta.
        if self.propagator == "cayley":
            for b in range(1, self.nbeads):
                sk = np.sqrt(self.o_nm_factor[b])
                square = (self.o_omegak[b] * dt / 2) ** 2
                pqk[b, 0, 0] = (1 - square) / (1 + square)
                pqk[b, 1, 1] = (1 - square) / (1 + square)
                pqk[b, 0, 1] = -(4 * square / dt * sk) / (1 + square)
                pqk[b, 1, 0] = dt / sk / (1 + square)
        else:  # exact propagator
            for b in range(1, self.nbeads):
                sk = np.sqrt(self.o_nm_factor[b])
                dto_omegak = self.o_omegak[b] * dt / sk
                c = np.cos(dto_omegak)
                s = np.sin(dto_omegak)
                pqk[b, 0, 0] = c
                pqk[b, 1, 1] = c
                pqk[b, 0, 1] = -s * self.o_omegak[b] * sk
                pqk[b, 1, 0] = s / (self.o_omegak[b] * sk)
        return pqk

    def get_nmm(self):
        """Returns dynamical mass factors, i.e. the scaling of normal mode
        masses that determine the path dynamics (but not statics)."""

        # also checks that the frequencies and the mode given in init are
        # consistent with the beads and ensemble

        dmf = np.ones(self.nbeads, float)
        if self.mode == "rpmd":
            if len(self.nm_freqs) > 0:
                warning("nm.frequencies will be ignored for RPMD mode.", verbosity.low)
        elif self.mode == "manual":
            if len(self.nm_freqs) != self.nbeads - 1:
                raise ValueError(
                    "Manual path mode requires (nbeads-1) frequencies, one for each internal mode of the path."
                )
            for b in range(1, self.nbeads):
                sk = self.omegak[b] / self.nm_freqs[b - 1]
                dmf[b] = sk ** 2
        elif self.mode == "pa-cmd":
            if len(self.nm_freqs) > 1:
                warning(
                    "Only the first element in nm.frequencies will be considered for PA-CMD mode.",
                    verbosity.low,
                )
            if len(self.nm_freqs) == 0:
                raise ValueError(
                    "PA-CMD mode requires the target frequency of all the internal modes."
                )
            for b in range(1, self.nbeads):
                sk = self.omegak[b] / self.nm_freqs[0]
                info(
                    " ".join(
                        [
                            "NM FACTOR",
                            str(b),
                            str(sk),
                            str(self.omegak[b]),
                            str(self.nm_freqs[0]),
                        ]
                    ),
                    verbosity.medium,
                )
                dmf[b] = sk ** 2
        elif self.mode == "wmax-cmd":
            if len(self.nm_freqs) > 2:
                warning(
                    "Only the first two element in nm.frequencies will be considered for WMAX-CMD mode.",
                    verbosity.low,
                )
            if len(self.nm_freqs) < 2:
                raise ValueError(
                    "WMAX-CMD mode requires [wmax, wtarget]. The normal modes will be scaled such that the first internal mode is at frequency wtarget and all the normal modes coincide at frequency wmax."
                )
            wmax = self.nm_freqs[0]
            wt = self.nm_freqs[1]
            for b in range(1, self.nbeads):
                sk = 1.0 / np.sqrt(
                    (wt) ** 2
                    * (1 + (wmax / self.omegak[1]) ** 2)
                    / (wmax ** 2 + (self.omegak[b]) ** 2)
                )
                dmf[b] = sk ** 2

        return dmf

    # define a function o_nm_factor so we have dynamical masses for the open case
    def get_o_nmm(self):
        """Returns dynamical mass factors, i.e. the scaling of normal mode
        masses that determine the path dynamics (but not statics)."""

        # also checks that the frequencies and the mode given in init are
        # consistent with the beads and ensemble

        dmf = np.ones(self.nbeads, float)
        if self.mode == "rpmd":
            if len(self.nm_freqs) > 0:
                warning("nm.frequencies will be ignored for RPMD mode.", verbosity.low)
        elif self.mode == "manual":
            if len(self.nm_freqs) != self.nbeads - 1:
                raise ValueError(
                    "Manual path mode requires (nbeads-1) frequencies, one for each internal mode of the path."
                )
            for b in range(1, self.nbeads):
                sk = self.o_omegak[b] / self.nm_freqs[b - 1]
                dmf[b] = sk ** 2
        elif self.mode == "pa-cmd":
            if len(self.nm_freqs) > 1:
                warning(
                    "Only the first element in nm.frequencies will be considered for PA-CMD mode.",
                    verbosity.low,
                )
            if len(self.nm_freqs) == 0:
                raise ValueError(
                    "PA-CMD mode requires the target frequency of all the internal modes."
                )
            for b in range(1, self.nbeads):
                sk = self.o_omegak[b] / self.nm_freqs[0]
                info(
                    " ".join(
                        [
                            "NM FACTOR",
                            str(b),
                            str(sk),
                            str(self.o_omegak[b]),
                            str(self.nm_freqs[0]),
                        ]
                    ),
                    verbosity.medium,
                )
                dmf[b] = sk ** 2
        elif self.mode == "wmax-cmd":
            if len(self.nm_freqs) > 2:
                warning(
                    "Only the first two element in nm.frequencies will be considered for WMAX-CMD mode.",
                    verbosity.low,
                )
            if len(self.nm_freqs) < 2:
                raise ValueError(
                    "WMAX-CMD mode requires [wmax, wtarget]. The normal modes will be scaled such that the first internal mode is at frequency wtarget and all the normal modes coincide at frequency wmax."
                )
            wmax = self.nm_freqs[0]
            wt = self.nm_freqs[1]
            for b in range(1, self.nbeads):
                sk = 1.0 / np.sqrt(
                    (wt) ** 2
                    * (1 + (wmax / self.o_omegak[1]) ** 2)
                    / (wmax ** 2 + (self.o_omegak[b]) ** 2)
                )
                dmf[b] = sk ** 2

        return dmf

    def get_dynm3(self):
        """Returns an array with the dynamical masses of individual atoms in the normal modes representation."""

        dm3 = np.zeros(self.beads.m3.shape, float)
        for b in range(self.nbeads):
            dm3[b] = self.beads.m3[b] * self.nm_factor[b]

        # dynamical masses for the open paths
        for j in self.open_paths:
            for a in range(3 * j, 3 * (j + 1)):
                for k in range(1, self.nbeads):
                    dm3[k, a] = self.beads.m3[k, a] * self.o_nm_factor[k]
        return dm3

    def get_vspring_and_fspring_B(self):
        """
        Calculates spring forces and potential for bosons.
        Evaluated using recursion relation from arXiv:1905.090.
        """

        if len(self.bosons) == 0:
            pass
        else:
            (E_k_N, V) = Evaluate_VB(self)

            P = self.nbeads

            F = np.zeros((P, 3 * self.natoms), float)

            for ind, l in enumerate(self.bosons):
                for j in range(P):
                    F[j, 3 * l : 3 * (l + 1)] = Evaluate_dVB(self, E_k_N, V, ind, j)

            return [V[-1], F]

    def get_fspring(self):
        """
        Returns the spring force. Required for numerical propagation in free_babstep().
        For distinguishable particles, simply transform fnm to Cartesian coordinates.
        For bosons, get the second element of vspring_and_fspring_B
        For a mixture of both, calculate separately and combine.
        """

        if len(self.bosons) == 0:
            return self.transform.nm2b(dstrip(self.fspringnm))
        elif len(self.bosons) is self.natoms:
            return self.vspring_and_fspring_B[1]
        else:
            # raise("@NormalModes: Implementing mixtures of B and D")
            f_distinguish = self.transform.nm2b(dstrip(self.fspringnm))
            # zero force for bosons
            for boson in self.bosons:
                f_distinguish[:, 3 * boson : (3 * boson + 3)] = 0.0

            return f_distinguish + self.vspring_and_fspring_B[1]

    def free_babstep(self):
        """
        Numerical propagator in Cartesian coordinates.
        So the propagation is done through a velocity verlet step with a time step that is
        self.nmts smaller than the one for the physical forces.
        All beads of all atoms are propagated in one step.
        Works for both distinguishable particles and bosons. Difference is in fspring.
        """

        if self.nbeads == 1:
            pass
        else:

            # Since the dynamics are done in Cartesian coordinates below (including all modes),
            # I need to revert centroid step done separately in qcstep
            self.qnm[0, :] -= (
                dstrip(self.pnm)[0, :] / dstrip(self.beads.m3)[0] * self.dt
            )

            # Free ring polymer dynamics are done with smaller time step detlat = dt/nmts
            dt = self.dt / dstrip(self.nmts)

            for j in range(0, dstrip(self.nmts)):

                self.beads.p += 0.5 * dt * self.fspring
                self.beads.q += dt * self.beads.p / dstrip(self.beads.m3)
                # The depend machinery will take care of automatically calculating
                # the forces at the updated positions.
                self.beads.p += 0.5 * dt * self.fspring

    def free_qstep(self):
        # !BH!: Should we update the comment here that now the propagator is either exact, NM or numerical, Cartesian?
        """Exact normal mode propagator for the free ring polymer.

        Note that the propagator works in mass scaled coordinates, so that the
        propagator matrix can be determined independently from the particular
        atom masses, and so the same propagator will work for all the atoms in
        the system. All the ring polymers are propagated at the same time by a
        matrix multiplication.

        Also note that the centroid coordinate is propagated in qcstep, so is
        not altered here.
        """

        if self.nbeads == 1:
            pass

        elif self.propagator == "bab":

            if len(self.open_paths) > 0:
                raise (
                    "@Normalmodes : Open path propagator not implemented for bosons. Feel free to implement it if you want to use it :) "
                )

            self.free_babstep()

        else:
            if len(self.bosons) > 0:
                raise (
                    "@Normalmodes : Bosonic forces not compatible right now with the exact or Cayley propagators."
                )

            pq = np.zeros((2, self.natoms * 3), float)
            sm = dstrip(self.beads.sm3)
            prop_pq = dstrip(self.prop_pq)
            o_prop_pq = dstrip(self.o_prop_pq)
            pnm = dstrip(self.pnm) / sm
            qnm = dstrip(self.qnm) * sm

            for k in range(1, self.nbeads):
                pq[0, :] = pnm[k]
                pq[1, :] = qnm[k]
                pq = np.dot(prop_pq[k], pq)
                qnm[k] = pq[1, :]
                pnm[k] = pq[0, :]

            for k in range(1, self.nbeads):
                pq[0, :] = pnm[k]
                pq[1, :] = qnm[k]
                qnm[k] = pq[1, :]
                pnm[k] = pq[0, :]

            # now for open paths we recover the initial conditions (that have not yet been overwritten)
            # and do open path propagation
            pq = np.zeros(2)
            for j in self.open_paths:
                for a in range(3 * j, 3 * (j + 1)):
                    for k in range(1, self.nbeads):
                        pq[0] = self.pnm[k, a] / sm[k, a]
                        pq[1] = self.qnm[k, a] * sm[k, a]
                        pq = np.dot(o_prop_pq[k], pq)
                        qnm[k, a] = pq[1]
                        pnm[k, a] = pq[0]
            self.pnm = pnm * sm
            self.qnm = qnm / sm
            # pq = np.zeros((2,self.natoms*3),float)
            # sm = dstrip(self.beads.sm3)[0]
            # prop_pq = dstrip(self.prop_pq)
            # for k in range(1,self.nbeads):
            #   pq[0,:] = dstrip(self.pnm)[k]/sm
            #   pq[1,:] = dstrip(self.qnm)[k]*sm
            #   pq = np.dot(prop_pq[k],pq)
            #   self.qnm[k] = pq[1,:]/sm
            #   self.pnm[k] = pq[0,:]*sm

    def get_kins(self):
        """Gets the MD kinetic energy for all the normal modes.

        Returns:
           A list of the kinetic energy for each NM.

        """
        # include the partially adiabatic CMD mass scaling
        pnm = dstrip(self.pnm) / dstrip(self.beads.sm3)
        kmd = 0.5 * (pnm ** 2).sum(axis=1) / dstrip(self.nm_factor)

        return kmd

    def get_kin(self):
        """Gets the total MD kinetic energy.

        Note that this does not correspond to the quantum kinetic energy estimate
        for the system.

        Returns:
           The sum of the kinetic energy of each NM in the path.
        """

        return self.kins.sum()

    def get_kstress(self):
        """Calculates the total MD kinetic stress tensor.

        Note that this does not correspond to the quantum kinetic stress tensor
        estimate for the system.

        Returns:
           The sum of the MD kinetic stress tensor contributions from each NM.
        """

        kmd = np.zeros((3, 3), float)
        sm = dstrip(self.beads.sm3[0])
        pnm = dstrip(self.pnm)
        nmf = dstrip(self.nm_factor)

        for b in range(self.nbeads):
            sp = pnm[b] / sm  # mass-scaled momentum of b-th NM

            for i in range(3):
                for j in range(3):
                    # computes the outer product of the p of various normal modes
                    # singling out Cartesian components to build the tensor
                    # also takes care of the possibility of having non-RPMD masses
                    kmd[i, j] += (
                        np.dot(sp[i : 3 * self.natoms : 3], sp[j : 3 * self.natoms : 3])
                        / nmf[b]
                    )

        return kmd

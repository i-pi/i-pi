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
from ipi.utils.array_backend import xp, to_numpy
from ipi.utils.messages import verbosity, warning, info
from ipi.utils.exchange import *

__all__ = ["NormalModes"]


class NormalModes:
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
          Trotter values, or economised ("eco") values fitted to reproduce
          harmonic radii of gyration. Depends on omegan, mode and nm_freqs.
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
        fft_threads=1,
        fft_float32=False,
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
        self.open_paths = np.asarray(open_paths, int)
        if bosons is None:
            bosons = np.zeros(0, int)
        self._bosons = depend_value(name="bosons", value=bosons)
        self._nmts = depend_value(name="nmts", value=nmts)
        self._dt = depend_value(name="dt", value=dt)
        self._mode = depend_value(name="mode", value=mode)
        self._transform_method = depend_value(
            name="transform_method", value=transform_method
        )
        self._propagator = depend_value(name="propagator", value=propagator)
        self._nm_freqs = depend_array(
            name="nm_freqs", value=np.asarray(freqs, float), host=True
        )
        self.fft_threads = fft_threads
        self.fft_float32 = fft_float32

    def copy(self, freqs=None):
        """Creates a new NormalModes object from the original.

        Returns:
           A NormalModes object with the same arrays as the original.
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
            self.fft_threads,
            self.fft_float32,
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

        # storage space for the propagator
        self.pq_buffer = np.zeros((2, self.natoms * 3), float)

        self.bosons = self.resolve_bosons()

        # eco springs are only defined for closed, distinguishable paths
        if self.mode == "eco":
            if len(self.bosons) > 0 or len(self.open_paths) > 0:
                raise ValueError("ECO mode cannot be used with bosons or open paths.")
            if len(self.nm_freqs) != 1:
                raise ValueError(
                    "ECO mode requires one frequency, the maximum physical frequency to be reproduced."
                )

        # stores a reference to the bound beads and ensemble objects
        self.ensemble = ensemble
        dpipe(motion._dt, self._dt)

        # sets up what's necessary to perform nm transformation.
        if self.nbeads == 1:  # classical trajectory! don't waste time doing anything!
            self.transform = nmtransform.nm_noop(nbeads=self.nbeads)
        elif self.transform_method == "fft":
            self.transform = nmtransform.nm_fft(
                nbeads=self.nbeads,
                natoms=self.natoms,
                open_paths=self.open_paths,
                n_threads=self.fft_threads,
                single_precision=self.fft_float32,
            )
        elif self.transform_method == "matrix":
            self.transform = nmtransform.nm_trans(
                nbeads=self.nbeads, open_paths=self.open_paths
            )

        # normal-mode representation of the path, synchronized with the bead
        # representation, which holds the initial positions and momenta
        self._qnm = depend_array(
            name="qnm", value=np.zeros((self.nbeads, 3 * self.natoms), float)
        )
        self._pnm = depend_array(
            name="pnm", value=np.zeros((self.nbeads, 3 * self.natoms), float)
        )
        dsynchro(
            self.beads._q,
            self._qnm,
            lambda: self.transform.nm2b(self.qnm),
            lambda: self.transform.b2nm(self.beads.q),
        )
        dsynchro(
            self.beads._p,
            self._pnm,
            lambda: self.transform.nm2b(self.pnm),
            lambda: self.transform.b2nm(self.beads.p),
        )

        # forces can be converted in nm representation, but here it makes no sense to set up a sync mechanism,
        # as they always get computed in the bead rep
        if self.forces is not None:
            self._fnm = depend_array(
                name="fnm",
                value=np.zeros((self.nbeads, 3 * self.natoms), float),
                func=(lambda: self.transform.b2nm(self.forces.f)),
                dependencies=[self.forces._f],
            )
        else:  # have a fall-back plan when we don't want to initialize a force mechanism, e.g. for ring-polymer initialization

            def no_forces():
                raise ValueError(
                    "Cannot access NM forces when initializing the NM object without providing a force reference!"
                )

            self._fnm = depend_array(
                name="fnm",
                value=np.zeros((self.nbeads, 3 * self.natoms), float),
                func=no_forces,
                dependencies=[],
            )

        # create path-frequencies related properties
        self._omegan = depend_value(
            name="omegan", func=self.get_omegan, dependencies=[self.ensemble._temp]
        )
        self._omegan2 = depend_value(
            name="omegan2", func=self.get_omegan2, dependencies=[self._omegan]
        )
        self._omegak = depend_array(
            name="omegak",
            value=np.zeros(self.beads.nbeads, float),
            func=self.get_omegak,
            dependencies=[self._omegan, self._nm_freqs, self._mode],
        )
        self._omegak2 = depend_array(
            name="omegak2",
            value=np.zeros(self.beads.nbeads, float),
            func=(lambda: self.omegak**2),
            dependencies=[self._omegak],
        )

        # Add o_omegak to calculate the freq in the case of open path
        self._o_omegak = depend_array(
            name="o_omegak",
            value=np.zeros(self.beads.nbeads, float),
            func=self.get_o_omegak,
            dependencies=[self._omegan],
        )

        # sets up "dynamical" masses -- mass-scalings to give the correct RPMD/CMD dynamics
        self._nm_factor = depend_array(
            name="nm_factor",
            value=np.zeros(self.nbeads, float),
            func=self.get_nmm,
            dependencies=[self._nm_freqs, self._mode],
        )
        # add o_nm_factor for the dynamical mass in the case of open paths
        self._o_nm_factor = depend_array(
            name="nmm",
            value=np.zeros(self.nbeads, float),
            func=self.get_o_nmm,
            dependencies=[self._nm_freqs, self._mode],
        )
        self._dynm3 = depend_array(
            name="dynm3",
            value=np.zeros((self.nbeads, 3 * self.natoms), float),
            func=self.get_dynm3,
            dependencies=[self._nm_factor, self.beads._m3],
        )
        self._dynomegak = depend_array(
            name="dynomegak",
            value=np.zeros(self.nbeads, float),
            func=self.get_dynwk,
            dependencies=[self._nm_factor, self._omegak],
        )

        self._dt = depend_value(name="dt", value=1.0)
        dpipe(self.motion._dt, self._dt)

        self._prop_pq = depend_array(
            name="prop_pq",
            value=np.zeros((self.beads.nbeads, 2, 2)),
            func=self.get_prop_pq,
            dependencies=[self._omegak, self._nm_factor, self._dt, self._propagator],
        )

        # mass-scaled propagator
        self._prop_pq_ms = depend_array(
            name="prop_pq_ms",
            value=np.zeros((2, 2, self.beads.nbeads, 3 * self.beads.natoms)),
            func=self.get_prop_pq_ms,
            dependencies=[self._prop_pq, self.beads._m3],
        )

        self._o_prop_pq = depend_array(
            name="o_prop_pq",
            value=np.zeros((self.beads.nbeads, 2, 2)),
            func=self.get_o_prop_pq,
            dependencies=[
                self._o_omegak,
                self._o_nm_factor,
                self._dt,
                self._propagator,
            ],
        )

        # mass-scaled propagator
        self._o_prop_pq_ms = depend_array(
            name="o_prop_pq_ms",
            value=np.zeros((2, 2, self.beads.nbeads, 3 * len(self.open_paths))),
            func=self.get_o_prop_pq_ms,
            dependencies=[self._o_prop_pq, self.beads._m3],
        )

        self.open_paths_coords = np.array(
            [[3 * i, 3 * i + 1, 3 * i + 2] for i in self.open_paths], int
        ).reshape(-1)

        # if the mass matrix is not the RPMD one, the MD kinetic energy can't be
        # obtained in the bead representation because the masses are all mixed up
        self._kins = depend_array(
            name="kins",
            value=np.zeros(self.nbeads, float),
            func=self.get_kins,
            dependencies=[self._pnm, self.beads._sm3, self._nm_factor],
        )
        self._kin = depend_value(
            name="kin", func=self.get_kin, dependencies=[self._kins]
        )
        self._kstress = depend_array(
            name="kstress",
            value=np.zeros((3, 3), float),
            func=self.get_kstress,
            dependencies=[self._pnm, self.beads._sm3, self._nm_factor],
        )

        self._vspring_and_fspring = depend_value(
            name="v_and_fs",
            value=[None, None],
            func=self.get_vspring_and_fspring,
            dependencies=[
                self._qnm,
                self._omegak,
                self._o_omegak,
                self.beads._m3,
                self.beads._q,
            ],
        )

        if len(self.bosons) > 0:
            self.exchange_potential = ExchangePotential(self.nbeads, len(self.bosons))
            self.exchange_potential.bind(self.beads, self.ensemble, self)
            self._vspring_and_fspring.add_dependency(
                self.exchange_potential._vspring_and_fspring
            )
        else:
            self.exchange_potential = None

        # just return split potential and force for ease of access
        self._vspring = depend_value(
            name="vspring",
            value=0.0,
            func=lambda: self.vspring_and_fspring[0],
            dependencies=[self._vspring_and_fspring],
        )

        self._fspring = depend_array(
            name="fspring",
            value=np.zeros((self.nbeads, 3 * self.natoms), float),
            func=lambda: self.vspring_and_fspring[1],
            dependencies=[self._vspring_and_fspring],
        )

        if self.propagator == "cayley" and not bool(
            xp.all(self.dynm3 == self.beads.m3)
        ):
            warning(
                "Cayley propagation has not been tested with dynamical masses "
                "different from the physical masses."
            )

    def resolve_bosons(self):
        if not isinstance(self.bosons, tuple):
            return self.bosons

        bosons_lst, id_mode = self.bosons
        if id_mode == "index":
            bosons_array = bosons_lst.astype(int)
        elif id_mode == "label":
            for latom in bosons_lst:
                if latom not in set(self.beads.names):
                    raise ValueError("Unknown atom label %s for boson" % latom)
            bosons_array = np.asarray(
                [
                    i
                    for i in range(self.beads.natoms)
                    if (self.beads.names[i] in bosons_lst)
                ]
            )
        else:
            raise ValueError(
                "Error resolving boson identifies using unknown method %s" % id_mode
            )

        if len(bosons_array) > 0 and (
            np.min(bosons_array) < 0 or np.max(bosons_array) >= self.beads.natoms
        ):
            raise ValueError("Invalid index for boson, got %s" % str(bosons_array))

        return bosons_array

    def get_omegan(self):
        """Returns the effective vibrational frequency for the interaction
        between replicas.
        """

        return (
            self.ensemble.temp * self.nbeads * units.Constants.kb / units.Constants.hbar
        )

    def get_omegan2(self):
        """Returns omegan**2."""

        return self.omegan**2

    def get_omegak(self):
        """Gets the normal mode frequencies.

        Returns:
           A list of the normal mode frequencies for the free ring polymer.
           The first element is the centroid frequency (0.0).
           For mode="eco", returns economised frequencies optimized to
           reproduce the radii of gyration of harmonic oscillators with
           frequencies up to nm_freqs[0].
        """

        if self.mode == "eco":
            if len(self.nm_freqs) != 1:
                raise ValueError(
                    "ECO mode requires one frequency, the maximum physical frequency to be reproduced."
                )
            # xmax = beta*hbar*omega_max, given that omegan = nbeads/(beta*hbar)
            xmax = self.nm_freqs[0] * self.nbeads / self.omegan
            # seeds the fit with the previous solution (the stale value held in
            # the depend array), which speeds up re-fits when temperature changes
            y0 = (
                to_numpy(self._omegak._value[1 : self.nbeads // 2 + 1])
                * self.nbeads
                / self.omegan
            )
            if not np.all(y0 > 0):  # zeros before the first evaluation
                y0 = None
            return xp.asarray(self.omegan * nmtransform.eco_eva(self.nbeads, xmax, y0))
        return xp.asarray(self.omegan * nmtransform.nm_eva(self.nbeads))

    def get_o_omegak(self):
        """Gets the normal mode frequencies for a open path.

        Returns:
           A list of the normal mode frequencies for the free polymer.
           The first element is the centroid frequency (0.0).
        """

        return xp.asarray(self.omegan * nmtransform.o_nm_eva(self.nbeads))

    def get_dynwk(self):
        """Gets the dynamical normal mode frequencies.

        Returns:
           A list of the scaled normal mode frequencies for the free ring polymer.
           The first element is the centroid frequency (0.0).
        """

        return self.omegak / xp.sqrt(self.nm_factor)

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

        return self._propagator_matrices(self.omegak, self.nm_factor)

    def _propagator_matrices(self, omegak, nm_factor):
        """Exact or Cayley propagators for normal modes with frequencies
        `omegak` and dynamical mass factors `nm_factor`, as an array of shape
        (nbeads, 2, 2). The centroid is propagated with velocity Verlet."""

        dt = self.dt
        # Note that the propagator uses mass-scaled momenta.
        sk = xp.sqrt(nm_factor[1:])
        wk = omegak[1:]
        if self.propagator == "cayley":
            square = (wk * dt / 2) ** 2
            diag = (1 - square) / (1 + square)
            p_from_q = -(4 * square / dt * sk) / (1 + square)
            q_from_p = dt / sk / (1 + square)
        else:  # exact propagator
            dtomegak = wk * dt / sk
            diag = xp.cos(dtomegak)
            p_from_q = -xp.sin(dtomegak) * wk * sk
            q_from_p = xp.sin(dtomegak) / (wk * sk)
        pqk = xp.stack(
            [xp.stack([diag, p_from_q], axis=-1), xp.stack([q_from_p, diag], axis=-1)],
            axis=-2,
        )
        centroid = xp.asarray([[[1.0, 0.0], [dt, 1.0]]])
        return xp.concat([centroid, pqk], axis=0)

    def get_prop_pq_ms(self):
        """Combines the propagator with the atom masses to make a single array
        that can be multiplied to propagate the normal modes dynamics"""

        pq_ms = xp.zeros((2, 2, self.nbeads, self.natoms * 3))
        pq_ms[...] = xp.permute_dims(self.prop_pq, (1, 2, 0))[:, :, :, None]
        pq_ms[0, 1] *= self.beads.m3
        pq_ms[1, 0] /= self.beads.m3

        return pq_ms

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

        return self._propagator_matrices(self.o_omegak, self.o_nm_factor)

    def get_o_prop_pq_ms(self):
        """Combines the propagator with the atom masses to make a single array
        that can be multiplied to propagate the normal modes dynamics. Open path version
        """

        pq_ms = xp.zeros((2, 2, self.nbeads, len(self.open_paths_coords)))
        pq_ms[...] = xp.permute_dims(self.o_prop_pq, (1, 2, 0))[:, :, :, None]
        pq_ms[0, 1] *= self.beads.m3[:, self.open_paths_coords]
        pq_ms[1, 0] /= self.beads.m3[:, self.open_paths_coords]

        return pq_ms

    def get_nmm(self):
        """Returns dynamical mass factors, i.e. the scaling of normal mode
        masses that determine the path dynamics (but not statics)."""

        # also checks that the frequencies and the mode given in init are
        # consistent with the beads and ensemble

        omegak = to_numpy(self.omegak)
        dmf = np.ones(self.nbeads, float)
        if self.mode == "rpmd":
            if len(self.nm_freqs) > 0:
                warning("nm.frequencies will be ignored for RPMD mode.", verbosity.low)
        elif self.mode == "eco":
            # eco changes the spring frequencies (statics), not the dynamical masses
            pass
        elif self.mode == "manual":
            if len(self.nm_freqs) != self.nbeads - 1:
                raise ValueError(
                    "Manual path mode requires (nbeads-1) frequencies, one for each internal mode of the path."
                )
            for b in range(1, self.nbeads):
                sk = omegak[b] / self.nm_freqs[b - 1]
                dmf[b] = sk**2
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
                sk = omegak[b] / self.nm_freqs[0]
                info(
                    " ".join(
                        [
                            "NM FACTOR",
                            str(b),
                            str(sk),
                            str(omegak[b]),
                            str(self.nm_freqs[0]),
                        ]
                    ),
                    verbosity.medium,
                )
                dmf[b] = sk**2
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
                    * (1 + (wmax / omegak[1]) ** 2)
                    / (wmax**2 + (omegak[b]) ** 2)
                )
                dmf[b] = sk**2

        return xp.asarray(dmf)

    # define a function o_nm_factor so we have dynamical masses for the open case
    def get_o_nmm(self):
        """Returns dynamical mass factors, i.e. the scaling of normal mode
        masses that determine the path dynamics (but not statics)."""

        # also checks that the frequencies and the mode given in init are
        # consistent with the beads and ensemble

        omegak = to_numpy(self.o_omegak)
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
                sk = omegak[b] / self.nm_freqs[b - 1]
                dmf[b] = sk**2
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
                sk = omegak[b] / self.nm_freqs[0]
                info(
                    " ".join(
                        [
                            "NM FACTOR",
                            str(b),
                            str(sk),
                            str(omegak[b]),
                            str(self.nm_freqs[0]),
                        ]
                    ),
                    verbosity.medium,
                )
                dmf[b] = sk**2
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
                    * (1 + (wmax / omegak[1]) ** 2)
                    / (wmax**2 + (omegak[b]) ** 2)
                )
                dmf[b] = sk**2

        return xp.asarray(dmf)

    def get_dynm3(self):
        """Returns an array with the dynamical masses of individual atoms in the normal modes representation."""

        m3 = self.beads.m3
        dm3 = m3 * self.nm_factor[:, None]

        # dynamical masses for the open paths
        if len(self.open_paths) > 0:
            opc = self.open_paths_coords
            dm3[1:, opc] = m3[1:, opc] * self.o_nm_factor[1:, None]
        return dm3

    def get_vspring_and_fspring(self):
        """Returns the total spring energy and spring force."""

        # classical simulation - do nothing!
        vspring, fspring = 0.0, xp.zeros_like(self.qnm)
        if self.nbeads == 1:
            return vspring, fspring

        if len(self.bosons) < self.natoms:
            # computes harmonic potential in normal-modes representation
            qnm = self.qnm
            sqnm = qnm * self.beads.sm3
            q2 = (sqnm**2).sum(axis=1)

            vspring += (self.omegak2 * q2).sum()
            # spring forces (this is only the closed path version -
            # open path correction is applied in the propagator)
            # TODO make it happen here, it'd be probably cleaner
            fspringnm = -self.beads.m3 * self.omegak[:, None] ** 2 * qnm

            # correction for open paths
            for j in self.open_paths:
                vspring += xp.sum(
                    self.beads.m[j]
                    * (self.o_omegak**2 - self.omegak**2)
                    * (
                        self.qnm[:, 3 * j] ** 2
                        + self.qnm[:, 3 * j + 1] ** 2
                        + self.qnm[:, 3 * j + 2] ** 2
                    )
                )

                # overrides open path forces
                fspringnm[:, 3 * j : 3 * (j + 1)] = (
                    -self.beads.m3[:, 3 * j : 3 * (j + 1)]
                    * self.o_omegak[:, None] ** 2
                    * qnm[:, 3 * j : 3 * (j + 1)]
                )

            # correction for bosons
            if len(self.bosons) > 0:
                for j in self.bosons:
                    vspring -= xp.sum(
                        self.beads.m[j]
                        * self.omegak**2
                        * (
                            self.qnm[:, 3 * j] ** 2
                            + self.qnm[:, 3 * j + 1] ** 2
                            + self.qnm[:, 3 * j + 2] ** 2
                        )
                    )

            vspring *= 0.5
            fspring += self.transform.nm2b(fspringnm)

        if len(self.bosons) > 0:
            vspring_b, fspring_b = self.exchange_potential.vspring_and_fspring

            vspring += vspring_b
            # overwrites the spring forces for the bosonic particles
            fspring = fspring.reshape(self.nbeads, self.natoms, 3)
            fspring[:, self.bosons, :] = fspring_b
            fspring = fspring.reshape(self.nbeads, 3 * self.natoms)

        return vspring, fspring

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
            self._qnm[0, :] -= self.pnm[0, :] / self.beads.m3[0] * self.dt

            # Free ring polymer dynamics are done with smaller time step detlat = dt/nmts
            dt = self.dt / self.nmts

            for j in range(0, self.nmts):
                self.beads._p += 0.5 * dt * self.fspring
                self.beads._q += dt * self.beads.p / self.beads.m3
                # The depend machinery will take care of automatically calculating
                # the forces at the updated positions.
                self.beads._p += 0.5 * dt * self.fspring

    def free_qstep(self):
        """Position propagator for the free ring polymer.

        The basic propagator is based on a normal mode transformation, and
        works for closed and open paths. Note that the propagator works in mass
        scaled coordinates, so that the propagator matrix can be determined
        independently from the particular atom masses, and so the same propagator
        will work for all the atoms in the system. All the ring polymers are
        propagated at the same time by a matrix multiplication.
        A special Cartesian propagator is used for bosonic exchange.

        Also note that the centroid coordinate is propagated in qcstep, so is
        not altered here.
        """

        if self.nbeads == 1:
            pass

        elif self.propagator == "bab":
            if len(self.open_paths) > 0:
                raise NotImplementedError(
                    """@Normalmodes : Open path propagator not implemented for bosons. 
                    Feel free to implement it if you want to use it :) """
                )

            self.free_babstep()

        else:
            if len(self.bosons) > 0:
                raise NotImplementedError(
                    "@Normalmodes : Bosonic forces not compatible right now with the exact or Cayley propagators."
                )

            """
            if len(self.open_paths) > 0:
                # if there are open paths, make a copy to preserve the original values
                pnm = pnm.copy()
                qnm = qnm.copy()

            # uses the buffer to apply the propagator in one go
            pq_buffer = self.pq_buffer

            for k in range(1, self.nbeads):
                # goes in mass-scaled coordinates and then applies
                pq_buffer[0] = pnm[k] / sm[k]
                pq_buffer[1] = qnm[k] * sm[k]
                pnm[k], qnm[k] = noddot(prop_pq[k], pq_buffer)                
            res1 = self.pnm[1:] * self.prop_pq_ms[0,0,1:] + self.qnm[1:] * self.prop_pq_ms[0,1,1:]
            res2 = self.pnm[1:] * self.prop_pq_ms[1,0,1:] + self.qnm[1:] * self.prop_pq_ms[1,1,1:]
            
            print("check, ", np.linalg.norm( pnm[1:] * sm[1:] - res1), 
            np.linalg.norm( qnm[1:] / sm[1:] - res2))
            # back to non-scaled coordinates, and update the actual arrays
            self.pnm[1:] = pnm[1:] * sm[1:]
            self.qnm[1:] = qnm[1:] / sm[1:]

            """

            # detach arrays
            pnm = self.pnm
            qnm = self.qnm

            prop_pq_ms = self.prop_pq_ms
            new_pnm = pnm[1:] * prop_pq_ms[0, 0, 1:] + qnm[1:] * prop_pq_ms[0, 1, 1:]
            new_qnm = pnm[1:] * prop_pq_ms[1, 0, 1:] + qnm[1:] * prop_pq_ms[1, 1, 1:]

            # now for open paths we recover the initial conditions (that have not yet been overwritten)
            # and do open path propagation. NB this will be slow, will probably need some optimization
            # to run large calculations with this
            if len(self.open_paths) > 0:
                o_prop_pq_ms = self.o_prop_pq_ms
                opc = self.open_paths_coords
                new_pnm[:, opc] = (
                    pnm[1:, opc] * o_prop_pq_ms[0, 0, 1:]
                    + qnm[1:, opc] * o_prop_pq_ms[0, 1, 1:]
                )
                new_qnm[:, opc] = (
                    pnm[1:, opc] * o_prop_pq_ms[1, 0, 1:]
                    + qnm[1:, opc] * o_prop_pq_ms[1, 1, 1:]
                )
                """
                pq = np.zeros(2)
                for j in self.open_paths:
                    for a in range(3 * j, 3 * (j + 1)):
                        for k in range(1, self.nbeads):
                            pq[0] = self.pnm[k, a] / sm[k, a]
                            pq[1] = self.qnm[k, a] * sm[k, a]
                            pq = np.dot(o_prop_pq[k], pq)
                            new_pnm[k, a] = pq[0] * sm[k,a]
                            new_qnm[k, a] = pq[1] / sm[k,a]
                """
            # updates the "real" vectors
            self._pnm[1:] = new_pnm
            self._qnm[1:] = new_qnm

    def get_kins(self):
        """Gets the MD kinetic energy for all the normal modes.

        Returns:
           A list of the kinetic energy for each NM.

        """
        # include the partially adiabatic CMD mass scaling
        pnm = self.pnm / self.beads.sm3
        kmd = 0.5 * (pnm**2).sum(axis=1) / self.nm_factor

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

        # computes mass and nm_scaling-scaled momenta
        spnm = (self.pnm / xp.sqrt(self.nm_factor[:, None])) / self.beads.sm3
        spnm = spnm.reshape(-1, 3)  # reshape for quick outer product
        return spnm.T @ spnm


dproperties(
    NormalModes,
    [
        "nmts",
        "dt",
        "mode",
        "transform_method",
        "propagator",
        "nm_freqs",
        "bosons",
        "qnm",
        "pnm",
        "fnm",
        "omegan",
        "omegan2",
        "omegak",
        "omegak2",
        "o_omegak",
        "nm_factor",
        "o_nm_factor",
        "dynm3",
        "dynomegak",
        "prop_pq",
        "prop_pq_ms",
        "o_prop_pq",
        "o_prop_pq_ms",
        "kins",
        "kin",
        "kstress",
        "exchange",
        "vspring",
        "vspring_and_fspring_bosons",
        "vspring_and_fspring_distinguishables",
        "vspring_and_fspring",
        "fspring",
        "fspringnm",
    ],
)

"""Contains the classes that evaluate forces on PI beads.

This contains both the class that gets the force acting on the beads,
and the class to compute individual components -- in case one wants to
use multiple force providers to get e.g. bonded and non-bonded interactions.
It is an extra layer between the dynamics (that only cares about TOTAL force)
and the driver (that only cares about a single bead).
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import time
import threading
from copy import deepcopy
from functools import partial

import numpy as np

from ipi.utils.softexit import softexit
from ipi.utils.messages import verbosity, warning, info
from ipi.utils import jit
from ipi.utils.depend import *
from ipi.utils.array_backend import xp, to_numpy
from ipi.utils.nmtransform import nm_rescale
from ipi.engine.beads import Beads
from ipi.engine.cell import Cell

__all__ = ["Forces", "ForceComponent"]


fbuid = 0


def _stack(arrays):
    """Stacks the arrays of the beads, as a view of a single one."""

    return arrays[0][None] if len(arrays) == 1 else xp.stack(arrays)


def _in_mts_level(ff, level):
    """Whether the force component `ff` contributes to the MTS level `level`."""

    return len(ff.mts_weights) > level and ff.mts_weights[level] != 0


def _scale(component, level=None):
    """The factor of a force component in the sums over the components: its
    weight, times its MTS weight at `level`, or the sum of its MTS weights if
    `level` is None. None if the component does not contribute."""

    if component.weight == 0:
        return None
    if level is None:
        return component.weight * component.mts_weights.sum()
    if _in_mts_level(component, level):
        return component.weight * component.mts_weights[level]
    return None


def _weighted_sum(mforces, mrpc, name, shape, level=None):
    """Sums the arrays `name` (f, pots or virs) of the force components that
    contribute (see `_scale`), expanded to the beads of the system and
    scaled. Returns zeros of shape `shape` if no component contributes."""

    total = None
    for component, rpc in zip(mforces, mrpc):
        scale = _scale(component, level)
        if scale is None:
            continue
        value = getattr(component, name)
        if not rpc.noop:
            value = rpc.b2tob1(value)
        # a new array, also for a single component
        value = scale * value
        total = value if total is None else total + value
    return xp.zeros(shape) if total is None else total


def _check_virial_symmetry(vir, ff):
    """Warns if a forcefield returns a full virial that is not symmetric,
    as only its upper triangular part is used."""

    ff.virial_checked = True
    vir = to_numpy(vir)
    lower = np.array([vir[1, 0], vir[2, 0], vir[2, 1]])
    upper = np.array([vir[0, 1], vir[0, 2], vir[1, 2]])
    if lower.any() and np.abs(lower - upper).max() > 1e-4 * np.abs(vir).max():
        warning(
            f"Forcefield '{ff.name}' returned a non-symmetric virial: "
            "only its upper triangular part is used",
            verbosity.low,
        )


class ForceBead:
    """Base force helper class.

    This is the object that computes forces for a single bead. This is the last
    layer before calling a forcefield object.

    Attributes:
       atoms: An Atoms object containing all the atom positions.
       cell: A Cell object containing the system box.
       ff: A forcefield object which can calculate the potential, virial
          and forces given an unit cell and atom positions of one replica
          of the system.
       uid: A unique id number identifying each of the different bead's
          forcefields.
       request: A dictionary containing information about the currently
          running job.
       _threadlock: Python handle used to lock the thread used to run the
          communication with the client code.
       _getallcount: An integer giving how many times the getall function has
          been called.

    Depend objects:
       ufvx: A list of the form [pot, f, vir, extra]. These quantities are calculated
          all at one time by the driver, so are collected together. Each separate
          object is then taken from the list. Taken from the results of all the
          beads of the force component.
       extra: A string containing some formatted output returned by the client. Depends on ufvx.
       pot: A float giving the potential energy of the system. Depends on ufvx.
       f: An array containing all the components of the force. Depends on ufvx.
       vir: An array containing the components of the virial tensor in upper
          triangular form, not divided by the volume. Depends on ufvx.
       request: a handle to the request that has been filed by the FF object
    """

    def __init__(self, component, index):
        """Initialises ForceBead.

        Args:
           component: The force component the bead belongs to, which collects
              the results of the requests of all its beads.
           index: The index of the bead in the component.
        """

        self.component = component
        # ufvx is a list [ u, f, vir, extra ] which holds the results of the force
        # calculation, taken from those of all the beads (`ForceComponent.results`)
        self.ufvx_ = depend_value(
            name="ufvx",
            func=partial(component.bead_results, index),
            dependencies=[component.results_],
        )
        self._threadlock = threading.Lock()
        self.request = None
        self._getallcount = 0

    def bind(self, atoms, cell, ff, output_maker):
        """Binds atoms, cell and a forcefield template to the ForceBead object.

        Args:
           atoms: The Atoms object from which the atom positions are taken.
           cell: The Cell object from which the system box is taken.
           ff: A forcefield object which can calculate the potential, virial
              and forces given an unit cell and atom positions of one replica
              of the system.
        """

        global fbuid  # assign a unique identifier to each forcebead object
        with self._threadlock:
            self.uid = fbuid
            fbuid += 1

        # stores a reference to the atoms and cell we are computing forces for
        self.atoms = atoms
        self.cell = cell
        self.ff = ff

        # potential and virial are to be extracted very simply from ufvx
        self.pot_ = depend_value(
            name="pot", func=self.get_pot, dependencies=[self.ufvx_]
        )

        self.vir_ = depend_array(
            name="vir",
            value=np.zeros((3, 3), float),
            func=self.get_vir,
            dependencies=[self.ufvx_],
        )

        self.f_ = depend_array(
            name="f",
            value=np.zeros(atoms.natoms * 3, float),
            func=self.get_f,
            dependencies=[self.ufvx_],
        )

        self.extra_ = depend_value(
            name="extra", func=self.get_extra, dependencies=[self.ufvx_]
        )

    def queue(self):
        """Sends the job to the interface queue directly.

        Allows the ForceBead object to ask for the ufvx list of each replica
        directly without going through the get_all function. This allows
        all the jobs to be sent at once, allowing them to be parallelized.

        Returns True if a calculation is running
        """

        # the compilation of a step (--jit) may be interrupted here, to compute
        # the forces eagerly (see ipi.engine.jit)
        jit.split(self)

        # only dispatches a request if the results of the beads are out of date
        with self._threadlock:
            is_tainted = self.component.results_.tainted()
            if self.request is None and is_tainted:
                self.request = self.ff.queue(self.atoms, self.cell, reqid=self.uid)
        return is_tainted

    def get_all(self):
        """Driver routine.

        Sends the atoms and cell to the forcefield, if they have not been sent,
        requesting that it calculates the potential, forces and virial tensor.
        This then waits until the driver is finished, and then returns the
        results of the request.

        Returns:
           A list of the form [potential, force, virial, extra].
        """

        # because we thread over many systems and outputs, we might get called
        # more than once. keep track of how many times we are called so we
        # can make sure to wait until the last call has returned before we release
        with self._threadlock:
            self._getallcount += 1

        # this is converting the distribution library requests into [ u, f, v ]  lists
        if self.request is None:
            self.queue()

        # sleeps until the request has been evaluated
        request = self.request
        latency = self.ff.latency
        while request["status"] != "Done":
            if request["status"] == "Exit" or softexit.triggered:
                # now, this is tricky. we are stuck here and we cannot return meaningful results.
                # if we return, we may as well output wrong numbers, or mess up things.
                # so we can only call soft-exit and wait until that is done. then kill the thread
                # we are in.
                softexit.trigger(
                    status="bad",
                    message=" @ FORCES : cannot return so will die off here",
                )
            request._event_done.wait(timeout=1.0)

        # print diagnostics about the elapsed time
        if verbosity.debug:
            info(
                " @forces: FF %s evaluated in %f (queue) and %f (dispatched) sec."
                % (
                    self.ff.name,
                    request["t_finished"] - request["t_queued"],
                    request["t_finished"] - request["t_dispatched"],
                )
            )

        # data has been collected, so the request can be released and a slot
        # freed up for new calculations
        result = request["result"]

        # reduce the reservation count (and wait for all calls to return)
        with self._threadlock:
            self._getallcount -= 1

        # releases just once, but wait for all requests to be complete
        if self._getallcount == 0:
            self.ff.release(request)
            self.request = None
        else:
            while self._getallcount > 0:
                time.sleep(latency)

        return result

    def get_pot(self):
        """Calls get_all routine of forcefield to update the potential.

        Returns:
           Potential energy.
        """

        return self.ufvx[0]

    def get_f(self):
        """Calls get_all routine of forcefield to update the force.

        Returns:
           An array containing all the components of the force.
        """

        return xp.asarray(self.ufvx[1])

    def get_vir(self):
        """Calls get_all routine of forcefield to update the virial.

        Returns:
           An array containing the virial in upper triangular form, not divided
           by the volume.
        """

        vir = self.ufvx[2]
        if not self.ff.virial_checked:
            _check_virial_symmetry(vir, self.ff)
        # upper triangular copy, leaving the forcefield result untouched
        return xp.triu(vir)

    def get_extra(self):
        """Calls get_all routine of forcefield to update the extras string.

        Returns:
           A string containing all formatted additional output that the
           client might have produced.
        """

        return self.ufvx[3]


dproperties(ForceBead, ["ufvx", "pot", "vir", "f", "extra"])


class ForceComponent:
    """Computes one component (e.g. bonded interactions) of the force.

    Deals with splitting the bead representation into
    separate replicas, and collecting the data from each replica.

    Attributes:
       natoms: An integer giving the number of atoms.
       nbeads: An integer giving the number of beads.
       name: The name of the forcefield.
       _forces: A list of the forcefield objects for all the replicas.
       weight: A float that will be used to weight the contribution of this
          forcefield to the total force.
       mts_weights: A list of floats that will be used to weight the
          contribution of this forcefield at each level of a MTS scheme
       ffield: A model to be used to create the forcefield objects for all
          the replicas of the system.

    Depend objects:
       results: The results of the requests of all the replicas, as [pots, f,
          virs, extras] with the arrays stacked, from which the ufvx list of
          each replica is taken. Depends on the positions and on the cell.
       f: An array containing the components of the force. Taken from results.
       pots: A list containing the potential energy for each system replica.
          Taken from results.
       virs: A list containing the virial tensor for each system replica.
          Taken from results.
       pot: The sum of the potential energy of the replicas.
       vir: The sum of the virial tensor of the replicas.
       extras: The extra quantities returned by the forcefield, as a dictionary
          that maps each key to the list of its values for the replicas, or, for
          the keys in `interpolate_extras`, to a numpy array of shape
          (nbeads, ...). Depends on each replica's ufvx list.
    """

    def __init__(
        self,
        ffield,
        nbeads=0,
        weight=1.0,
        name="",
        mts_weights=None,
        interpolate_extras=None,
        epsilon=-0.001,
    ):
        """Initializes ForceComponent

        Args:
           ffield: A model to be used to create the forcefield objects for all
              the replicas of the system.
           nbeads: The number of replicas.
           weight: A relative weight to be given to the values obtained with this
              forcefield. When the contribution of all the forcefields is
              combined to give a total force, the contribution of this forcefield
              will be weighted by this factor. The combination is a weighted sum.
           name: The name of the forcefield.
           mts_weights: Weight of forcefield at each mts level. Defaults to [1],
              i.e. a single (outer) level with unit weight, as in the input.
           interpolate_extras: A list of properties that should be treated as physical quantities,
              converted to numpy arrays and treated with ring polymer contraction. If
              different force components have this field, they will also be summed with
              the respective weight like a forces object.
        """

        self.ffield = ffield
        self.name = name
        self.nbeads = nbeads
        self.weight_ = depend_value(name="weight", value=weight)
        if mts_weights is None:
            self.mts_weights = np.ones(1)
        else:
            self.mts_weights = np.asarray(mts_weights)
        if interpolate_extras is None:
            self.interpolate_extras = []
        else:
            self.interpolate_extras = interpolate_extras
        self.epsilon = epsilon

    def bind(self, beads, cell, fflist, output_maker):
        """Binds beads, cell and force to the forcefield.

        Takes the beads, cell objects and makes them members of the forcefield.
        Also takes the force object and copies it once for each replica of the
        system, then binds each replica to one of the copies so that the force
        calculation can be parallelized. Creates the objects that will
        hold the data that the driver returns and the dependency network.

        Args:
           beads: Beads object from which the bead positions are taken.
           cell: Cell object from which the system box is taken.
           fflist: A list of forcefield objects to use to calculate the potential,
              forces and virial for each replica.
        """

        # stores a copy of the number of atoms and of beads
        self.natoms = beads.natoms
        if self.nbeads != beads.nbeads:
            raise ValueError(
                "Binding together a Beads and a ForceBeads objects with different numbers of beads"
            )

        # creates an array of force objects, which are bound to the beads
        # and the cell
        if self.ffield not in fflist:
            raise ValueError(
                "Force component name '"
                + self.ffield
                + "' is not in the forcefields list"
            )

        self.ff = fflist[self.ffield]

        # the results of all the beads are collected at once, with their arrays
        # stacked, and those of each bead are taken from them
        self.results_ = depend_value(
            name="results", func=self.results_gather, dependencies=[beads.q_, cell.h_]
        )

        self._forces = []
        self.beads = beads
        for b in range(self.nbeads):
            new_force = ForceBead(self, b)
            new_force.bind(beads[b], cell, self.ff, output_maker=output_maker)
            self._forces.append(new_force)

        # f is a big array which assembles the forces on individual beads
        self.f_ = depend_array(
            name="f",
            value=np.zeros((self.nbeads, 3 * self.natoms)),
            func=self.f_gather,
            dependencies=[self.results_],
        )

        # collection of pots, virs and extras from individual beads
        self.pots_ = depend_array(
            name="pots",
            value=np.zeros(self.nbeads, float),
            func=self.pot_gather,
            dependencies=[self.results_],
        )
        self.virs_ = depend_array(
            name="virs",
            value=np.zeros((self.nbeads, 3, 3), float),
            func=self.vir_gather,
            dependencies=[self.results_],
        )
        self.extras_ = depend_value(
            name="extras",
            value={},
            func=self.extra_gather,
            dependencies=[self._forces[b].extra_ for b in range(self.nbeads)],
        )

        # total potential and total virial
        self.pot_ = depend_value(
            name="pot", func=(lambda: self.pots.sum()), dependencies=[self.pots_]
        )
        self.vir_ = depend_array(
            name="vir",
            func=self.get_vir,
            value=np.zeros((3, 3)),
            dependencies=[self.virs_],
        )

    def results_gather(self):
        """Collects the results of the requests of all the beads, as [pots, f,
        virs, extras], with the arrays stacked."""

        self.queue()
        results = [b.get_all() for b in self._forces]
        return [
            _stack([xp.asarray(r[0]) for r in results]),
            _stack([xp.asarray(r[1]) for r in results]),
            _stack([xp.asarray(r[2]) for r in results]),
            [r[3] for r in results],
        ]

    def bead_results(self, index):
        """Returns the results of the bead `index`, as [pot, f, vir, extra]."""

        pots, f, virs, extras = self.results
        return [pots[index], f[index], virs[index], extras[index]]

    def queue(self):
        """Submits all the required force calculations to the interface.

        Returns True if calculations are running, False if nothing is tainted
        """

        # this should be called in functions which access u,v,f for ALL the beads,
        # before accessing them. it is basically pre-queueing so that the
        # distributed-computing magic can work

        is_tainted = False
        for b in range(self.nbeads):
            is_tainted = self._forces[b].queue() or is_tainted
        return is_tainted

    def pot_gather(self):
        """Obtains the potential energy for each replica.

        Returns:
           A list of the potential energy of each replica of the system.
        """

        return self.results[0]

    def extra_gather(self):
        """Obtains the extra string information for each replica.

        Returns:
           A list of the extra dictionaries for each replica of the system.
        """

        self.queue()

        # converts a list of dictionaries to a dictionary of lists
        fc_extra = {}
        for e in self._forces[0].extra:
            fc_extra[e] = []

        for b in self._forces:
            for e in b.extra:
                if not e in fc_extra:
                    raise KeyError(
                        "Extras mismatch between beads in the same force component, key: "
                        + e
                    )
                fc_extra[e].append(b.extra[e])

        # interpolate_extras should be numerical, thus can be converted to numpy arrays.
        # we enforce the type and numpy will raise an error if not.
        for e in self.interpolate_extras:
            try:
                fc_extra[e] = np.asarray(fc_extra[e], dtype=float)
            except KeyError:
                raise KeyError(
                    "interpolate_extras required "
                    + e
                    + " to promote, but was not found among extras "
                    + str(list(fc_extra.keys()))
                )
            except:
                raise Exception(
                    "interpolate_extras has to be numerical to be treated as a physical quantity. It is not -- check the quantity that is being passed."
                )
        return fc_extra

    def vir_gather(self):
        """Obtains the virial for each replica.

        Returns:
           A list of the virial of each replica of the system.
        """

        virs = self.results[2]
        if not self.ff.virial_checked:
            _check_virial_symmetry(virs[0], self.ff)
        return xp.triu(virs)

    def f_gather(self):
        """Obtains the force vector for each replica.

        Returns:
           An array with all the components of the force. Row i gives the force
           array for replica i of the system.
        """

        return self.results[1]

    def get_vir(self):
        """Sums the virial of each replica.

        Not the actual system virial, as it has not been divided by either the
        number of beads or the cell volume.

        Returns:
            Virial sum.
        """

        return self.virs.sum(axis=0)


dproperties(
    ForceComponent, ["weight", "f", "pots", "pot", "virs", "vir", "extras", "results"]
)


class ScaledForceComponent:
    """A force component that is another one (`baseforce`) multiplied by
    `scaling`: its energies, forces, virials and the extras that are
    interpolated as physical quantities are those of the base component times
    the scaling, e.g. for the bias of Hamiltonian replica exchange.

    Depend objects:
       scaling: The scaling factor. The base component is not evaluated when
          it is zero.
       f, pots, virs, extras: The scaled quantities of the base component.
       pot, vir: Their sums over the replicas.
    """

    def __init__(self, baseforce, scaling=1):
        self.bf = baseforce
        self.name = baseforce.name
        self.nbeads = baseforce.nbeads
        self.ffield = baseforce.ffield
        self.scaling_ = depend_value(name="scaling", value=scaling)
        self.f_ = depend_array(
            name="f",
            func=lambda: (
                self.scaling * self.bf.f
                if self.scaling != 0
                else xp.zeros((self.bf.nbeads, 3 * self.bf.natoms))
            ),
            value=np.zeros((self.bf.nbeads, 3 * self.bf.natoms)),
            dependencies=[self.bf.f_, self.scaling_],
        )
        self.pots_ = depend_array(
            name="pots",
            func=self.get_pots,
            value=np.zeros(self.bf.nbeads),
            dependencies=[self.bf.pots_, self.scaling_],
        )
        self.virs_ = depend_array(
            name="virs",
            func=lambda: self.scaling * self.bf.virs,
            value=np.zeros((self.bf.nbeads, 3, 3)),
            dependencies=[self.bf.virs_, self.scaling_],
        )
        self.extras_ = depend_value(
            name="extras",
            func=self.get_extras,
            value={},
            dependencies=[self.bf.extras_, self.scaling_],
        )

        # total potential and total virial
        self.pot_ = depend_value(
            name="pot", func=(lambda: self.pots.sum()), dependencies=[self.pots_]
        )
        self.vir_ = depend_array(
            name="vir",
            func=self.get_vir,
            value=np.zeros((3, 3)),
            dependencies=[self.virs_],
        )

        # pipes weight from the force, since the scaling is applied on top of that
        self.weight_ = depend_value(name="weight", value=0)
        dpipe(self.bf.weight_, self.weight_)
        self.mts_weights = self.bf.mts_weights
        self.interpolate_extras = self.bf.interpolate_extras

    def get_pots(self):
        return (
            self.scaling * self.bf.pots
            if self.scaling != 0
            else xp.zeros(self.bf.nbeads)
        )

    def get_extras(self):
        """Extras of the base component, scaling those that are interpolated
        as physical quantities."""

        # the extras are numpy arrays, and the scaling may be a backend scalar
        scaling = float(self.scaling)
        if scaling == 0:
            return {}
        extras = dict(self.bf.extras)
        for e in self.interpolate_extras:
            extras[e] = scaling * extras[e]
        return extras

    def get_vir(self):
        """Sums the virial of each replica.

        Not the actual system virial, as it has not been divided by either the
        number of beads or the cell volume.

        Returns:
            Virial sum.
        """

        return self.virs.sum(axis=0)

    def queue(self):
        pass  # this should be taken care of when the force/potential/etc is accessed


dproperties(
    ScaledForceComponent,
    ["weight", "scaling", "f", "pots", "pot", "virs", "vir", "extras"],
)


class MTSForces:
    """The forces and virials of the components that act at a level of a
    multiple time stepping scheme, weighted by their weights and by their MTS
    weights at that level.

    Depend objects:
       f: The forces on the beads.
       virs: The virials of the beads.
       vir: Their sum.
    """

    def __init__(self, parent, level):
        self.nbeads = parent.nbeads
        self.natoms = parent.natoms
        self.mforces = parent.mforces
        self.mrpc = parent.mrpc
        self.level = level

        self.f_ = depend_array(
            name="f",
            value=np.zeros((self.nbeads, 3 * self.natoms)),
            func=self.get_forces_mts,
        )

        self.virs_ = depend_array(
            name="virs", value=np.zeros((self.nbeads, 3, 3)), func=self.get_virs_mts
        )

        self.vir_ = depend_array(
            name="vir",
            value=np.zeros((3, 3)),
            func=self.get_vir_mts,
            dependencies=[self.virs_],
        )

        for ff in self.mforces:
            # add dependencies from the forces that actually depend on this
            if _in_mts_level(ff, level):
                ddeps(self.f_, [ff.f_, ff.weight_])
                ddeps(self.virs_, [ff.virs_, ff.weight_])

    def queue_mts(self):
        """Submits all the required force calculations to the forcefields."""

        level = self.level
        for ff in self.mforces:
            # do not queue forces which have zero weight
            if _in_mts_level(ff, level) and ff.weight != 0:
                ff.queue()

    def get_forces_mts(self):
        """Fetches ONLY the forces associated with a given MTS level."""

        self.queue_mts()
        return _weighted_sum(
            self.mforces, self.mrpc, "f", (self.nbeads, 3 * self.natoms), self.level
        )

    def get_vir_mts(self):
        """Fetches ONLY the total virial associated with a given MTS level."""
        return self.virs.sum(axis=0)

    def get_virs_mts(self):
        """Fetches ONLY the total virial associated with a given MTS level."""

        self.queue_mts()
        return _weighted_sum(
            self.mforces, self.mrpc, "virs", (self.nbeads, 3, 3), self.level
        )


dproperties(MTSForces, ["f", "vir", "virs"])


class Forces:
    """Class that gathers all the forces together.
    Collects many forcefield instances and parallelizes getting the forces
    in a PIMD environment.

    Attributes:
       natoms: An integer giving the number of atoms.
       nbeads: An integer giving the number of beads.
       nforces: An integer giving the number of ForceBeads objects.
       mforces: A list of all the forcefield objects.
       mbeads: A list of all the beads objects. Some of these may be contracted
          ring polymers, with a smaller number of beads than of the simulation.
       mrpc: A list of the objects containing the functions required to
          contract the ring polymers of the different forcefields.

    Depend objects:
       f: An array containing the components of the force. Depends on each
          replica's ufvx list.
       pots: A list containing the potential energy for each system replica.
          Depends on each replica's ufvx list.
       virs: A list containing the virial tensor for each system replica.
          Depends on each replica's ufvx list.
       extras: The extras of the components, as a dictionary. The keys in the
          `interpolate_extras` of a component are treated as the energy: contracted
          or expanded to the beads of the system, weighted and summed over the
          components, as numpy arrays of shape (nbeads, ...). "raw" holds the raw
          strings of the components, concatenated for each bead. The other keys
          map to the list of values of each bead (the last component that returns
          a key wins), and are dropped for contracted components.
       pot: The sum of the potential energy of the replicas.
       vir: The sum of the virial tensor of the replicas.
    """

    def __init__(self):
        self.bound = False
        # auxiliary systems for the finite differences of the fourth-order term
        self.aux = {}

    def add_component(self, nbeads, nrpc, nforces):
        self.mrpc.append(nrpc)
        self.mbeads.append(nbeads)
        self.mforces.append(nforces)
        self.nforces += 1
        ddeps(self.f_, [nforces.f_])
        ddeps(self.pots_, [nforces.pots_])
        ddeps(self.virs_, [nforces.virs_])
        ddeps(self.extras_, [nforces.extras_])
        for dep in (self.f_, self.pots_, self.virs_, self.extras_):
            ddeps(dep, [nforces.weight_])

    def bind(self, beads, cell, fcomponents, fflist, open_paths, output_maker):
        """Binds beads, cell and forces to the forcefield.


        Args:
           beads: Beads object from which the bead positions are taken.
           cell: Cell object from which the system box is taken.
           fcomponents: A list of different objects for each force type.
              For example, if ring polymer contraction is being used,
              then there may be separate forces for the long and short
              range part of the potential.
           fflist: A list of forcefield objects to use to calculate the potential,
              forces and virial for each force type. To clarify: fcomponents are the
              names and parameters of forcefields that are active for a certain
              system. fflist contains the overall list of force providers,
              and one typically has just one per force kind.
        """

        self.natoms = beads.natoms
        self.nbeads = beads.nbeads
        self.beads = beads
        self.cell = cell
        self.bound = True
        self.nforces = len(fcomponents)

        # prepares force components
        self.fcomp = fcomponents

        self.ff = fflist
        self.open_paths = open_paths
        self.output_maker = output_maker

        # fflist should be a dictionary of forcefield objects
        self.mforces = []
        self.mbeads = []
        self.mrpc = []

        # a "function factory" to generate functions to automatically update
        # contracted paths
        def make_rpc(rpc, beads):
            if rpc.noop:
                # the positions are copied when stored in the contracted beads
                return lambda: beads.q
            return lambda: rpc.b1tob2(beads.q)

        # creates new force objects, possibly acting on contracted path
        # representations. note that this new object is always created even if no contraction is required.
        for fc in self.fcomp:
            # creates an automatically-updated contracted beads object
            newb = fc.nbeads
            # if the number of beads for this force component is unspecified,
            # assume full force evaluation
            if newb == 0 or newb > beads.nbeads:
                newb = beads.nbeads
            newforce = ForceComponent(
                ffield=fc.ffield,
                name=fc.name,
                nbeads=newb,
                weight=fc.weight,
                mts_weights=fc.mts_weights,
                interpolate_extras=fc.interpolate_extras,
                epsilon=fc.epsilon,
            )
            newbeads = Beads(beads.natoms, newb)
            newrpc = nm_rescale(beads.nbeads, newb, open_paths=self.open_paths)

            # the beads positions for this force components are obtained
            # automatically, when needed, as a contraction of the full beads
            dsetfunc(newbeads.q_, make_rpc(newrpc, beads), [beads.q_])

            # now we create a new forcecomponent which is bound to newbeads!
            newforce.bind(newbeads, cell, fflist, output_maker=self.output_maker)

            # adds information we will later need to the appropriate lists.
            self.mbeads.append(newbeads)
            self.mforces.append(newforce)
            self.mrpc.append(newrpc)

        # now must expose an interface that gives overall forces
        self.f_ = depend_array(
            name="f",
            value=np.zeros((self.nbeads, 3 * self.natoms)),
            func=self.f_combine,
            dependencies=[ff.f_ for ff in self.mforces],
        )

        # collection of pots and virs from individual ff objects
        self.pots_ = depend_array(
            name="pots",
            value=np.zeros(self.nbeads, float),
            func=self.pot_combine,
            dependencies=[ff.pots_ for ff in self.mforces],
        )

        # must take care of the virials!
        self.virs_ = depend_array(
            name="virs",
            value=np.zeros((self.nbeads, 3, 3), float),
            func=self.vir_combine,
            dependencies=[ff.virs_ for ff in self.mforces],
        )

        self.extras_ = depend_value(
            name="extras",
            value={},
            func=self.extra_combine,
            dependencies=[ff.extras_ for ff in self.mforces],
        )

        # total potential and total virial
        self.pot_ = depend_value(
            name="pot", func=(lambda: self.pots.sum()), dependencies=[self.pots_]
        )

        self.vir_ = depend_array(
            name="vir",
            func=self.get_vir,
            value=np.zeros((3, 3)),
            dependencies=[self.virs_],
        )

        # SC forces and potential
        self.alpha_ = depend_value(name="alpha", value=0.0)

        # The number of MTS levels
        self.nmtslevels_ = depend_value(
            name="nmtslevels", value=0, func=self.get_nmtslevels
        )

        if len(self.mforces) > 0:
            self.mts_forces = [MTSForces(self, i) for i in range(self.nmtslevels)]

        # This will be piped from normalmodes
        self.omegan2_ = depend_value(name="omegan2", value=0)

        # The Suzuki-Chin difference potential
        self.potssc_ = depend_array(
            name="potssc",
            value=np.zeros(self.nbeads, float),
            dependencies=[
                self.beads.m_,
                self.f_,
                self.pots_,
                self.alpha_,
                self.omegan2_,
            ],
            func=self.get_potssc,
        )

        self.potsc_ = depend_value(
            name="potsc",
            dependencies=[self.potssc_],
            func=(lambda: self.potssc.sum()),
        )

        # The coefficients of the physical and the |f|^2 terms
        self.coeffsc_part_1_ = depend_array(
            name="coeffsc_part_1",
            value=np.zeros((self.nbeads, 1), float),
            func=self.get_coeffsc_part_1,
        )

        self.coeffsc_part_2_ = depend_array(
            name="coeffsc_part_2",
            value=np.zeros((self.nbeads, 1), float),
            dependencies=[self.alpha_, self.omegan2_],
            func=self.get_coeffsc_part_2,
        )

        # A list that contains the high order component of the force and the virial
        self.fvir_4th_order_ = depend_value(
            name="fvir_4th_order",
            value=[None, None],
            dependencies=[self.beads.m_, self.f_, self.pots_, self.alpha_],
            func=self.fvir_4th_order_combine,
        )

        # The high order component of the Suzuki-Chin force.
        self.f_4th_order_ = depend_array(
            name="f_4th_order",
            value=np.zeros((self.nbeads, 3 * self.natoms), float),
            dependencies=[self.fvir_4th_order_],
            func=(lambda: self.fvir_4th_order[0]),
        )

        self.fsc_part_1_ = depend_array(
            name="fsc_part_1",
            value=np.zeros((self.nbeads, 3 * self.natoms), float),
            dependencies=[self.coeffsc_part_1_, self.f_],
            func=self.get_fsc_part_1,
        )

        self.fsc_part_2_ = depend_array(
            name="fsc_part_2",
            value=np.zeros((self.nbeads, 3 * self.natoms), float),
            dependencies=[self.coeffsc_part_2_, self.f_4th_order_],
            func=self.get_fsc_part_2,
        )

        self.fsc_ = depend_array(
            name="fsc",
            value=np.zeros((self.nbeads, 3 * self.natoms), float),
            dependencies=[self.fsc_part_1_, self.fsc_part_2_],
            func=self.get_fsc,
        )

        # The high order component of the Suzuki-Chin virial.
        self.vir_4th_order_ = depend_array(
            name="vir_4th_order",
            value=np.zeros((self.nbeads, 3, 3), float),
            dependencies=[self.fvir_4th_order_],
            func=(lambda: self.fvir_4th_order[1]),
        )

        self.virssc_part_1_ = depend_array(
            name="virssc_part_1",
            value=np.zeros((self.nbeads, 3, 3), float),
            dependencies=[self.coeffsc_part_1_, self.virs_],
            func=self.get_virssc_part_1,
        )

        self.virssc_part_2_ = depend_array(
            name="virssc_part_2",
            value=np.zeros((self.nbeads, 3, 3), float),
            dependencies=[self.coeffsc_part_2_, self.vir_4th_order_],
            func=self.get_virssc_part_2,
        )

        self.virssc_ = depend_array(
            name="virssc",
            value=np.zeros((self.nbeads, 3, 3), float),
            dependencies=[self.virssc_part_1_, self.virssc_part_2_],
            func=self.get_virssc,
        )

        self.virsc_ = depend_value(
            name="virsc",
            dependencies=[self.virssc_],
            func=(lambda: self.virssc.sum(axis=0)),
        )

        # Add dependencies from the force weights, that are applied here when the total
        # force is assembled from its components

        weights = [fc.weight_ for fc in self.mforces]
        ddeps(self.f_, weights)
        ddeps(self.pots_, weights)
        ddeps(self.virs_, weights)
        ddeps(self.extras_, weights)

    def clone(self, beads, cell):
        """Duplicates the force object, so that it can be used to compute forces
        using the same forcefields, but with a separate beads and cell object
        (that can e.g. be created by cloning existing beads and cell objects).
        This is useful whenever one wants to compute the same forces and energies
        without touching the physical system, e.g. for use in internal loops of
        geometry optimizers, when attempting Monte Carlo moves or for property
        calculation.

        Args:
           beads: Beads object that the clone should be bound to
           cell: Cell object that the clone should be bound to

        Returns: The copy of the Forces object
        """

        if not self.bound:
            raise ValueError("Cannot copy a forces object that has not yet been bound.")
        if not type(beads) is Beads:
            raise ValueError(
                "The 'beads' argument must be provided and be a Beads object."
            )
        if not type(cell) is Cell:
            raise ValueError(
                "The 'cell' argument must be provided and be a Cell object."
            )

        nforce = Forces()
        nbeads = beads
        ncell = cell
        nforce.bind(
            nbeads, ncell, self.fcomp, self.ff, self.open_paths, self.output_maker
        )
        return nforce

    def dump_state(self):
        """
        Stores a (deep) copy of the internal state of a force object.
        See `transfer_forces` for an explanation of the usage and
        logic.
        """

        force_data = []
        for k in range(len(self.mforces)):
            mself = self.mforces[k]
            # access the actual data values

            # forces that are not up to date (e.g. those of components with
            # zero weight, that are not computed) are not stored
            results = mself.results_
            force_data.append(
                (
                    # will also need bead positions, see _load_state
                    deepcopy(mself.beads.q),
                    None if results.tainted() else deepcopy(results._load()),
                )
            )
        return force_data

    def load_state(self, state):
        """
        Loads a previously-saved state of a force object. See
        `transfer_forces` for an explanation of the usage and
        logic.
        """

        if len(state) != len(self.mforces):
            raise ValueError(
                "Cannot load force state from an force with a different numbers of components"
            )

        for k in range(len(self.mforces)):
            mself = self.mforces[k]
            q, results = state[k]

            # this is VERY subtle. beads in this force component are
            # obtained as a contraction, and so are computed automatically.
            # when we set the main q, these get marked as tainted.
            # then we copy the force value, and set the force as untainted.
            # next time we touch the main q, the tainting does not get
            # propagated, because the contracted q is already marked as tainted,
            # so the force does not get updated. we can fix this by copying
            # the value of the contracted bead, so that it's marked as NOT
            # tainted - it should not be as it's an internal of the force and
            # therefore get copied
            mself.beads.q_.set(q, manual=False)
            if results is None:
                mself.results_.taint()
            elif len(results[3]) != mself.nbeads:
                raise ValueError(
                    "Cannot copy forces between objects with different numbers of beads for the "
                    + str(k)
                    + "th component"
                )
            else:
                mself.results_.set(results, manual=False)

    def transfer_forces(self, refforce):
        """Low-level function copying over the value of a second force object,
        triggering updates but un-tainting this force depends themselves.

        We have noted that in some corner cases it is necessary to copy only
        the values of updated forces rather than the full depend object, in order to
        avoid triggering a repeated call to the client code that is potentially
        very costly. This happens routinely in geometry relaxation routines,
        for example.

        The transfer can also be implemented in a delayed way, by using
        separately set_state and dump_state. For instance, one can

        old_state = force.dump_state()

        do-something-that-changes-force

        force.load_state(old_state)
        """

        self.load_state(refforce.dump_state())

    def transfer_forces_manual(
        self, new_q, new_v, new_forces, new_x=None, vir=np.zeros((3, 3))
    ):
        """Manual (and flexible) version of the transfer forces function.
        Instead of passing a force object, list with vectors are passed
        Expected shape and sizes:
          - new_q list of length equal to number of force type, containing the beads positions
          - new_v list of length equal to number of force type, containing the beads potential energy
          - new_f list of length equal to number of force type, containing the beads forces
          - new_x list of length equal to number of force type, containing the beads extras
        """
        msg = "Inconsistent dimensions inside transfer_forces_manual"
        assert len(self.mforces) == len(new_q), msg
        assert len(self.mforces) == len(new_v), msg
        assert len(self.mforces) == len(new_forces), msg
        if new_x == None:
            new_x = [{"raw": [None] * self.nbeads}] * len(self.mforces)
            info("WARNING: No extras information has been passed.", verbosity.debug)

        assert len(self.mforces) == len(new_x), msg

        for k in range(len(self.mforces)):
            mv = np.array(to_numpy(new_v[k]))
            mf = np.array(to_numpy(new_forces[k]))
            mq = new_q[k]
            mextra = new_x[k]
            mself = self.mforces[k]
            assert mq.shape == mf.shape, msg
            assert mq.shape[0] == mv.shape[0], msg
            assert mself.nbeads == mv.shape[0], msg
            assert mself.nbeads == mq.shape[0], msg

            mself.beads.q_.set(mq, manual=False)
            extras = [{key: mextra[key][b] for key in mextra} for b in range(mself.nbeads)]
            virs = np.tile(vir, (mself.nbeads, 1, 1))
            mself.results_.set(
                [xp.asarray(mv), xp.asarray(mf), xp.asarray(virs), extras], manual=False
            )

    def queue(self):
        """Submits all the required force calculations to the forcefields."""

        is_tainted = False
        for ff in self.mforces:
            if ff.weight != 0:  # do not compute forces which have zero weight
                is_tainted = ff.queue() or is_tainted
        return is_tainted

    def get_vir(self):
        """Sums the virial of each forcefield.

        Not the actual system virial.

        Returns:
            Virial sum.
        """

        return self.virs.sum(axis=0)

    def pots_component(self, index, weighted=True, interpolate=True):
        """Fetches the index^th component of the total potential."""

        return self._component(index, "pots", (), weighted, interpolate)

    def forces_component(self, index, weighted=True, interpolate=True):
        """Fetches the index^th component of the total force."""

        return self._component(index, "f", (3 * self.natoms,), weighted, interpolate)

    def virs_component(self, index, weighted=True, interpolate=True):
        """Fetches the index^th component of the virial."""

        return self._component(index, "virs", (3, 3), weighted, interpolate)

    def _component(self, index, name, shape, weighted, interpolate):
        """Returns a copy of the array `name` of the index^th force component,
        whose values for each bead have shape `shape`, multiplied by the
        weight of the component if `weighted`, and expanded to the beads of
        the system if `interpolate`."""

        component, rpc = self.mforces[index], self.mrpc[index]
        if weighted and component.weight == 0:
            nbeads = rpc.nbeads1 if interpolate else rpc.nbeads2
            return xp.zeros((nbeads, *shape))
        value = xp.asarray(getattr(component, name), copy=True)
        if weighted:
            value *= component.weight
        if interpolate:
            value = rpc.b2tob1(value)
        return value

    def _expanded_fvirs(self, index):
        """Returns the forces and the virials of the index^th force component,
        expanded to the beads of the system."""

        rpc, component = self.mrpc[index], self.mforces[index]
        return rpc.b2tob1(component.f), rpc.b2tob1(component.virs)

    def extras_component(self, index):
        """Fetches extras that are computed for one specific force component.

        Does not attempt to apply weights or interpolate, always returns raw stuff.
        """

        return self.mforces[index].extras

    def forcesvirs_4th_order(self, index):
        """Fetches the 4th order |f^2| correction to the force vector and the virial associated with a given component."""

        if self.nbeads % 2 != 0:
            raise ValueError(
                "Suzuki-Chin factorization requires an even number of beads"
            )

        # calculates the finite displacement, along f / m.
        fm = self.f / self.beads.m3
        eps = self.mforces[index].epsilon
        foverm = xp.sqrt((fm * fm).sum() / (self.nbeads * self.natoms))
        if float(foverm) > 1e-20:
            delta = abs(eps) / foverm
        else:  # defaults to eps if otherwise we'd get an 1/0
            delta = abs(eps)
        dq = delta * fm

        # stores the force component.
        fbase, vbase = self._expanded_fvirs(index)

        # uses a fwd difference if epsilon > 0.
        if self.mforces[index].epsilon > 0.0:
            # gives an error if RPC is used with a fwd difference.
            # The problem is that the finite difference is computed as [f(q + e.f) - f(q)].e^-1,
            # where f(q + e.f) is computed on a smaller rimg polymer of P / 2 beads which is made
            # by displacing the odd beads of the original ring polymer. Upon contraction,  it yields
            # a ring polymer that is different from the contracted one on which f(q) was computed.
            # This makes the FD incorrect. A fix could be to ALWAYS compute the odd and the even
            # beads on two different ring polymers of P / 2 beads, but centered difference + RPC
            # seems to be a neater solution. Anyway, the cost of a CNT-DIFF in comparison to a FWD-DIFF
            # is marginal when RPC + MTS is used.

            if self.mforces[index].nbeads != self.nbeads:
                raise ValueError(
                    "High-order PIMD with ring-polymer contraction needs a centered "
                    "finite difference (a negative epsilon)"
                )

            # for the case of alpha = 0, only odd beads are displaced.
            if self.alpha == 0:
                # we use an aux force evaluator with half the number of beads.
                dbeads, dforces = self.aux_system(self.nbeads // 2)

                # displaces odd beads only.
                dbeads.q = self.beads.q[1::2] - dq[1::2]
                fminus, vminus = dforces._expanded_fvirs(index)

                # calculates the finite difference on the odd beads (the
                # copies on the even beads are zeroed below)
                f_4th_order = 2.0 * (xp.repeat(fminus, 2, axis=0) - fbase) / delta
                v_4th_order = 2.0 * (xp.repeat(vminus, 2, axis=0) - vbase) / delta

            # For the case of alpha != 0, all the beads are displaced.
            else:
                # we use an aux force evaluator with the same number of beads.
                dbeads, dforces = self.aux_system(self.nbeads)

                # displaces the beads.
                dbeads.q = self.beads.q + dq
                fplus, vplus = dforces._expanded_fvirs(index)

                # calculates the finite difference.
                f_4th_order = 2.0 * (fbase - fplus) / delta
                v_4th_order = 2.0 * (vbase - vplus) / delta

        # uses a centered difference for epsilon  < 0.
        if self.mforces[index].epsilon < 0.0:
            # we use an aux force evaluator with the same number of beads.
            dbeads, dforces = self.aux_system(self.nbeads)

            # for the case of alpha = 0, only odd beads are displaced, on a ring whose
            # first and second halves hold the + and - displacements; contracted
            # components would mix them, and displace the full ring instead (the
            # terms of its even beads are zeroed below)
            if self.alpha == 0 and self.mrpc[index].noop:
                # the first half of the aux beads are fwd displaced while the second half are bkwd displaced configurations.
                dbeads.q_[: self.nbeads // 2] = self.beads.q[1::2] + dq[1::2]
                dbeads.q_[-self.nbeads // 2 :] = self.beads.q[1::2] - dq[1::2]
                fplusminus, vplusminus = dforces._expanded_fvirs(index)

                # calculates the centred difference, 2 (f- - f+) / (2 delta), on
                # the odd beads (the copies on the even beads are zeroed below)
                h = self.nbeads // 2
                df = xp.repeat(fplusminus[h:] - fplusminus[:h], 2, axis=0)
                dv = xp.repeat(vplusminus[h:] - vplusminus[:h], 2, axis=0)
                f_4th_order = df / delta
                v_4th_order = dv / delta

            # otherwise, all the beads are displaced.
            else:
                # displaces the beads.
                dbeads.q = self.beads.q + dq
                fplus, vplus = dforces._expanded_fvirs(index)
                dbeads.q = self.beads.q - dq
                fminus, vminus = dforces._expanded_fvirs(index)

                # calculates the centred difference, 2 (f- - f+) / (2 delta).
                f_4th_order = (fminus - fplus) / delta
                v_4th_order = (vminus - vplus) / delta

        # the strain derivative of |f|^2 / m also contains 2 (f / m) x f, summed over
        # the components (Kapil et al., JCTC 15, 3237 (2019), Eq. 9), in the upper
        # triangle like the virials
        u = fm.reshape((self.nbeads, -1, 3))
        fk = fbase.reshape((self.nbeads, -1, 3))
        v_4th_order = v_4th_order + 2.0 * xp.triu(xp.swapaxes(u, 1, 2) @ fk)
        if self.alpha == 0:
            # only the odd beads have non-zero coefficients
            odd = (xp.arange(self.nbeads) % 2)[:, None]
            f_4th_order = f_4th_order * odd
            v_4th_order = v_4th_order * odd[:, :, None]

        # returns the 4th order |f^2| correction.
        return [f_4th_order, v_4th_order]

    def aux_system(self, nbeads):
        """Returns the beads and the forces of an auxiliary system with nbeads
        beads, on the current cell, for the fourth-order finite differences."""

        if nbeads not in self.aux:
            beads, cell = self.beads.clone(nbeads), self.cell.clone()
            self.aux[nbeads] = (beads, cell, self.clone(beads, cell))
        beads, cell, forces = self.aux[nbeads]
        cell.h = self.cell.h
        return beads, forces

    def get_nmtslevels(self):
        """Returns the total number of mts levels."""

        nm = len(self.mforces[0].mts_weights)
        if all(len(x.mts_weights) == nm for x in self.mforces):
            return nm
        else:
            raise ValueError(
                "The mts_weights of all the force components are not the same."
            )

    def f_combine(self):
        """Obtains the total force vector."""

        self.queue()
        return _weighted_sum(
            self.mforces, self.mrpc, "f", (self.nbeads, 3 * self.natoms)
        )

    def fvir_4th_order_combine(self):
        """Obtains the total fourth order |f^2| correction to the force vector and the virial."""

        rf = xp.zeros((self.nbeads, 3 * self.natoms))
        rv = xp.zeros((self.nbeads, 3, 3))

        for k, component in enumerate(self.mforces):
            scale = _scale(component)
            # the finite differences are skipped for the components that do
            # not contribute
            if scale is not None and scale != 0:
                f, v = self.forcesvirs_4th_order(k)
                rf += scale * f
                rv += scale * v
        return [rf, rv]

    def pot_combine(self):
        """Obtains the potential energy for each forcefield."""

        self.queue()
        return _weighted_sum(self.mforces, self.mrpc, "pots", (self.nbeads,))

    def extra_combine(self):
        """Combines the extra dictionaries for each forcefield for each bead."""

        self.queue()

        re = {}
        for component, rpc in zip(self.mforces, self.mrpc):
            scale = _scale(component)
            if scale is None:
                continue
            # combines the extras from the different force components
            for e, v in component.extras.items():
                if e in component.interpolate_extras:
                    # extras that are tagged as interpolate_extras are treated exactly as if they were an energy/force/stress
                    # (extras are stored as numpy arrays)
                    v = scale * to_numpy(v if rpc.noop else rpc.b2tob1(v))
                    if e in re:
                        # if multiple forcefields have the same "promoted extras" these get summed
                        re[e] += v
                    else:
                        # the first is just added to the dict, so we don't have to worry about extras having different shapes
                        re[e] = v
                else:
                    # other extras are not touched, just accummulated. if there is a contraction, they are dropped because there is no meaningful way to do a contraction
                    if self.nbeads != component.nbeads:
                        warning(
                            "Extra field '"
                            + e
                            + "' cannot be contracted unless interpreted as physical properties. Will just drop it",
                            verbosity.high,
                        )
                    else:
                        if e == "raw":
                            # concatenates raw outputs
                            if e in re:
                                # concatenation must happen bead per bead
                                for ib in range(self.nbeads):
                                    re["raw"][ib] += v[ib]
                            else:
                                # a copy, as the strings of the other components are appended
                                re["raw"] = list(v)
                        else:
                            if e in re:
                                warning(
                                    "Extra field '"
                                    + e
                                    + "' appears in multiple forcefields will be overwritten unless interpreted as physical property",
                                    verbosity.high,
                                )
                            re[e] = v
        return re

    def vir_combine(self):
        """Obtains the virial tensor for each forcefield."""

        self.queue()
        return _weighted_sum(self.mforces, self.mrpc, "virs", (self.nbeads, 3, 3))

    def get_potssc(self):
        """Obtains Suzuki-Chin contribution to the potential."""
        if self.nbeads % 2 != 0:
            raise ValueError(
                "Suzuki-Chin factorization requires an even number of beads"
            )

        # this evaluates the square forces contribution to the SC potential (only the difference with the Trotter potential is returned)
        return self.coeffsc_part_1[:, 0] * self.pots + self.coeffsc_part_2[:, 0] * (
            self.f / self.beads.m3 * self.f
        ).sum(axis=1)

    def get_fsc_part_1(self):
        """Obtains the linear component of Suzuki-Chin correction to the force."""

        return self.coeffsc_part_1 * self.f

    def get_fsc_part_2(self):
        """Obtains the quadratic component of Suzuki-Chin correction to the force."""

        return self.coeffsc_part_2 * self.f_4th_order

    def get_virssc_part_1(self):
        """Obtains the linear component of Suzuki-Chin correction to the virial."""
        return self.coeffsc_part_1.reshape((self.beads.nbeads, 1, 1)) * self.virs

    def get_virssc_part_2(self):
        """Obtains the quadratic component of Suzuki-Chin correction to the virial."""
        return (
            self.coeffsc_part_2.reshape((self.beads.nbeads, 1, 1)) * self.vir_4th_order
        )

    def get_fsc(self):
        """Obtains the total Suzuki-Chin correction to the force."""
        return self.fsc_part_1 + self.fsc_part_2

    def get_virssc(self):
        """Obtains the total Suzuki-Chin correction to the force."""
        return self.virssc_part_1 + self.virssc_part_2

    def get_coeffsc_part_1(self):
        """Obtains the coefficients of the linear part of the Suzuki-Chin correction."""

        rc = np.zeros((self.beads.nbeads, 1))
        rc[0::2] = -1.0 / 3.0
        rc[1::2] = 1.0 / 3.0
        return xp.asarray(rc)

    def get_coeffsc_part_2(self):
        """Obtains the coefficients of the linear part of the Suzuki-Chin correction."""

        rc = np.zeros((self.beads.nbeads, 1))
        rc[0::2] = self.alpha / self.omegan2 / 9.0
        rc[1::2] = (1.0 - self.alpha) / self.omegan2 / 9.0
        return xp.asarray(rc)


dproperties(
    Forces,
    [
        "f",
        "pots",
        "virs",
        "extras",
        "pot",
        "vir",
        "alpha",
        "nmtslevels",
        "omegan2",
        "potssc",
        "potsc",
        "coeffsc_part_1",
        "coeffsc_part_2",
        "fvir_4th_order",
        "f_4th_order",
        "fsc_part_1",
        "fsc_part_2",
        "fsc",
        "vir_4th_order",
        "virssc_part_1",
        "virssc_part_2",
        "virssc",
        "virsc",
    ],
)

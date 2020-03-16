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
import sys
import threading
from copy import deepcopy

import numpy as np

from ipi.utils.softexit import softexit
from ipi.utils.messages import verbosity, warning, info
from ipi.utils.depend import *
from ipi.utils.nmtransform import nm_rescale
from ipi.engine.beads import Beads


__all__ = ['Forces', 'ForceComponent']


fbuid = 0


class ForceBead(dobject):

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
       ufvx: A list of the form [pot, f, vir]. These quantities are calculated
          all at one time by the driver, so are collected together. Each separate
          object is then taken from the list. Depends on the atom positions and
          the system box.
       extra: A string containing some formatted output returned by the client. Depends on ufvx.
       pot: A float giving the potential energy of the system. Depends on ufvx.
       f: An array containing all the components of the force. Depends on ufvx.
       fx: A slice of f containing only the x components of the forces.
       fy: A slice of f containing only the y components of the forces.
       fz: A slice of f containing only the z components of the forces.
       vir: An array containing the components of the virial tensor in upper
          triangular form, not divided by the volume. Depends on ufvx.
       request: a handle to the request that has been filed by the FF object
    """

    def __init__(self):
        """Initialises ForceBead."""
        dself = dd(self)

        # ufvx is a list [ u, f, vir, extra ]  which stores the results of the force calculation
        dself.ufvx = depend_value(name="ufvx", func=self.get_all)
        self._threadlock = threading.Lock()
        self.request = None
        self._getallcount = 0

    def bind(self, atoms, cell, ff):
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
        dself = dd(self)

        # ufv depends on the atomic positions and on the cell
        dself.ufvx.add_dependency(dd(self.atoms).q)
        dself.ufvx.add_dependency(dd(self.cell).h)

        # potential and virial are to be extracted very simply from ufv
        dself.pot = depend_value(name="pot", func=self.get_pot,
                                 dependencies=[dself.ufvx])

        dself.vir = depend_array(name="vir", value=np.zeros((3, 3), float),
                                 func=self.get_vir,
                                 dependencies=[dself.ufvx])

        # NB: the force requires a bit more work, to define shortcuts to xyz
        # slices without calculating the force at this point.
        fbase = np.zeros(atoms.natoms * 3, float)
        dself.f = depend_array(name="f", value=fbase, func=self.get_f,
                               dependencies=[dself.ufvx])

        dself.extra = depend_value(name="extra", func=self.get_extra,
                                   dependencies=[dself.ufvx])

        dself.fx = depend_array(name="fx", value=fbase[0:3 * atoms.natoms:3])
        dself.fy = depend_array(name="fy", value=fbase[1:3 * atoms.natoms:3])
        dself.fz = depend_array(name="fz", value=fbase[2:3 * atoms.natoms:3])
        dcopy(dself.f, dself.fx)
        dcopy(dself.f, dself.fy)
        dcopy(dself.f, dself.fz)

    def queue(self):
        """Sends the job to the interface queue directly.

        Allows the ForceBead object to ask for the ufvx list of each replica
        directly without going through the get_all function. This allows
        all the jobs to be sent at once, allowing them to be parallelized.
        """

        with self._threadlock:
            if self.request is None and dd(self).ufvx.tainted():
                self.request = self.ff.queue(self.atoms, self.cell, reqid=self.uid)

    def get_all(self):
        """Driver routine.

        When one of the force, potential or virial are called, this sends the
        atoms and cell to the client code, requesting that it calculates the
        potential, forces and virial tensor. This then waits until the
        driver is finished, and then returns the ufvx list.

        Returns:
           A list of the form [potential, force, virial, extra].
        """

        # because we thread over many systems and outputs, we might get called
        # more than once. keep track of how many times we are called so we
        # can make sure to wait until the last call has returned before we release
        with self._threadlock:
            self._getallcount += 1

        # this is converting the distribution library requests into [ u, f, v ]  lists
        t_start = time.time()
        if self.request is None:
            self.request = self.queue()

        # sleeps until the request has been evaluated
        while self.request["status"] != "Done":
            if self.request["status"] == "Exit" or softexit.triggered:
                # now, this is tricky. we are stuck here and we cannot return meaningful results.
                # if we return, we may as well output wrong numbers, or mess up things.
                # so we can only call soft-exit and wait until that is done. then kill the thread
                # we are in.
                softexit.trigger(" @ FORCES : cannot return so will die off here")
                while softexit.exiting:
                    time.sleep(self.ff.latency)
                sys.exit()
            time.sleep(self.ff.latency)
        # print diagnostics about the elapsed time
        info("# forcefield %s evaluated in %f (queue) and %f (dispatched) sec." % (self.ff.name, self.request["t_finished"] - self.request["t_queued"], self.request["t_finished"] - self.request["t_dispatched"]), verbosity.debug)

        # data has been collected, so the request can be released and a slot
        # freed up for new calculations
        result = self.request["result"]

        # reduce the reservation count (and wait for all calls to return)
        with self._threadlock:
            self._getallcount -= 1

        # releases just once, but wait for all requests to be complete
        if self._getallcount == 0:
            self.ff.release(self.request)
            self.request = None
        else:
            while self._getallcount > 0:
                time.sleep(self.ff.latency)

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

        return dstrip(self.ufvx[1])

    def get_vir(self):
        """Calls get_all routine of forcefield to update the virial.

        Returns:
           An array containing the virial in upper triangular form, not divided
           by the volume.
        """

        vir = dstrip(self.ufvx[2])
        vir[1, 0] = 0.0
        vir[2, 0:2] = 0.0
        return vir

    def get_extra(self):
        """Calls get_all routine of forcefield to update the extras string.

        Returns:
           A string containing all formatted additional output that the
           client might have produced.
        """

        return self.ufvx[3]


class ForceComponent(dobject):
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
       f: An array containing the components of the force. Depends on each
          replica's ufvx list.
       pots: A list containing the potential energy for each system replica.
          Depends on each replica's ufvx list.
       virs: A list containing the virial tensor for each system replica.
          Depends on each replica's ufvx list.
       pot: The sum of the potential energy of the replicas.
       vir: The sum of the virial tensor of the replicas.
       extras: Strings containing some formatted output returned by the client.
          Depends on each replica's ufvx list.
    """

    def __init__(self, ffield, nbeads=0, weight=1.0, name="", mts_weights=None, epsilon=-0.001):
        """Initializes ForceComponent

        Args:
           ffield: A model to be used to create the forcefield objects for all
              the replicas of the system.
           nbeads: The number of replicas.
           weight: A relative weight to be given to the values obtained with this
              forcefield. When the contribution of all the forcefields is
              combined to give a total force, the contribution of this forcefield
              will be weighted by this factor.
           name: The name of the forcefield.
           mts_weights: Weight of forcefield at each mts level.
        """

        self.ffield = ffield
        self.name = name
        self.nbeads = nbeads
        self.weight = depend_value(name="weight", value=weight)
        if mts_weights is None:
            self.mts_weights = np.asarray([])
        else:
            self.mts_weights = np.asarray(mts_weights)
        self.epsilon = epsilon

    def bind(self, beads, cell, fflist):
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
        dself = dd(self)

        # stores a copy of the number of atoms and of beads
        self.natoms = beads.natoms
        if (self.nbeads != beads.nbeads):
            raise ValueError("Binding together a Beads and a ForceBeads objects with different numbers of beads")

        # creates an array of force objects, which are bound to the beads
        # and the cell
        if not self.ffield in fflist:
            raise ValueError("Force component name '" + self.ffield + "' is not in the forcefields list")

        self.ff = fflist[self.ffield]

        self._forces = [];
        self.beads = beads
        for b in range(self.nbeads):
            new_force = ForceBead()
            new_force.bind(beads[b], cell, self.ff)
            self._forces.append(new_force)

        # f is a big array which assembles the forces on individual beads
        dself.f = depend_array(name="f",
                               value=np.zeros((self.nbeads, 3 * self.natoms)),
                               func=self.f_gather,
                               dependencies=[dd(self._forces[b]).f for b in
                                             range(self.nbeads)])

        # collection of pots and virs from individual beads
        dself.pots = depend_array(name="pots", value=np.zeros(self.nbeads, float),
                                  func=self.pot_gather,
                                  dependencies=[dd(self._forces[b]).pot for b in
                                                range(self.nbeads)])
        dself.virs = depend_array(name="virs", value=np.zeros((self.nbeads, 3, 3), float),
                                  func=self.vir_gather,
                                  dependencies=[dd(self._forces[b]).vir for b in
                                                range(self.nbeads)])
        dself.extras = depend_value(name="extras", value=np.zeros(self.nbeads, float),
                                    func=self.extra_gather,
                                    dependencies=[dd(self._forces[b]).extra for b in
                                                  range(self.nbeads)])

        # total potential and total virial
        dself.pot = depend_value(name="pot", func=(lambda: self.pots.sum()),
                                 dependencies=[dself.pots])
        dself.vir = depend_array(name="vir", func=self.get_vir, value=np.zeros((3, 3)),
                                 dependencies=[dself.virs])

    def queue(self):
        """Submits all the required force calculations to the interface."""

        # this should be called in functions which access u,v,f for ALL the beads,
        # before accessing them. it is basically pre-queueing so that the
        # distributed-computing magic can work
        for b in range(self.nbeads):
            self._forces[b].queue()

    def pot_gather(self):
        """Obtains the potential energy for each replica.

        Returns:
           A list of the potential energy of each replica of the system.
        """

        self.queue()
        return np.array([b.pot for b in self._forces], float)

    def extra_gather(self):
        """Obtains the potential energy for each replica.

        Returns:
           A list of the potential energy of each replica of the system.
        """

        self.queue()
        return [b.extra for b in self._forces]

    def vir_gather(self):
        """Obtains the virial for each replica.

        Returns:
           A list of the virial of each replica of the system.
        """

        self.queue()
        return np.array([b.vir for b in self._forces], float)

    def f_gather(self):
        """Obtains the force vector for each replica.

        Returns:
           An array with all the components of the force. Row i gives the force
           array for replica i of the system.
        """

        newf = np.zeros((self.nbeads, 3 * self.natoms), float)
        self.queue()
        for b in range(self.nbeads):
            newf[b] = dstrip(self._forces[b].f)

        return newf

    def get_vir(self):
        """Sums the virial of each replica.

        Not the actual system virial, as it has not been divided by either the
        number of beads or the cell volume.

        Returns:
            Virial sum.
        """

        vir = np.zeros((3, 3))
        for v in dstrip(self.virs):
            vir += v
        return vir


class ScaledForceComponent(dobject):
    def __init__(self, baseforce, scaling=1):

        self.bf = baseforce
        self.name = baseforce.name
        self.ffield = baseforce.ffield
        dself = dd(self)
        dself.scaling = depend_value(name="scaling", value=scaling)
        dself.f = depend_array(name="f",
                               func=lambda: self.scaling * self.bf.f if scaling != 0 else np.zeros((self.bf.nbeads, 3 * self.bf.natoms)),
                               value=np.zeros((self.bf.nbeads, 3 * self.bf.natoms)),
                               dependencies=[dd(self.bf).f, dself.scaling])
        dself.pots = depend_array(name="pots",
                                  func=lambda: self.scaling * self.bf.pots if scaling != 0 else np.zeros(self.bf.nbeads),
                                  value=np.zeros(self.bf.nbeads),
                                  dependencies=[dd(self.bf).pots, dself.scaling])
        dself.virs = depend_array(name="virs", func=lambda: self.scaling * self.bf.virs,
                                  value=np.zeros((self.bf.nbeads, 3, 3)),
                                  dependencies=[dd(self.bf).virs, dself.scaling])
        dself.extras = depend_array(name="extras", func=lambda: self.bf.extras,
                                    value=np.zeros(self.bf.nbeads),
                                    dependencies=[dd(self.bf).extras])

        # total potential and total virial
        dself.pot = depend_value(name="pot", func=(lambda: self.pots.sum()),
                                 dependencies=[dself.pots])
        dself.vir = depend_array(name="vir", func=self.get_vir, value=np.zeros((3, 3)),
                                 dependencies=[dself.virs])

        # pipes weight from the force, since the scaling is applied on top of that
        self.weight = depend_value(name="weight", value=0)
        dpipe(dd(self.bf).weight, dself.weight)
        self.mts_weights = self.bf.mts_weights

    def get_vir(self):
        """Sums the virial of each replica.

        Not the actual system virial, as it has not been divided by either the
        number of beads or the cell volume.

        Returns:
            Virial sum.
        """

        vir = np.zeros((3, 3))
        for v in dstrip(self.virs):
            vir += v
        return vir

    def queue(self):
        pass  # this should be taken care of when the force/potential/etc is accessed


class Forces(dobject):

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
       extras: A list containing the "extra" strings for each replica.
       pot: The sum of the potential energy of the replicas.
       vir: The sum of the virial tensor of the replicas.
    """

    def __init__(self):
        self.bound = False
        self.dforces = None
        self.dbeads = None
        self.dcell = None

    def add_component(self, nbeads, nrpc, nforces):
        self.mrpc.append(nrpc)
        self.mbeads.append(nbeads)
        self.mforces.append(nforces)
        self.nforces += 1
        dself = dd(self)
        dforces = dd(nforces)
        dself.f.add_dependency(dforces.f)
        dself.pots.add_dependency(dforces.pots)
        dself.virs.add_dependency(dforces.virs)
        dself.extras.add_dependency(dforces.extras)

    def bind(self, beads, cell, fcomponents, fflist, open_paths):
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
        self.fcomp = fcomponents
        self.ff = fflist
        self.open_paths = open_paths

        # fflist should be a dictionary of forcefield objects
        self.mforces = []
        self.mbeads = []
        self.mrpc = []

        dself = dd(self)

        # a "function factory" to generate functions to automatically update
        # contracted paths
        def make_rpc(rpc, beads):
            return lambda: rpc.b1tob2(dstrip(beads.q))

        # creates new force objects, possibly acting on contracted path
        # representations. note that this new object is always created even if no contraction is required.
        for fc in self.fcomp:

            # creates an automatically-updated contracted beads object
            newb = fc.nbeads
            # if the number of beads for this force component is unspecified,
            # assume full force evaluation
            if newb == 0 or newb > beads.nbeads: newb = beads.nbeads
            newforce = ForceComponent(ffield=fc.ffield, name=fc.name, nbeads=newb, weight=fc.weight, mts_weights=fc.mts_weights, epsilon=fc.epsilon)
            newbeads = Beads(beads.natoms, newb)
            newrpc = nm_rescale(beads.nbeads, newb, open_paths=self.open_paths)

            # the beads positions for this force components are obtained
            # automatically, when needed, as a contraction of the full beads
            dd(newbeads).q._func = make_rpc(newrpc, beads)
            for b in newbeads:
                # must update also indirect access to the beads coordinates
                dd(b).q._func = dd(newbeads).q._func

            # makes newbeads.q depend from beads.q
            dd(beads).q.add_dependant(dd(newbeads).q)

            # now we create a new forcecomponent which is bound to newbeads!
            newforce.bind(newbeads, cell, fflist)

            # adds information we will later need to the appropriate lists.
            self.mbeads.append(newbeads)
            self.mforces.append(newforce)
            self.mrpc.append(newrpc)

        # now must expose an interface that gives overall forces
        dself.f = depend_array(name="f", value=np.zeros((self.nbeads, 3 * self.natoms)),
                               func=self.f_combine,
                               dependencies=[dd(ff).f for ff in self.mforces])

        # collection of pots and virs from individual ff objects
        dself.pots = depend_array(name="pots", value=np.zeros(self.nbeads, float),
                                  func=self.pot_combine,
                                  dependencies=[dd(ff).pots for ff in self.mforces])

        # must take care of the virials!
        dself.virs = depend_array(name="virs", value=np.zeros((self.nbeads, 3, 3), float),
                                  func=self.vir_combine,
                                  dependencies=[dd(ff).virs for ff in self.mforces])

        dself.extras = depend_value(name="extras", value=np.zeros(self.nbeads, float),
                                    func=self.extra_combine,
                                    dependencies=[dd(ff).extras for ff in self.mforces])

        # total potential and total virial
        dself.pot = depend_value(name="pot", func=(lambda: self.pots.sum()),
                                 dependencies=[dself.pots])

        dself.vir = depend_array(name="vir", func=self.get_vir, value=np.zeros((3, 3)),
                                 dependencies=[dself.virs])

        # SC forces and potential
        dself.alpha = depend_value(name="alpha", value=0.0)

        # The number of MTS levels
        dself.nmtslevels = depend_value(name="nmtslevels", value=0, func=self.get_nmtslevels)

        # This will be piped from normalmodes
        dself.omegan2 = depend_value(name="omegan2", value=0)

        # The Suzuki-Chin difference potential
        dself.potssc = depend_array(name="potssc", value=np.zeros(self.nbeads, float),
                                    dependencies=[dd(self.beads).m, dself.f, dself.pots, dself.alpha,
                                                  dself.omegan2],
                                    func=self.get_potssc)

        dself.potsc = depend_value(name="potsc",
                                   dependencies=[dself.potssc],
                                   func=(lambda: self.potssc.sum()))

        # The coefficients of the physical and the |f|^2 terms
        dself.coeffsc_part_1 = depend_array(name="coeffsc_part_1", value=np.zeros((self.nbeads, 1), float),
                                            func=self.get_coeffsc_part_1)

        dself.coeffsc_part_2 = depend_array(name="coeffsc_part_2", value=np.zeros((self.nbeads, 1), float),
                                            dependencies=[dself.alpha, dself.omegan2], func=self.get_coeffsc_part_2)

        # A list that contains the high order component of the force and the virial
        dself.fvir_4th_order = depend_value(name="fvir_4th_order", value=[None, None],
                                            dependencies=[dd(self.beads).m, dself.f, dself.pots],
                                            func=self.fvir_4th_order_combine)

        # The high order component of the Suzuki-Chin force.
        dself.f_4th_order = depend_array(name="f_4th_order", value=np.zeros((self.nbeads, 3 * self.natoms), float),
                                         dependencies=[dself.fvir_4th_order],
                                         func=(lambda: self.fvir_4th_order[0]))

        dself.fsc_part_1 = depend_array(name="fsc_part_1", value=np.zeros((self.nbeads, 3 * self.natoms), float),
                                        dependencies=[dself.coeffsc_part_1, dself.f],
                                        func=self.get_fsc_part_1)

        dself.fsc_part_2 = depend_array(name="fsc_part_2", value=np.zeros((self.nbeads, 3 * self.natoms), float),
                                        dependencies=[dself.coeffsc_part_2, dself.f_4th_order],
                                        func=self.get_fsc_part_2)

        dself. fsc = depend_array(name="fsc", value=np.zeros((self.nbeads, 3 * self.natoms), float),
                                  dependencies=[dself.fsc_part_1, dself.fsc_part_2],
                                  func=self.get_fsc)

        # The high order component of the Suzuki-Chin virial.
        dself.virs_4th_order = depend_array(name="vir_4th_order", value=np.zeros((self.nbeads, 3, 3), float),
                                            dependencies=[dself.fvir_4th_order],
                                            func=(lambda: self.fvir_4th_order[1]))

        dself.virssc_part_1 = depend_array(name="virssc_part_1", value=np.zeros((self.nbeads, 3, 3), float),
                                           dependencies=[dself.coeffsc_part_1, dself.virs],
                                           func=self.get_virssc_part_1)

        dself.virssc_part_2 = depend_array(name="virssc_part_2", value=np.zeros((self.nbeads, 3, 3), float),
                                           dependencies=[dself.coeffsc_part_2, dself.virs_4th_order],
                                           func=self.get_virssc_part_2)

        dself.virssc = depend_array(name="virssc", value=np.zeros((self.nbeads, 3, 3), float),
                                    dependencies=[dself.virssc_part_1, dself.virssc_part_2],
                                    func=self.get_virssc)

        dself.virsc = depend_value(name="potsc",
                                   dependencies=[dself.potssc],
                                   func=(lambda: np.sum(self.virssc, axis=0)))

        # Add dependencies from the force weights, that are applied here when the total
        # force is assembled from its components

        for fc in self.mforces:
            dself.f.add_dependency(dd(fc).weight)
            dself.pots.add_dependency(dd(fc).weight)
            dself.virs.add_dependency(dd(fc).weight)

    def copy(self, beads=None, cell=None):
        """ Returns a copy of this force object that can be used to compute forces,
        e.g. for use in internal loops of geometry optimizers, or for property
        calculation.

        Args:
           beads: Optionally, bind this to a different beads object than the one
              this Forces is currently bound
           cell: Optionally, bind this to a different cell object

        Returns: The copy of the Forces object
        """

        if not self.bound: raise ValueError("Cannot copy a forces object that has not yet been bound.")
        nforce = Forces()
        nbeads = beads
        if nbeads is None: nbeads = self.beads
        ncell = cell
        if cell is None: ncell = self.cell
        nforce.bind(nbeads, ncell, self.fcomp, self.ff, self.open_paths)
        return nforce

    def transfer_forces(self, refforce):
        """Low-level function copying over the value of a second force object,
        triggering updates but un-tainting this force depends themselves.

        We have noted that in some corner cases it is necessary to copy only
        the values of updated forces rather than the full depend object, in order to
        avoid triggering a repeated call to the client code that is potentially
        very costly. This happens routinely in geometry relaxation routines, for example."""

        if len(self.mforces) != len(refforce.mforces):
            raise ValueError("Cannot copy forces between objects with different numbers of components")

        for k in xrange(len(self.mforces)):
            mreff = refforce.mforces[k]
            mself = self.mforces[k]
            if mreff.nbeads != mself.nbeads:
                raise ValueError("Cannot copy forces between objects with different numbers of beads for the " + str(k) + "th component")

            # this is VERY subtle. beads in this force component are
            # obtained as a contraction, and so are computed automatically.
            # when we set the master q, these get marked as tainted.
            # then we copy the force value, and set the force as untainted.
            # next time we touch the master q, the tainting does not get
            # propagated, because the contracted q is already marked as tainted,
            # so the force does not get updated. we can fix this by copying
            # the value of the contracted bead, so that it's marked as NOT
            # tainted - it should not be as it's an internal of the force and
            # therefore get copied
            dd(mself.beads).q.set(mreff.beads.q, manual=False)
            for b in xrange(mself.nbeads):
                dfkbref = dd(mreff._forces[b])
                dfkbself = dd(mself._forces[b])

                dfkbself.ufvx.set(deepcopy(dfkbref.ufvx._value), manual=False)
                dfkbself.ufvx.taint(taintme=False)

    def run(self):
        """Makes the socket start looking for driver codes.

        Tells the interface code to start the thread that looks for
        connection from the driver codes in a loop. Until this point no
        jobs can be queued.
        """

        for ff in self.mforces:
            ff.run()

    def stop(self):
        """Makes the socket stop looking for driver codes.

        Tells the interface code to stop the thread that looks for
        connection from the driver codes in a loop. After this point no
        jobs can be queued.
        """

        for ff in self.mforces:
            ff.stop()

    def queue(self):
        """Submits all the required force calculations to the forcefields."""

        for ff in self.mforces:
            if ff.weight != 0:  # do not compute forces which have zero weight
                ff.queue()

    def get_vir(self):
        """Sums the virial of each forcefield.

        Not the actual system virial.

        Returns:
            Virial sum.
        """

        vir = np.zeros((3, 3))
        for v in dstrip(self.virs):
            vir += v
        return vir

    def pots_component(self, index, weighted=True):
        """Fetches the index^th component of the total potential."""
        if weighted:
            if self.mforces[index].weight != 0:
                return self.mforces[index].weight * self.mrpc[index].b2tob1(self.mforces[index].pots)
            else:
                return 0
        else:
            return self.mrpc[index].b2tob1(self.mforces[index].pots)

    def forces_component(self, index, weighted=True):
        """ Fetches the index^th component of the total force."""
        if weighted:
            if self.mforces[index].weight != 0:
                return self.mforces[index].weight * self.mrpc[index].b2tob1(dstrip(self.mforces[index].f))
            else:
                return np.zeros((self.nbeads, self.natoms * 3), float)
        else:
            return self.mrpc[index].b2tob1(dstrip(self.mforces[index].f))

    def queue_mts(self, level):
        """Submits all the required force calculations to the forcefields."""

        for ff in self.mforces:
            # forces with no MTS specification are applied at the outer level
            if ((len(ff.mts_weights) == 0 and level == 0) or
                (len(ff.mts_weights) > level
                 and ff.mts_weights[level] != 0
                 and ff.weight != 0)):
                # do not queue forces which have zero weight
                ff.queue()

    def forces_mts(self, level):
        """ Fetches ONLY the forces associated with a given MTS level."""

        self.queue_mts(level)
        fk = np.zeros((self.nbeads, 3 * self.natoms))
        for index in range(len(self.mforces)):
            # forces with no MTS specification are applied at the outer level
            if ((len(self.mforces[index].mts_weights) == 0 and level == 0) or
                (len(self.mforces[index].mts_weights) > level
                 and self.mforces[index].mts_weights[level] != 0
                 and self.mforces[index].weight != 0)):
                fk += self.mforces[index].weight * self.mforces[index].mts_weights[level] * self.mrpc[index].b2tob1(dstrip(self.mforces[index].f))
        return fk

    def forcesvirs_4th_order(self, index):
        """ Fetches the 4th order |f^2| correction to the force vector and the virial associated with a given component."""

        # gives an error is number of beads is not even.
        if self.nbeads % 2 != 0:
            warning("ERROR: Suzuki-Chin factorization requires even number of beads!")
            exit()

        # calculates the finite displacement.
        fbase = dstrip(self.f)
        eps = self.mforces[index].epsilon
        delta = np.abs(eps) / np.sqrt((fbase / self.beads.m3 * fbase / self.beads.m3).sum() / (self.nbeads * self.natoms))
        dq = delta * fbase / self.beads.m3

        # stores the force component.
        fbase = self.mrpc[index].b2tob1(dstrip(self.mforces[index].f))
        vbase = np.zeros((self.nbeads, 3, 3), float)
        mvirs = dstrip(self.mforces[index].virs)
        for i in range(3):
            for j in range(3):
                vbase[:, i, j] += self.mrpc[index].b2tob1(mvirs[:, i, j])

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
            # is marginal in when RPC + MTS is used.

            if self.mforces[index].nbeads != self.nbeads:
                warning("ERROR: high order PIMD + RPC works with a centered finite difference only! (Uness you find an elegant solution :))")
                exit()

            # for the case of alpha = 0, only odd beads are displaced.
            if self.alpha == 0:

                # we use an aux force evaluator with half the number of beads.
                if self.dforces is None:
                    self.dbeads = self.beads.copy(self.nbeads / 2)
                    self.dcell = self.cell.copy()
                    self.dforces = self.copy(self.dbeads, self.dcell)

                self.dcell.h = self.cell.h

                f_4th_order = fbase * 0.0
                v_4th_order = vbase * 0.0

                # displaces odd beads only.
                self.dbeads.q = dstrip(self.beads.q)[1::2] - dq[1::2]

                # calculates the force.
                fminus = self.dforces.mrpc[index].b2tob1(dstrip(self.dforces.mforces[index].f))

                # calculates the virial.
                vminus = np.zeros((self.nbeads / 2, 3, 3), float)
                dmvirs = dstrip(self.dforces.mforces[index].virs)
                for i in range(3):
                    for j in range(3):
                        vminus[:, i, j] += self.dforces.mrpc[index].b2tob1(dmvirs[:, i, j])

                # calculates the finite difference.
                f_4th_order[1::2] = 2.0 * (fminus - fbase[1::2]) / delta
                v_4th_order[1::2] = 2.0 * (vminus - vbase[1::2]) / delta

            # For the case of alpha != 0, all the beads are displaced.
            else:

                # we use an aux force evaluator with the same number of beads.
                if self.dforces is None:
                    self.dbeads = self.beads.copy()
                    self.dcell = self.cell.copy()
                    self.dforces = self.copy(self.dbeads, self.dcell)

                self.dcell.h = self.cell.h

                f_4th_order = fbase * 0.0
                v_4th_order = vbase * 0.0

                # displaces the beads.
                self.dbeads.q = self.beads.q + dq

                # calculates the force.
                fplus = self.dforces.mrpc[index].b2tob1((dstrip(self.dforces.mforces[index].f)))

                # calculates the virial.
                vplus = np.zeros((self.nbeads, 3, 3), float)
                dmvirs = dstrip(self.dforces.mforces[index].virs)
                for i in range(3):
                    for j in range(3):
                        vplus[:, i, j] += self.dforces.mrpc[index].b2tob1(dmvirs[:, i, j])

                # calculates the finite difference.
                f_4th_order = 2.0 * (fbase - fplus) / delta
                v_4th_order = 2.0 * (vbase - vplus) / delta

        # uses a centered difference for epsilon  < 0.
        if self.mforces[index].epsilon < 0.0:

                # we use an aux force evaluator with the same number of beads.
            if self.dforces is None:
                self.dbeads = self.beads.copy()
                self.dcell = self.cell.copy()
                self.dforces = self.copy(self.dbeads, self.dcell)

            self.dcell.h = self.cell.h

            f_4th_order = fbase * 0.0
            v_4th_order = vbase * 0.0

            # for the case of alpha = 0, only odd beads are displaced.
            if self.alpha == 0:

                # the first half of the aux beads are fwd displaced while the second half are bkwd displaced configurations.
                self.dbeads.q[:self.nbeads / 2] = dstrip(self.beads.q)[1::2] + dq[1::2]
                self.dbeads.q[-self.nbeads / 2:] = dstrip(self.beads.q)[1::2] - dq[1::2]

                # calculates the forces.
                fplusminus = self.dforces.mrpc[index].b2tob1(dstrip(self.dforces.mforces[index].f))

                # calculates the virial.
                vplusminus = np.zeros((self.nbeads, 3, 3), float)
                dmvirs = dstrip(self.dforces.mforces[index].virs)
                for i in range(3):
                    for j in range(3):
                        vplusminus[:, i, j] += self.dforces.mrpc[index].b2tob1(dmvirs[:, i, j])

                # calculates the finite difference.
                for k in range(self.nbeads / 2):
                    j = 2 * k + 1
                    f_4th_order[j] = 2.0 * (fplusminus[self.nbeads / 2 + k] - fplusminus[k]) / 2.0 / delta
                    v_4th_order[j] = 2.0 * (vplusminus[self.nbeads / 2 + k] - vplusminus[k]) / 2.0 / delta

            # For the case of alpha != 0, all the beads are displaced.
            else:
                # displaces the beads.
                self.dbeads.q = self.beads.q + dq

                # calculates the forces.
                fplus = self.dforces.mrpc[index].b2tob1(dstrip(self.dforces.mforces[index].f))

                # calculates the virial.
                vplus = np.zeros((self.nbeads, 3, 3), float)
                dmvirs = dstrip(self.dforces.mforces[index].virs)
                for i in range(3):
                    for j in range(3):
                        vplus[:, i, j] += self.dforces.mrpc[index].b2tob1(dmvirs[:, i, j])

                # displaces the beads.
                self.dbeads.q = self.beads.q - dq

                # calculates the forces.
                fminus = self.dforces.mrpc[index].b2tob1(dstrip(self.dforces.mforces[index].f))

                # calculates the virial.
                vminus = np.zeros((self.nbeads, 3, 3), float)
                dmvirs = dstrip(self.dforces.mforces[index].virs)
                for i in range(3):
                    for j in range(3):
                        vminus[:, i, j] += self.dforces.mrpc[index].b2tob1(dmvirs[:, i, j])

                # calculates the finite difference.
                f_4th_order = 2.0 * (fminus - fplus) / 2.0 / delta
                v_4th_order = 2.0 * (vminus - vplus) / 2.0 / delta

        # returns the 4th order |f^2| correction.
        return [f_4th_order, v_4th_order]

    def vir_mts(self, level):
        """ Fetches ONLY the total virial associated with a given MTS level."""
        return np.sum(self.virs_mts(level), axis=0)

    def virs_mts(self, level):
        """ Fetches ONLY the total virial associated with a given MTS level."""

        self.queue_mts(level)
        rp = np.zeros((self.beads.nbeads, 3, 3), float)
        for index in range(len(self.mforces)):
            if (len(self.mforces[index].mts_weights) > level
                and self.mforces[index].mts_weights[level] != 0 and
                self.mforces[index].weight != 0):
                dv = np.zeros((self.beads.nbeads, 3, 3), float)
                for i in range(3):
                    for j in range(3):
                        dv[:, i, j] += self.mrpc[index].b2tob1(self.mforces[index].virs[:, i, j])
                rp += self.mforces[index].weight * self.mforces[index].mts_weights[level] * dv
        return rp

    def get_nmtslevels(self):
        """ Returns the total number of mts levels."""

        nm = len(self.mforces[0].mts_weights)
        if all(len(x.mts_weights) == nm for x in self.mforces):
            return nm
        else:
            raise ValueError("The mts_weights of all the force components are not the same.")

    def f_combine(self):
        """Obtains the total force vector."""

        self.queue()
        rf = np.zeros((self.nbeads, 3 * self.natoms), float)
        for k in range(self.nforces):
            # "expand" to the total number of beads the forces from the
            # contracted one
            if self.mforces[k].weight != 0:
                rf += self.mforces[k].weight * self.mforces[k].mts_weights.sum() * self.mrpc[k].b2tob1(dstrip(self.mforces[k].f))
        return rf

    def fvir_4th_order_combine(self):
        """Obtains the total fourth order |f^2| correction to the force vector and the virial."""

        rf = np.zeros((self.nbeads, 3 * self.natoms), float)
        rv = np.zeros((self.nbeads, 3, 3), float)

        for k in range(self.nforces):
            if self.mforces[k].weight != 0 and self.mforces[k].mts_weights.sum() != 0:
                fv = dstrip(self.forcesvirs_4th_order(k))
                rf += self.mforces[k].weight * self.mforces[k].mts_weights.sum() * fv[0]
                rv += self.mforces[k].weight * self.mforces[k].mts_weights.sum() * fv[1]
        return [rf, rv]

    def pot_combine(self):
        """Obtains the potential energy for each forcefield."""

        self.queue()
        rp = np.zeros(self.nbeads, float)
        for k in range(self.nforces):
            # "expand" to the total number of beads the potentials from the
            # contracted one
            if self.mforces[k].weight != 0:
                rp += self.mforces[k].weight * self.mforces[k].mts_weights.sum() * self.mrpc[k].b2tob1(self.mforces[k].pots)
        return rp

    def extra_combine(self):
        """Obtains the potential energy for each forcefield."""

        self.queue()
        rp = ["" for b in range(self.nbeads)]
        for k in range(self.nforces):
            # "expand" to the total number of beads the potentials from the
            # contracted one
            for b in range(self.nbeads):
                rp[b] += self.mforces[k].extras[b]
        return rp

    def vir_combine(self):
        """Obtains the virial tensor for each forcefield."""

        self.queue()
        rp = np.zeros((self.nbeads, 3, 3), float)
        for k in range(self.nforces):
            if self.mforces[k].weight != 0:
                virs = dstrip(self.mforces[k].virs)
                # "expand" to the total number of beads the virials from the
                # contracted one, element by element
                for i in range(3):
                    for j in range(3):
                        rp[:, i, j] += self.mforces[k].weight * self.mforces[k].mts_weights.sum() * self.mrpc[k].b2tob1(virs[:, i, j])
        return rp

    def get_potssc(self):
        """Obtains Suzuki-Chin contribution to the potential."""
        if self.nbeads % 2 != 0:
            warning("ERROR: Suzuki-Chin factorization requires even number of beads!")
            exit()

        # this evaluates the square forces contribution to the SC potential (only the difference with the Trotter potential is returned)
        return self.coeffsc_part_1.T * dstrip(self.pots) + self.coeffsc_part_2.T * np.sum(dstrip(self.f) / self.beads.m3 * dstrip(self.f), axis=1)

    def get_fsc_part_1(self):
        """Obtains the linear component of Suzuki-Chin correction to the force."""

        return self.coeffsc_part_1 * dstrip(self.f)

    def get_fsc_part_2(self):
        """Obtains the quadratic component of Suzuki-Chin correction to the force."""

        return self.coeffsc_part_2 * dstrip(self.f_4th_order)

    def get_virssc_part_1(self):
        """Obtains the linear component of Suzuki-Chin correction to the virial."""
        return self.coeffsc_part_1.reshape((self.beads.nbeads, 1, 1)) * dstrip(self.virs)

    def get_virssc_part_2(self):
        """Obtains the quadratic component of Suzuki-Chin correction to the virial."""
        return self.coeffsc_part_2.reshape((self.beads.nbeads, 1, 1)) * dstrip(self.virs_4th_order)

    def get_fsc(self):
        """Obtains the total Suzuki-Chin correction to the force."""
        return dstrip(self.fsc_part_1) + dstrip(self.fsc_part_2)

    def get_virssc(self):
        """Obtains the total Suzuki-Chin correction to the force."""
        return dstrip(self.virssc_part_1) + dstrip(self.virssc_part_2)

    def get_coeffsc_part_1(self):
        """Obtains the coefficients of the linear part of the Suzuki-Chin correction."""

        rc = np.zeros(self.beads.nbeads)
        rc[0::2] = -1.0 / 3.0
        rc[1::2] = 1.0 / 3.0
        return np.asmatrix(rc).T

    def get_coeffsc_part_2(self):
        """Obtains the coefficients of the linear part of the Suzuki-Chin correction."""

        rc = np.zeros(self.beads.nbeads)
        rc[0::2] = (self.alpha / self.omegan2 / 9.0)
        rc[1::2] = ((1.0 - self.alpha) / self.omegan2 / 9.0)
        return np.asmatrix(rc).T

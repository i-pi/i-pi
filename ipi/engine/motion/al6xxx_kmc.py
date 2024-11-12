"""Highly specialized Kinetic Monte Carlo class for aluminum 6xxx alloys.
Could probably be generalized to something more general.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import pickle
import threading
import numpy as np
import collections

from ipi.engine.motion import Motion, GeopMotion
from ipi.utils.depend import dstrip, depend_value
from ipi.engine.cell import Cell
from ipi.utils.units import Constants
import ipi.utils.io as io


class AlKMC(Motion):
    """Stepper for a KMC for Al-6xxx alloys.

    Gives the standard methods and attributes needed in all the
    dynamics classes.

    Attributes:
        beads: A beads object giving the atoms positions.
        cell: A cell object giving the system box.
        forces: A forces object giving the virial and the forces acting on
            each bead.
        prng: A random number generator object.
        nm: An object which does the normal modes transformation.

    Depend objects:
        econs: The conserved energy quantity appropriate to the given
            ensemble. Depends on the various energy terms which make it up,
            which are different depending on the ensemble.he
        temp: The system temperature.
        dt: The timestep for the algorithms.
        ntemp: The simulation temperature. Will be nbeads times higher than
            the system temperature as PIMD calculations are done at this
            effective classical temperature.
    """

    def __init__(
        self,
        mode,
        optimizer,
        nstep,
        a0,
        ncell,
        nvac,
        nsi,
        nmg,
        neval,
        diffusion_barrier_al,
        diffusion_prefactor_al,
        diffusion_barrier_mg,
        diffusion_prefactor_mg,
        diffusion_barrier_si,
        diffusion_prefactor_si,
        idx=[],
        tottime=0,
        ecache_file="",
        qcache_file="",
        thermostat=None,
        barostat=None,
        fixcom=False,
        fixatoms_dof=None,
        nmts=None,
        max_cache_len=1000,
    ):
        """Initialises a "dynamics" motion object.

        Args:
            dt: The timestep of the simulation algorithms.
            fixcom: An optional boolean which decides whether the centre of mass
                motion will be constrained or not. Defaults to False.
        """

        # This will generate a lattice model based on a primitive FCC cell. the lattice is represented in three ways:
        # 1. as a string in which each lattice site is identified by a letter
        # 2. by a list of the lattice sites in 3D space, in 1-1 mapping with the letters
        # 3. by a list of the atoms, whose lattice position is indicated by an integer

        self.nstep = nstep
        self.ncell = ncell
        self.nvac = nvac
        self.nsi = nsi
        self.nmg = nmg
        self.nsites = self.ncell**3
        self.natoms = self.nsites - self.nvac
        self.neval = neval
        self.diffusion_barrier_al = diffusion_barrier_al
        self.diffusion_prefactor_al = diffusion_prefactor_al
        if diffusion_barrier_mg > 0:
            self.diffusion_barrier_mg = diffusion_barrier_mg
        else:
            self.diffusion_barrier_mg = diffusion_barrier_al
        if diffusion_barrier_si > 0:
            self.diffusion_barrier_si = diffusion_barrier_si
        else:
            self.diffusion_barrier_si = diffusion_barrier_al
        if diffusion_prefactor_mg > 0:
            self.diffusion_prefactor_mg = diffusion_prefactor_mg
        else:
            self.diffusion_prefactor_mg = diffusion_prefactor_al
        if diffusion_prefactor_si > 0:
            self.diffusion_prefactor_si = diffusion_prefactor_si
        else:
            self.diffusion_prefactor_si = diffusion_prefactor_al
        self.barriers = {
            "A": self.diffusion_barrier_al,
            "M": self.diffusion_barrier_mg,
            "S": self.diffusion_barrier_si,
        }
        self.prefactors = {
            "A": self.diffusion_prefactor_al,
            "M": self.diffusion_prefactor_mg,
            "S": self.diffusion_prefactor_si,
        }

        self.a0 = a0
        cell = np.zeros((3, 3))
        cell[0] = [0.7071067811865475, 0.35355339059327373, 0.35355339059327373]
        cell[1] = [0.0, 0.6123724356957945, 0.20412414523193154]
        cell[2] = [0.0, 0.0, 0.5773502691896258]
        self.scell = self.a0 * cell
        self.dcell = Cell()
        self.dcell.h = self.scell * self.ncell

        print("LATTICE PARAM ", self.a0)
        # this is the list of lattice sites, in 3D coordinates
        ix, iy, iz = np.meshgrid(
            list(range(self.ncell)),
            list(range(self.ncell)),
            list(range(self.ncell)),
            indexing="ij",
        )
        self.sites = np.dot(
            np.asarray([ix.flatten(), iy.flatten(), iz.flatten()]).T, self.scell.T
        )
        print(len(self.sites), self.nsites, "###")
        # now we build list of nearest neighbors (fcc-lattice hardcoded!)
        self.neigh = np.zeros((self.nsites, 12), int)
        nneigh = np.zeros(self.nsites, int)
        # could be done in a more analytic way but whatever, I'm too lazy
        a02 = 1.01 * 0.5 * self.a0**2  # perhaps 1.01 it is not enough, must check!
        for i in range(self.nsites):  # determines the connectivity of the lattice
            rij = self.sites.copy().flatten()
            for j in range(self.nsites):
                rij[3 * j : 3 * j + 3] -= self.sites[i]
            self.dcell.array_pbc(rij)
            rij.shape = (self.nsites, 3)
            for j in range(i):
                if np.dot(rij[j], rij[j]) < a02:  # found nearest neighbor
                    self.neigh[i, nneigh[i]] = j
                    self.neigh[j, nneigh[j]] = i
                    nneigh[i] += 1
                    nneigh[j] += 1

        self.idx = idx

        # the KMC step is variable and so it cannot be stored as proper timing
        self._dt = depend_value(name="dt", value=0.0)
        self.fixatoms_dof = np.asarray([])
        self.fixcom = True
        self.optimizer = [None] * self.neval
        # geop should not trigger exit if there is early convergence, but just carry on.
        # we hard-code this option to avoid early-termination that would be hard to debug for a user
        optimizer["exit_on_convergence"] = False
        for i in range(self.neval):
            # geometry optimizer should not have *any* hystory dependence
            self.optimizer[i] = GeopMotion(
                fixcom=fixcom, fixatoms_dof=fixatoms_dof, **optimizer
            )  # mode="cg", ls_options={"tolerance": 1, "iter": 20,  "step": 1e-3, "adaptive": 0.0}, tolerances={"energy": 1e-7, "force": 1e-2, "position": 1e-4}, ) #!TODO: set the geop parameters properly

        # dictionary of previous energy evaluations - kind of tricky to use this with the omaker thingie
        self.ecache_file = ecache_file
        self.qcache_file = qcache_file
        self.max_cache_len = max_cache_len  # default 1000; user modification allowed
        try:
            ff = open(self.ecache_file, "rb")
            self.ecache = pickle.load(ff)
            ff.close()
            ff = open(self.qcache_file, "rb")
            self.qcache = pickle.load(ff)
            ff.close()
            print("Loaded %d cached energies" % (len(self.ecache)))
        except (OSError, ValueError, NameError):
            print(
                "Couldn't load cache files "
                + self.ecache_file
                + ","
                + self.qcache_file
                + " - resetting"
            )
            self.ecache = collections.OrderedDict()
            self.qcache = collections.OrderedDict()
        self.ncache = len(self.ecache)
        self.ncache_stored = self.ncache
        self.struct_count = self.ncache

        # no TS evaluation implemented yet
        self.tscache = {}  # collections.OrderedDict()
        self.tottime = tottime

    def bind(self, ens, beads, nm, cell, bforce, prng, omaker):
        """Binds ensemble beads, cell, bforce, and prng to the dynamics.

        This takes a beads object, a cell object, a forcefield object and a
        random number generator object and makes them members of the ensemble.
        It also then creates the objects that will hold the data needed in the
        ensemble algorithms and the dependency network. Note that the conserved
        quantity is defined in the init, but as each ensemble has a different
        conserved quantity the dependencies are defined in bind.

        Args:
            beads: The beads object from whcih the bead positions are taken.
            nm: A normal modes object used to do the normal modes transformation.
            cell: The cell object from which the system box is taken.
            bforce: The forcefield object from which the force and virial are
                taken.
            prng: The random number generator object which controls random number
                generation.
        """

        super(AlKMC, self).bind(ens, beads, nm, cell, bforce, prng, omaker)

        # todo make these optional and restarted
        self.kmcfile = self.output_maker.get_output("KMC_AL6XXX")

        # this is the index for the atoms, self.idx[i] indicates the lattice site of atom i.
        # atoms are in the order Si1 Si2 ... Mg1 Mg2 .... Al1 Al2 .... Vac1 Vac2 ...
        f_restart = True
        if self.idx is None or len(self.idx) == 0:
            f_restart = False
            idx = np.asarray(
                list(range(self.ncell**3)), int
            )  # initialize random distribution of atoms
            self.prng.shuffle(idx)
            self.idx = idx

        # initialize state based on the index
        # this is a string indicating the state of the system. each character corresponds to a lattice site
        # and can be S=Si, M=Mg, A=Al, V=Vac.
        # create a random string
        names = np.asarray(
            list(
                "S" * self.nsi
                + "M" * self.nmg
                + (self.natoms - self.nsi - self.nmg) * "A"
                + "V" * self.nvac
            )
        )
        state = names.copy()

        # this maps the string to random sites
        state[self.idx] = names  # backshuffle!
        state = "".join(state)

        # reverse lookup index [i.e. ridx[i] gives the index of the atoms at site i]
        self.ridx = np.zeros(self.nsites, int)
        self.ridx[self.idx] = list(range(self.nsites))

        self.state = np.asarray(list(state))
        print("".join(self.state))

        print("CHECKING INITIAL ASSIGNMENTS")
        for i in range(self.nsites):
            if self.ridx[i] < self.natoms:
                print(self.beads.names[self.ridx[i]], self.state[i])
            else:
                print("V", self.state[i])
            if self.idx[self.ridx[i]] != i:
                print(
                    "inconsistent site string for atom ", i, " and site ", self.ridx[i]
                )

        if not f_restart:
            self.beads.q[0, :] = self.sites[
                self.idx
            ].flatten()  # also sets global q so we can use it further down

        self.dbeads = [None] * self.neval
        self.dforces = [None] * self.neval
        self.dnm = [None] * self.neval
        self.dens = [None] * self.neval
        self.dbias = [None] * self.neval
        for i in range(self.neval):
            self.dbeads[i] = beads.clone()
            self.dforces[i] = bforce.copy(self.dbeads[i], self.dcell)
            self.dnm[i] = nm.copy()
            self.dens[i] = ens.copy()
            self.dnm[i].bind(
                self.dens[i], self, beads=self.dbeads[i], forces=self.dforces[i]
            )
            self.dbias[i] = ens.bias.copy(self.dbeads[i], self.dcell)
            self.dens[i].bind(
                self.dbeads[i],
                self.dnm[i],
                self.dcell,
                self.dforces[i],
                self.dbias[i],
                output_maker=self.output_maker,
            )
            self.optimizer[i].bind(
                self.dens[i],
                self.dbeads[i],
                self.dnm[i],
                self.dcell,
                self.dforces[i],
                prng,
                omaker,
            )
        self.feval = np.ones(self.neval, int)
        self._threadlock = threading.Lock()

    # threaded geometry optimization
    def geop_thread(self, ieval, nstr, nevent, ostr=None):
        self.optimizer[ieval].reset()
        ipot = self.dforces[ieval].pot

        for i in range(self.nstep):
            # print "geop ", i, self.dforces[ieval].pot
            self.optimizer[ieval].step(i)
        newq = dstrip(self.dbeads[ieval].q[0]).copy()
        newpot = self.dforces[ieval].pot

        # print "geop ", self.nstep, self.dforces[ieval].pot
        with self._threadlock:
            self.ecache[nstr] = newpot
            self.qcache[nstr] = newq
            if self.ncache == self.max_cache_len:
                self.ecache.popitem(last=False)
                self.qcache.popitem(last=False)
            else:
                self.ncache += 1
            nevent[2] = self.ecache[nstr]
            nevent[3] = self.qcache[nstr]
            self.struct_count += 1

        # launches TS calculation
        # if not ostr is None:
        #    self.ts_thread(ieval, ostr, nstr, nevent)
        # else:
        #    with self._threadlock:
        #        self.feval[ieval] = 1
        with self._threadlock:
            self.feval[ieval] = 1

        with self._threadlock:
            print("Finished ", nstr)
            print("Energy, initial - TS - final: ", ipot, nevent[-1], newpot)

    # threaded ts evaluation
    def ts_thread(self, ieval, ostr, nstr, nevent, setfev=1):
        # computes TS energy by linearly interpolating initial & final state
        # interpolates between two structures considering PBCs
        qstart = self.qcache[ostr]
        qend = self.qcache[nstr].copy()
        # finds atom matching assuming initial and final states differ only by a vacancy swap and are based on unique_idx lists
        midx = self.unique_idx_match(list(ostr), list(nstr))[0 : self.natoms]
        qend.shape = (self.natoms, 3)
        qend[:] = qend[midx]
        qend.shape = 3 * self.natoms

        qts = qend - qstart
        self.dcell.array_pbc(qts)  # use pbc!
        qts *= 0.5
        qts += qstart
        self.dbeads[ieval].q[0, :] = qts

        tspot = self.dforces[ieval].pot
        io.print_file(
            "xyz",
            self.beads[0],
            self.dcell,
            self.tslist,
            title=("START  "),
            key="positions",
            dimension="length",
            units="angstrom",
            cell_units="angstrom",
        )
        io.print_file(
            "xyz",
            self.dbeads[ieval][0],
            self.dcell,
            self.tslist,
            title=(
                "TS: %s %s  Energy: %15.8e  %15.8e "
                % (ostr, nstr, tspot, tspot - self.forces.pot)
            ),
            key="positions",
            dimension="length",
            units="angstrom",
            cell_units="angstrom",
        )
        self.dbeads[ieval].q[0] = qend
        io.print_file(
            "xyz",
            self.dbeads[ieval][0],
            self.dcell,
            self.tslist,
            title=("END  "),
            key="positions",
            dimension="length",
            units="angstrom",
            cell_units="angstrom",
        )
        # self.tslist.flush()

        with self._threadlock:
            # sets the tscache for both FW and BW transition (uses a dictionary to be super-safe and lazy, although of course this could be just a list)
            self.tscache[ostr][nstr] = tspot
            if nstr not in self.tscache:
                self.tscache[nstr] = {}
            self.tscache[nstr][ostr] = tspot
            self.feval[ieval] = 1
            nevent[4] = tspot

    def find_eval(self, ethreads):
        # finds first free evaluator
        while (
            self.feval.sum() == 0
        ):  # if all evaluators are busy, wait for one to get free
            for st in ethreads:
                st.join(1e-2)
                if st is None or not st.is_alive():
                    break
        with self._threadlock:
            # finds free evaluator
            for e in range(self.neval):
                if self.feval[e] == 1:
                    ieval = e
                    self.feval[ieval] = 0
                    break
        return ieval

    def unique_idx(self, state):
        # generates a starting lattice configuration that corresponds to a given state vector (a string of atomic types)
        # makes sure that the same occupation string corresponds to the same atom positions,
        # even though the actual atom indices might be different basically, makes the configurations
        # independent of same-atom permutations
        ksi = 0
        kmg = self.nsi
        kal = self.nsi + self.nmg
        kvac = self.natoms
        k = 0
        idx = np.zeros(self.nsites, int)
        for s in state:
            if s == "S":
                idx[ksi] = k
                ksi += 1
            elif s == "M":
                idx[kmg] = k
                kmg += 1
            elif s == "A":
                idx[kal] = k
                kal += 1
            elif s == "V":
                idx[kvac] = k
                kvac += 1
            k += 1
        return idx

    def unique_idx_match(self, state_1, state_2):
        # finds the best matching between the atoms in states state_1 and state_2,
        # assuming there is only one site swap
        # (should raise an error otherwise but it's not trivial without doing too many extra checks)
        # this basically says what is the motion of atoms from state 1 to state 2. This is useful
        # to interpolate between atomic coordinates in the two states, e.g. to get the TS geometry

        # gets the unique mapping of atoms from the state to the site ids
        uid_1 = self.unique_idx(state_1)
        ru1 = np.zeros(self.nsites, int)
        # this is the reverse map. what is the atom index that sits in a given site?
        ru1[uid_1] = np.asarray(list(range(self.nsites)), int)

        uid_2 = self.unique_idx(state_2)
        ru2 = np.zeros(self.nsites, int)
        ru2[uid_2] = np.asarray(
            list(range(self.nsites)), int
        )  # this says which atom is in a given site in u2

        iu12 = ru2[uid_1]
        iu21 = ru1[uid_2]
        for i in range(self.natoms):
            if iu12[i] >= self.natoms:
                # print "found vacancy swap 1->2", i, u1[i], u2[iu12[i]], iu12[i]
                i1vac2 = i
            if iu21[i] >= self.natoms:
                i2vac1 = i
                # print "found vacancy swap 2->1", i, u2[i], u1[iu21[i]], iu21[i]
        iu12[i1vac2] = i2vac1

        return iu12

    def step(self, step=None):
        kT = Constants.kb * self.ensemble.temp
        # computes current energy (if not already stored)
        ostr = "".join(
            self.state
        )  # this is a unique id string that charactrizes the current state
        self.tscache[ostr] = {}
        if ostr not in self.ecache:
            self.dbeads[0].q[0, :] = self.sites[self.unique_idx(self.state)].flatten()
            rv = [0, 0, 0, 0, 0]
            self.geop_thread(0, ostr, rv)
            # self.beads.q[0,:] = self.dbeads[0].q[0,:] # also updates current position
            # self.forces.transfer_forces(self.dforces[0]) # forces have already been computed here...

        ecurr = self.ecache[ostr]

        # enumerates possible reactive events (vacancy swaps)
        levents = []
        ethreads = [None] * self.neval
        # loops over the vacancy
        for ivac in range(self.natoms, self.natoms + self.nvac):
            svac = self.idx[ivac]  # lattice site associated with this vacancy
            if self.state[svac] != "V":
                raise IndexError(
                    "Something got screwed and a vacancy state is not a vacancy anymore!"
                )
            # loops over the neighbors of the selected vacancy
            for sneigh in self.neigh[svac]:
                # if the neighbor is a vacancy, move on. does not make sense to swap two vacancies!
                if self.state[sneigh] == "V":
                    continue

                # creates a new state vector with swapped atoms-vacancy and the associated label string
                nstate = self.state.copy()
                nstate[svac], nstate[sneigh] = self.state[sneigh], self.state[svac]
                nstr = "".join(
                    nstate
                )  # this is the string that corresponds to the new state
                if nstr not in self.ecache:
                    # new state, must compute!
                    # creates a swapped index
                    nidx = self.idx.copy()
                    if self.ridx[svac] != ivac or self.idx[ivac] != svac:
                        raise IndexError(
                            "Something got screwed and the index does not correspond anymore to site occupancy"
                        )

                    # ivac = self.ridx[svac]
                    ineigh = self.ridx[
                        sneigh
                    ]  # gets index of atom associated with the neighboring site
                    nidx[ivac], nidx[ineigh] = self.idx[ineigh], self.idx[ivac]

                    ieval = self.find_eval(ethreads)
                    # launches evaluator
                    self.dbeads[ieval].q[0, :] = self.sites[
                        self.unique_idx(nstate)
                    ].flatten()

                    nevent = [svac, sneigh, 0.0, 0.0, 0.0]

                    # runs a geometry optimization
                    # self.geop_thread(ieval=ieval, nstr=nstr, nevent=nevent)
                    st = threading.Thread(
                        target=self.geop_thread,
                        name=str(ieval),
                        kwargs={
                            "ieval": ieval,
                            "nstr": nstr,
                            "nevent": nevent,
                            "ostr": ostr,
                        },
                    )
                    st.daemon = True
                    st.start()
                    ethreads[ieval] = st
                else:
                    print(
                        "Found state ",
                        nstr,
                        " retrieving cached energy ",
                        self.ecache[nstr],
                    )

                    # fetch energy from previous calculation
                    nevent = [svac, sneigh, self.ecache[nstr], self.qcache[nstr], 0.0]

                    # EVALUATION OF TS ENERGY IS DISABLED FOR THE MOMENT...
                    # we might still need to compute the TS energy!
                    # if not nstr in self.tscache[ostr]:
                    #   print "Computing TS"
                    #   ieval = self.find_eval(ethreads)
                    #    st = threading.Thread(target=self.ts_thread, name=str(ieval), kwargs={"ieval":ieval, "ostr": ostr, "nstr":nstr, "nevent" : nevent})
                    #    st.daemon = True
                    #    st.start()
                    #    ethreads[ieval] = st
                    # else:
                    #    print "Found TS"
                    #    nevent[3] = self.tscache[ostr][nstr]
                nevent[-1] = self.state[sneigh]
                levents.append(nevent)
        # wait for all evaluators to finish
        for st in ethreads:
            while st is not None and st.is_alive():
                st.join(2)

        print(
            "Computed ",
            len(levents),
            " possible reactions. Cache len ",
            len(self.ecache),
        )

        # get list of rates
        rates = np.zeros(len(levents), float)
        crates = np.zeros(len(levents), float)
        cdf = 0.0
        for i in range(len(levents)):
            # print ("Barrier, naive: %f, static: %f" % (0.5*(ecurr + levents[i][2]) + self.diffusion_barrier_al, levents[i][4]))

            ets = (
                0.5 * (ecurr + levents[i][2]) + self.barriers[levents[i][-1]]
            )  # naive heuristic for the ts energy
            print("Event ", i, levents[i][-1], ecurr, ">>", ets, ">>", levents[i][2])
            rates[i] = self.prefactors[levents[i][-1]] * np.exp(-(ets - ecurr) / kT)
            cdf += rates[i]
            crates[i] = cdf

        # KMC selection rule based on the rate
        fpick = self.prng.u * cdf
        isel = 0
        while fpick > crates[isel]:
            isel += 1
        dt = -1.0 / cdf * np.log(1.0 - self.prng.u)
        print(("Time spent %12.5e at %s nrg %12.5e" % (dt, ostr, ecurr)))
        print("Selected event ", isel, " with rate ", rates[isel], " / ", cdf)

        iev = levents[isel]  # levents[self.prng.randint(len(levents))]
        svac, sneigh = iev[0], iev[1]
        ivac, ineigh = self.ridx[svac], self.ridx[sneigh]

        # does the swap (never reject, for the moment)
        self.state[svac], self.state[sneigh] = self.state[sneigh], self.state[svac]
        self.ridx[svac], self.ridx[sneigh] = self.ridx[sneigh], self.ridx[svac]
        self.idx[ivac], self.idx[ineigh] = self.idx[ineigh], self.idx[ivac]

        # we got a new configuration but the residence time is linked to the previous configuration so we output that
        self.kmcfile.write(
            "%12.5e  %12.5e  %18.11e  %s\n" % (self.tottime, dt, ecurr, ostr)
        )
        self.kmcfile.force_flush()
        self.tottime += dt
        self.ensemble.time += dt  # updates time counter
        print("Finishing step at ", "".join(self.state))

        # updates the positions
        self.cell.h = self.dcell.h

        uidx = self.unique_idx(self.state)
        ruidx = np.zeros(self.nsites, int)
        ruidx[uidx] = list(range(self.nsites))

        self.sites[self.unique_idx(self.state)]
        oldq = dstrip(self.beads.q[0]).copy()

        newq = np.zeros(self.nsites * 3, float)
        # we want continuity (modulo PBC jumps, that we'll take care of later...)
        for i in range(self.nsites):
            # in which site sits atom i?
            # isite = self.idx[i]       #

            # which atom sits in this site in the unique-mapped structure?
            iuid = ruidx[self.idx[i]]
            newq[3 * i : 3 * i + 3] = iev[3][3 * iuid : 3 * (iuid + 1)]
        newq -= oldq
        self.cell.array_pbc(newq)
        self.beads.q[0] += newq

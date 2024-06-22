"""Contains the classes that connect the driver to the python code.

ForceField objects are force providers, i.e. they are the abstraction
layer for a driver that gets positions and returns forces (and energy).
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import time
import threading
import json
import sys

import numpy as np

from ipi.utils.softexit import softexit
from ipi.utils.messages import info, verbosity, warning
from ipi.interfaces.sockets import InterfaceSocket
from ipi.utils.depend import dstrip
from ipi.utils.io import read_file
from ipi.utils.units import unit_to_internal
from ipi.utils.distance import vector_separation

try:
    import plumed
except ImportError:
    plumed = None


class ForceRequest(dict):
    """An extension of the standard Python dict class which only has a == b
    if a is b == True, rather than if the elements of a and b are identical.

    Standard dicts are checked for equality if elements have the same value.
    Here I only care if requests are instances of the very same object.
    This is useful for the `in` operator, which uses equality to test membership.
    """

    def __eq__(self, y):
        """Overwrites the standard equals function."""
        return self is y


class ForceField:
    """Base forcefield class.

    Gives the standard methods and quantities needed in all the forcefield
    classes.

    Attributes:
        pars: A dictionary of the parameters needed to initialize the forcefield.
            Of the form {'name1': value1, 'name2': value2, ... }.
        name: The name of the forcefield.
        latency: A float giving the number of seconds the socket will wait
            before updating the client list.
        offset: A float giving a constant value that is subtracted from the return
            value of the forcefield
        requests: A list of all the jobs to be given to the client codes.
        dopbc: A boolean giving whether or not to apply the periodic boundary
            conditions before sending the positions to the client code.
        _thread: The thread on which the socket polling loop is being run.
        _doloop: A list of booleans. Used to decide when to stop running the
            polling loop.
        _threadlock: Python handle used to lock the thread held in _thread.
    """

    def __init__(
        self,
        latency=1.0,
        offset=0.0,
        name="",
        pars=None,
        dopbc=False,
        active=np.array([-1]),
        threaded=False,
    ):
        """Initialises ForceField.

        Args:
            latency: The number of seconds the socket will wait before updating
                the client list.
            offset: A constant offset subtracted from the energy value given by the
                client.
            name: The name of the forcefield.
            pars: A dictionary used to initialize the forcefield, if required.
                Of the form {'name1': value1, 'name2': value2, ... }.
            dopbc: Decides whether or not to apply the periodic boundary conditions
                before sending the positions to the client code.
            active: Indexes of active atoms in this forcefield
        """

        if pars is None:
            self.pars = {}
        else:
            self.pars = pars

        self.name = name
        self.latency = latency
        self.offset = offset
        self.requests = []
        self.dopbc = dopbc
        self.active = active
        self.iactive = None
        self.threaded = threaded
        self._thread = None
        self._doloop = [False]
        self._threadlock = threading.Lock()

    def bind(self, output_maker=None):
        """Binds the FF, at present just to allow for
        managed output"""

        self.output_maker = output_maker

    def queue(self, atoms, cell, reqid=-1, template=None):
        """Adds a request.

        Note that the pars dictionary need to be sent as a string of a
        standard format so that the initialisation of the driver can be done.

        Args:
            atoms: An Atoms object giving the atom positions.
            cell: A Cell object giving the system box.
            pars: An optional dictionary giving the parameters to be sent to the
                driver for initialisation. Defaults to {}.
            reqid: An optional integer that identifies requests of the same type,
               e.g. the bead index
            template: a dict giving a base model for the request item -
               e.g. to add entries that are not needed for the base class execution

        Returns:
            A dict giving the status of the request of the form {'pos': An array
            giving the atom positions folded back into the unit cell,
            'cell': Cell object giving the system box, 'pars': parameter string,
            'result': holds the result as a list once the computation is done,
            'status': a string labelling the status of the calculation,
            'id': the id of the request, usually the bead number, 'start':
            the starting time for the calculation, used to check for timeouts.}.
        """

        par_str = " "

        if self.pars is not None:
            for k, v in list(self.pars.items()):
                par_str += k + " : " + str(v) + " , "
        else:
            par_str = " "

        pbcpos = dstrip(atoms.q).copy()

        # Indexes come from input in a per atom basis and we need to make a per atom-coordinate basis
        # Reformat indexes for full system (default) or piece of system
        # active atoms do not change but we only know how to build this array once we get the positions once
        if self.iactive is None:
            if self.active[0] == -1:
                activehere = np.arange(len(pbcpos))
            else:
                activehere = np.array(
                    [[3 * n, 3 * n + 1, 3 * n + 2] for n in self.active]
                )

            # Reassign active indexes in order to use them
            activehere = activehere.flatten()

            # Perform sanity check for active atoms
            if len(activehere) > len(pbcpos) or activehere[-1] > (len(pbcpos) - 1):
                raise ValueError("There are more active atoms than atoms!")

            self.iactive = activehere

        if self.dopbc:
            cell.array_pbc(pbcpos)

        if template is None:
            template = {}
        template.update(
            {
                "id": reqid,
                "pos": pbcpos,
                "active": self.iactive,
                "cell": (dstrip(cell.h).copy(), dstrip(cell.ih).copy()),
                "pars": par_str,
                "result": None,
                "status": "Queued",
                "start": -1,
                "t_queued": time.time(),
                "t_dispatched": 0,
                "t_finished": 0,
            }
        )

        newreq = ForceRequest(template)

        with self._threadlock:
            self.requests.append(newreq)

        if not self.threaded:
            self.poll()

        return newreq

    def poll(self):
        """Polls the forcefield object to check if it has finished."""

        with self._threadlock:
            for r in self.requests:
                if r["status"] == "Queued":
                    r["t_dispatched"] = time.time()
                    r["result"] = [
                        0.0 - self.offset,
                        np.zeros(len(r["pos"]), float),
                        np.zeros((3, 3), float),
                        {"raw": ""},
                    ]
                    r["status"] = "Done"
                    r["t_finished"] = time.time()

    def _poll_loop(self):
        """Polling loop.

        Loops over the different requests, checking to see when they have
        finished.
        """

        info(" @ForceField: Starting the polling thread main loop.", verbosity.low)
        while self._doloop[0]:
            time.sleep(self.latency)
            if len(self.requests) > 0:
                self.poll()

    def release(self, request):
        """Removes a request from the evaluation queue.

        Args:
            request: The id of the job to release.
        """

        """Frees up a request."""

        with self._threadlock:
            if request in self.requests:
                try:
                    self.requests.remove(request)
                except ValueError:
                    print("failed removing request", id(request), " ", end=" ")
                    print(
                        [id(r) for r in self.requests], "@", threading.currentThread()
                    )
                    raise

    def stop(self):
        """Dummy stop method."""

        self._doloop[0] = False
        for r in self.requests:
            r["status"] = "Exit"

    def start(self):
        """Spawns a new thread.

        Splits the main program into two threads, one that runs the polling loop
        which updates the client list, and one which gets the data.

        Raises:
            NameError: Raised if the polling thread already exists.
        """

        if self._thread is not None:
            raise NameError("Polling thread already started")

        if self.threaded:
            self._doloop[0] = True
            self._thread = threading.Thread(
                target=self._poll_loop, name="poll_" + self.name
            )
            self._thread.daemon = True
            self._thread.start()
            softexit.register_thread(self._thread, self._doloop)
        softexit.register_function(self.softexit)

    def softexit(self):
        """Takes care of cleaning up upon softexit"""

        self.stop()

    def update(self):
        """Makes updates to the potential that only need to be triggered
        upon completion of a time step."""

        pass


class FFSocket(ForceField):
    """Interface between the PIMD code and a socket for a single replica.

    Deals with an individual replica of the system, obtaining the potential
    force and virial appropriate to this system. Deals with the distribution of
    jobs to the interface.

    Attributes:
        socket: The interface object which contains the socket through which
            communication between the forcefield and the driver is done.
    """

    def __init__(
        self,
        latency=1.0,
        offset=0.0,
        name="",
        pars=None,
        dopbc=True,
        active=np.array([-1]),
        threaded=True,
        interface=None,
    ):
        """Initialises FFSocket.

        Args:
           latency: The number of seconds the socket will wait before updating
              the client list.
           name: The name of the forcefield.
           pars: A dictionary used to initialize the forcefield, if required.
              Of the form {'name1': value1, 'name2': value2, ... }.
           dopbc: Decides whether or not to apply the periodic boundary conditions
              before sending the positions to the client code.
           interface: The object used to create the socket used to interact
              with the client codes.
        """

        # a socket to the communication library is created or linked
        super(FFSocket, self).__init__(
            latency, offset, name, pars, dopbc, active, threaded
        )
        if interface is None:
            self.socket = InterfaceSocket()
        else:
            self.socket = interface
        self.socket.requests = self.requests
        self.socket.offset = self.offset

    def poll(self):
        """Function to check the status of the client calculations."""

        self.socket.poll()

    def start(self):
        """Spawns a new thread."""

        self.socket.open()
        super(FFSocket, self).start()

    def stop(self):
        """Closes the socket and the thread."""

        super(FFSocket, self).stop()
        if self._thread is not None:
            # must wait until loop has ended before closing the socket
            self._thread.join()
        self.socket.close()


class FFEval(ForceField):
    """General class for models that provide a self.evaluate(request)
    to compute the potential, force and virial.
    """

    def poll(self):
        """Polls the forcefield checking if there are requests that should
        be answered, and if necessary evaluates the associated forces and energy."""

        # We have to be thread-safe, as in multi-system mode this might get
        # called by many threads at once.
        with self._threadlock:
            for r in self.requests:
                if r["status"] == "Queued":
                    r["status"] = "Running"
                    r["t_dispatched"] = time.time()
                    self.evaluate(r)
                    r["result"][0] -= self.offset  # subtract constant offset

    def evaluate(self, request):
        request["result"] = [
            0.0,
            np.zeros(len(request["pos"]), float),
            np.zeros((3, 3), float),
            {"raw": ""},
        ]
        request["status"] = "Done"


class FFLennardJones(FFEval):
    """Basic fully pythonic force provider.

    Computes LJ interactions without minimum image convention, cutoffs or
    neighbour lists. Parallel evaluation with threads.

    Attributes:
        parameters: A dictionary of the parameters used by the driver. Of the
            form {'name': value}.
        requests: During the force calculation step this holds a dictionary
            containing the relevant data for determining the progress of the step.
            Of the form {'atoms': atoms, 'cell': cell, 'pars': parameters,
                         'status': status, 'result': result, 'id': bead id,
                         'start': starting time}.
    """

    def __init__(
        self,
        latency=1.0e-3,
        offset=0.0,
        name="",
        pars=None,
        dopbc=False,
        threaded=False,
    ):
        """Initialises FFLennardJones.

        Args:
           pars: Optional dictionary, giving the parameters needed by the driver.
        """

        # check input - PBCs are not implemented here
        if dopbc:
            raise ValueError(
                "Periodic boundary conditions are not supported by FFLennardJones."
            )

        # a socket to the communication library is created or linked
        super(FFLennardJones, self).__init__(
            latency, offset, name, pars, dopbc=dopbc, threaded=threaded
        )
        self.epsfour = float(self.pars["eps"]) * 4
        self.sixepsfour = 6 * self.epsfour
        self.sigma2 = float(self.pars["sigma"]) * float(self.pars["sigma"])

    def evaluate(self, r):
        """Just a silly function evaluating a non-cutoffed, non-pbc and
        non-neighbour list LJ potential."""

        q = r["pos"].reshape((-1, 3))
        nat = len(q)

        v = 0.0
        f = np.zeros(q.shape)
        for i in range(1, nat):
            dij = q[i] - q[:i]
            rij2 = (dij**2).sum(axis=1)

            x6 = (self.sigma2 / rij2) ** 3
            x12 = x6**2

            v += (x12 - x6).sum()
            dij *= (self.sixepsfour * (2.0 * x12 - x6) / rij2)[:, np.newaxis]
            f[i] += dij.sum(axis=0)
            f[:i] -= dij

        v *= self.epsfour

        r["result"] = [v, f.reshape(nat * 3), np.zeros((3, 3), float), {"raw": ""}]
        r["status"] = "Done"


class FFdmd(FFEval):
    """Pythonic force provider.

    Computes DMD forces as in Bowman, .., Brown JCP 2003 DOI: 10.1063/1.1578475. It is a time dependent potential.
    Here extended for periodic systems and for virial term calculation.

    Attributes:
        parameters: A dictionary of the parameters used by the driver. Of the
            form {'name': value}.
        requests: During the force calculation step this holds a dictionary
            containing the relevant data for determining the progress of the step.
            Of the form {'atoms': atoms, 'cell': cell, 'pars': parameters,
                         'status': status, 'result': result, 'id': bead id,
                         'start': starting time}.
    """

    def __init__(
        self,
        latency=1.0e-3,
        offset=0.0,
        name="",
        coupling=None,
        freq=0.0,
        dtdmd=0.0,
        dmdstep=0,
        pars=None,
        dopbc=False,
        threaded=False,
    ):
        """Initialises FFdmd.

        Args:
           pars: Optional dictionary, giving the parameters needed by the driver.
        """

        # a socket to the communication library is created or linked
        super(FFdmd, self).__init__(
            latency, offset, name, pars, dopbc=dopbc, threaded=threaded
        )

        if coupling is None:
            raise ValueError("Must provide the couplings for DMD.")
        if freq is None:
            raise ValueError(
                "Must provide a frequency for the periodically oscillating potential."
            )
        if dtdmd is None:
            raise ValueError(
                "Must provide a time step for the periodically oscillating potential."
            )
        self.coupling = coupling
        self.freq = freq
        self.dtdmd = dtdmd
        self.dmdstep = dmdstep

    def evaluate(self, r):
        """Evaluating dmd: pbc: YES,
        cutoff: NO, neighbour list: NO."""

        q = r["pos"].reshape((-1, 3))
        nat = len(q)
        cell_h, cell_ih = r["cell"]

        if len(self.coupling) != int(nat * (nat - 1) / 2):
            raise ValueError("Coupling matrix size mismatch")

        v = 0.0
        f = np.zeros(q.shape)
        vir = np.zeros((3, 3), float)
        # must think and check handling of time step
        periodic = np.sin(self.dmdstep * self.freq * self.dtdmd)
        # MR: the algorithm below has been benchmarked against explicit loop implementation
        for i in range(1, nat):
            # MR's first implementation:
            #            dij = q[i] - q[:i]
            #            rij = np.sqrt((dij ** 2).sum(axis=1))
            # KF's implementation:
            dij, rij = vector_separation(cell_h, cell_ih, q[i], q[:i])
            cij = self.coupling[i * (i - 1) // 2 : i * (i + 1) // 2]
            prefac = np.dot(
                cij, rij
            )  # for each i it has the distances to all indexes previous
            v += np.sum(prefac) * periodic
            nij = np.copy(dij)
            nij *= -(cij / rij)[:, np.newaxis]  # magic line...
            f[i] += nij.sum(axis=0) * periodic
            f[:i] -= nij * periodic  # everything symmetric
            # virial:
            fij = nij * periodic
            for j in range(i):
                for cart1 in range(3):
                    for cart2 in range(3):
                        vir[cart1][cart2] += fij[j][cart1] * dij[j][cart2]
            # MR 2021: The virial looks correct and produces stable NPT simulations. It was not bullet-proof benchmarked, though.
            #          Because this is "out of equilibrium" I still did not find a good benchmark. Change cell and look at variation of energy only for this term?

        r["result"] = [v, f.reshape(nat * 3), vir, ""]
        r["status"] = "Done"

    def dmd_update(self):
        """Updates time step when a full step is done. Can only be called after implementation goes into smotion mode..."""
        self.dmdstep += 1


class FFDebye(FFEval):
    """Debye crystal harmonic reference potential

    Computes a harmonic forcefield.

    Attributes:
       parameters: A dictionary of the parameters used by the driver. Of the
          form {'name': value}.
       requests: During the force calculation step this holds a dictionary
          containing the relevant data for determining the progress of the step.
          Of the form {'atoms': atoms, 'cell': cell, 'pars': parameters,
                       'status': status, 'result': result, 'id': bead id,
                       'start': starting time}.
    """

    def __init__(
        self,
        latency=1.0,
        offset=0.0,
        name="",
        H=None,
        xref=None,
        vref=0.0,
        pars=None,
        dopbc=False,
        threaded=False,
    ):
        """Initialises FFDebye.

        Args:
           pars: Optional dictionary, giving the parameters needed by the driver.
        """

        # a socket to the communication library is created or linked
        # NEVER DO PBC -- forces here are computed without.
        super(FFDebye, self).__init__(latency, offset, name, pars, dopbc=False)

        if H is None:
            raise ValueError("Must provide the Hessian for the Debye crystal.")
        if xref is None:
            raise ValueError(
                "Must provide a reference configuration for the Debye crystal."
            )

        self.H = H
        self.xref = xref
        self.vref = vref

        eigsys = np.linalg.eigh(self.H)
        info(
            " @ForceField: Hamiltonian eigenvalues: " + " ".join(map(str, eigsys[0])),
            verbosity.medium,
        )

    def evaluate(self, r):
        """A simple evaluator for a harmonic Debye crystal potential."""

        q = r["pos"]
        n3 = len(q)
        if self.H.shape != (n3, n3):
            raise ValueError("Hessian size mismatch")
        if self.xref.shape != (n3,):
            raise ValueError("Reference structure size mismatch")

        d = q - self.xref
        mf = np.dot(self.H, d)

        r["result"] = [
            self.vref + 0.5 * np.dot(d, mf),
            -mf,
            np.zeros((3, 3), float),
            {"raw": ""},
        ]
        r["status"] = "Done"
        r["t_finished"] = time.time()


class FFPlumed(FFEval):
    """Direct PLUMED interface

    Computes forces from a PLUMED input.

    Attributes:
        parameters: A dictionary of the parameters used by the driver. Of the
            form {'name': value}.
        requests: During the force calculation step this holds a dictionary
            containing the relevant data for determining the progress of the step.
            Of the form {'atoms': atoms, 'cell': cell, 'pars': parameters,
                      'status': status, 'result': result, 'id': bead id,
                      'start': starting time}.
    """

    def __init__(
        self,
        latency=1.0e-3,
        offset=0.0,
        name="",
        pars=None,
        dopbc=False,
        threaded=False,
        init_file="",
        plumeddat="",
        plumedstep=0,
        plumed_extras=[],
    ):
        """Initialises FFPlumed.

        Args:
           pars: Optional dictionary, giving the parameters needed by the driver.
        """

        # a socket to the communication library is created or linked
        if plumed is None:
            raise ImportError(
                "Cannot find plumed libraries to link to a FFPlumed object/"
            )
        super(FFPlumed, self).__init__(
            latency, offset, name, pars, dopbc=False, threaded=threaded
        )
        self.plumed = plumed.Plumed()
        self.plumeddat = plumeddat
        self.plumedstep = plumedstep
        self.plumed_extras = plumed_extras
        self.init_file = init_file

        if self.init_file.mode == "xyz":
            infile = open(self.init_file.value, "r")
            myframe = read_file(self.init_file.mode, infile)
            myatoms = myframe["atoms"]
            mycell = myframe["cell"]
            myatoms.q *= unit_to_internal("length", self.init_file.units, 1.0)
            mycell.h *= unit_to_internal("length", self.init_file.units, 1.0)

        self.natoms = myatoms.natoms
        self.plumed.cmd("setRealPrecision", 8)  # i-PI uses double precision
        self.plumed.cmd("setMDEngine", "i-pi")
        self.plumed.cmd("setPlumedDat", self.plumeddat)
        self.plumed.cmd("setNatoms", self.natoms)
        timeunit = 2.4188843e-05  # atomic time to ps
        self.plumed.cmd("setMDTimeUnits", timeunit)
        # given we don't necessarily call plumed once per step, so time does not make
        # sense, we set the time step so that time in plumed is a counter of the number of times
        # called
        self.plumed.cmd("setTimestep", 1 / timeunit)
        self.plumed.cmd(
            "setMDEnergyUnits", 2625.4996
        )  # Pass a pointer to the conversion factor between the energy unit used in your code and kJ mol-1
        self.plumed.cmd(
            "setMDLengthUnits", 0.052917721
        )  # Pass a pointer to the conversion factor between the length unit used in your code and nm
        self.plumedrestart = False
        if self.plumedstep > 0:
            # we are restarting, signal that PLUMED should continue
            self.plumedrestart = True
            self.plumed.cmd("setRestart", 1)
        self.plumed.cmd("init")

        self.plumed_data = {}
        for x in plumed_extras:
            rank = np.zeros(1, dtype=np.int_)
            self.plumed.cmd(f"getDataRank {x}", rank)
            if rank[0] > 1:
                raise ValueError("Cannot retrieve varibles with rank > 1")
            shape = np.zeros(rank[0], dtype=np.int_)
            if shape[0] > 1:
                raise ValueError("Cannot retrieve varibles with size > 1")
            self.plumed.cmd(f"getDataShape {x}", shape)
            self.plumed_data[x] = np.zeros(shape, dtype=np.double)
            self.plumed.cmd(f"setMemoryForData {x}", self.plumed_data[x])

        self.charges = dstrip(myatoms.q) * 0.0
        self.masses = dstrip(myatoms.m)
        self.lastq = np.zeros(3 * self.natoms)
        self.system_force = None  # reference to physical force calculator

    def evaluate(self, r):
        """A wrapper function to call the PLUMED evaluation routines
        and return forces."""

        if self.natoms != len(r["pos"]) / 3:
            raise ValueError(
                "Size of atom array changed after initialization of FFPlumed"
            )

        v = 0.0
        f = np.zeros((self.natoms, 3))
        vir = np.zeros((3, 3))

        self.lastq[:] = r["pos"]
        # for the moment these are set to dummy values taken from an init file.
        # linking with the current value in simulations is non-trivial, as masses
        # are not expected to be the force evaluator's business, and charges are not
        # i-PI's business.
        self.plumed.cmd("setStep", self.plumedstep)
        self.plumed.cmd("setCharges", self.charges)
        self.plumed.cmd("setMasses", self.masses)

        # these instead are set properly. units conversion is done on the PLUMED side
        self.plumed.cmd("setBox", r["cell"][0].T.copy())
        pos = r["pos"].reshape(-1, 3)

        if self.system_force is not None:
            f[:] = dstrip(self.system_force.f).reshape((-1, 3))
            vir[:] = -dstrip(self.system_force.vir)
            self.plumed.cmd("setEnergy", dstrip(self.system_force.pot))

        self.plumed.cmd("setPositions", pos)
        self.plumed.cmd("setForces", f)
        self.plumed.cmd("setVirial", vir)
        self.plumed.cmd("prepareCalc")
        self.plumed.cmd("performCalcNoUpdate")

        bias = np.zeros(1, float)
        self.plumed.cmd("getBias", bias)
        v = bias[0]
        f = f.flatten()
        vir *= -1

        if self.system_force is not None:
            # plumed increments the value of the force, here we need only the correction term
            f[:] -= dstrip(self.system_force.f).flatten()
            vir[:] -= -dstrip(self.system_force.vir)

        extras = {"raw": ""}
        for x in self.plumed_data:
            extras[str(x)] = self.plumed_data[x].copy()

        # nb: the virial is a symmetric tensor, so we don't need to transpose
        r["result"] = [v, f, vir, extras]
        r["status"] = "Done"

    def mtd_update(self, pos, cell):
        """Makes updates to the potential that only need to be triggered
        upon completion of a time step."""

        self.plumedstep += 1
        f = np.zeros((self.natoms, 3))
        vir = np.zeros((3, 3))

        self.plumed.cmd("setStep", self.plumedstep)
        self.plumed.cmd("setCharges", self.charges)
        self.plumed.cmd("setMasses", self.masses)
        rpos = pos.reshape((-1, 3))
        self.plumed.cmd("setPositions", rpos)
        self.plumed.cmd("setBox", cell.T.copy())
        if self.system_force is not None:
            f[:] = dstrip(self.system_force.f).reshape((-1, 3))
            vir[:] = -dstrip(self.system_force.vir)
            self.plumed.cmd("setEnergy", dstrip(self.system_force.pot))
        self.plumed.cmd("setForces", f)
        self.plumed.cmd("setVirial", vir)
        self.plumed.cmd("prepareCalc")
        self.plumed.cmd("performCalcNoUpdate")
        self.plumed.cmd("update")

        return True


class FFYaff(FFEval):
    """Use Yaff as a library to construct a force field"""

    def __init__(
        self,
        latency=1.0,
        offset=0.0,
        name="",
        threaded=False,
        yaffpara=None,
        yaffsys=None,
        yafflog="yaff.log",
        rcut=18.89726133921252,
        alpha_scale=3.5,
        gcut_scale=1.1,
        skin=0,
        smooth_ei=False,
        reci_ei="ewald",
        pars=None,
        dopbc=False,
    ):
        """Initialises FFYaff and enables a basic Yaff force field.

        Args:

           yaffpara: File name of the Yaff parameter file

           yaffsys: File name of the Yaff system file

           yafflog: File name to which Yaff will write some information about the system and the force field

           pars: Optional dictionary, giving the parameters needed by the driver.

           **kwargs: All keyword arguments that can be provided when generating
                     a Yaff force field; see constructor of FFArgs in Yaff code

        """

        from yaff import System, ForceField, log
        import codecs
        import locale
        import atexit

        # a socket to the communication library is created or linked
        super(FFYaff, self).__init__(
            latency, offset, name, pars, dopbc, threaded=threaded
        )

        # A bit weird to use keyword argument for a required argument, but this
        # is also done in the code above.
        if yaffpara is None:
            raise ValueError("Must provide a Yaff parameter file.")

        if yaffsys is None:
            raise ValueError("Must provide a Yaff system file.")

        self.yaffpara = yaffpara
        self.yaffsys = yaffsys
        self.rcut = rcut
        self.alpha_scale = alpha_scale
        self.gcut_scale = gcut_scale
        self.skin = skin
        self.smooth_ei = smooth_ei
        self.reci_ei = reci_ei
        self.yafflog = yafflog

        # Open log file
        logf = open(yafflog, "w")
        # Tell Python to close the file when the script exits
        atexit.register(logf.close)

        # Redirect Yaff log to file
        log._file = codecs.getwriter(locale.getpreferredencoding())(logf)

        self.system = System.from_file(self.yaffsys)
        self.ff = ForceField.generate(
            self.system,
            self.yaffpara,
            rcut=self.rcut,
            alpha_scale=self.alpha_scale,
            gcut_scale=self.gcut_scale,
            skin=self.skin,
            smooth_ei=self.smooth_ei,
            reci_ei=self.reci_ei,
        )

        log._active = False

    def evaluate(self, r):
        """Evaluate the energy and forces with the Yaff force field."""

        q = r["pos"]
        nat = len(q) / 3
        rvecs = r["cell"][0]

        self.ff.update_rvecs(np.ascontiguousarray(rvecs.T, dtype=np.float64))
        self.ff.update_pos(q.reshape((nat, 3)))
        gpos = np.zeros((nat, 3))
        vtens = np.zeros((3, 3))
        e = self.ff.compute(gpos, vtens)

        r["result"] = [e, -gpos.ravel(), -vtens, {"raw": ""}]
        r["status"] = "Done"


class FFsGDML(FFEval):
    """A symmetric Gradient Domain Machine Learning (sGDML) force field.
    Chmiela et al. Sci. Adv., 3(5), e1603015, 2017; Nat. Commun., 9(1), 3887, 2018.
    http://sgdml.org/doc/
    https://github.com/stefanch/sGDML
    """

    def __init__(
        self,
        latency=1.0,
        offset=0.0,
        name="",
        threaded=False,
        sGDML_model=None,
        pars=None,
        dopbc=False,
    ):
        """Initialises FFsGDML

        Args:

           sGDML_model: Filename contaning the sGDML model

        """

        # a socket to the communication library is created or linked
        super(FFsGDML, self).__init__(
            latency, offset, name, pars, dopbc, threaded=threaded
        )

        from ipi.utils.units import unit_to_user

        # --- Load sGDML package ---
        try:
            from sgdml.predict import GDMLPredict
            from sgdml import __version__

            info(" @ForceField: Using sGDML version " + __version__, verbosity.low)
        except ImportError:
            raise ValueError(
                "ERROR: sGDML package not located. Install it via: pip install sgdml"
            )

        # A bit weird to use keyword argument for a required argument, but this
        # is also done in the code above.
        if sGDML_model is None:
            raise ValueError("Must provide a sGDML model file.")

        if dopbc is True:
            raise ValueError("Must set PBCs to False.")

        self.sGDML_model = sGDML_model

        # --- Load sGDML model file. ---
        try:
            self.model = dict(np.load(self.sGDML_model, allow_pickle=True))
            info(
                " @ForceField: sGDML model " + self.sGDML_model + " loaded",
                verbosity.medium,
            )
        except ValueError:
            raise ValueError(
                "ERROR: Reading sGDML model " + self.sGDML_model + " file failed."
            )

        info(
            " @ForceField: IMPORTANT: It is always assumed that the units in"
            + " the provided model file are in Angstroms and kcal/mol.",
            verbosity.low,
        )

        # --- Units ---
        transl_units_names = {
            "Ang": "angstrom",
            "bohr": "atomic_unit",
            "millihartree": "milliatomic_unit",
            "eV": "electronvolt",
            "kcal/mol": "kilocal/mol",
        }
        if "r_unit" in self.model and "e_unit" in self.model:
            distanceUnits = transl_units_names["{}".format(self.model["r_unit"])]
            energyUnits = transl_units_names["{}".format(self.model["e_unit"])]
            info(
                " @ForceField: The units used in the sGDML model are: "
                "distance --> {}, energy --> {}".format(distanceUnits, energyUnits),
                verbosity.low,
            )
        else:
            info(
                " @ForceField: IMPORTANT: Since the sGDML model doesn't have units for distance and energy it is "
                "assumed that the units are in Angstroms and kilocal/mol, respectively.",
                verbosity.low,
            )
            distanceUnits = "angstrom"
            energyUnits = "kilocal/mol"

        # --- Constants ---
        self.bohr_to_distanceUnits = unit_to_user("length", distanceUnits, 1)
        self.energyUnits_to_hartree = unit_to_internal("energy", energyUnits, 1)
        self.forceUnits_to_hartreebohr = (
            self.bohr_to_distanceUnits * self.energyUnits_to_hartree
        )

        # --- Creates predictor ---
        self.predictor = GDMLPredict(self.model)

        info(
            " @ForceField: Optimizing parallelization settings for sGDML FF.",
            verbosity.medium,
        )
        self.predictor.prepare_parallel(n_bulk=1)

    def evaluate(self, r):
        """Evaluate the energy and forces."""

        E, F = self.predictor.predict(r["pos"] * self.bohr_to_distanceUnits)

        r["result"] = [
            E[0] * self.energyUnits_to_hartree,
            F.flatten() * self.forceUnits_to_hartreebohr,
            np.zeros((3, 3), float),
            {"raw": ""},
        ]
        r["status"] = "Done"
        r["t_finished"] = time.time()


class FFCommittee(ForceField):
    """Combines multiple forcefields into a single forcefield object that consolidates
    individual components. Provides the infrastructure to run a simulation based on a
    committee of potentials, and implements the weighted baseline method."""

    def __init__(
        self,
        latency=1.0,
        offset=0.0,
        name="",
        pars=None,
        dopbc=True,
        active=np.array([-1]),
        threaded=True,
        fflist=[],
        ffweights=[],
        alpha=1.0,
        baseline_name="",
        baseline_uncertainty=-1.0,
        active_thresh=0.0,
        active_out=None,
        parse_json=False,
    ):
        # force threaded mode as otherwise it cannot have threaded children
        super(FFCommittee, self).__init__(
            latency=latency,
            offset=offset,
            name=name,
            pars=pars,
            dopbc=dopbc,
            active=active,
            threaded=True,
        )
        if len(fflist) == 0:
            raise ValueError(
                "Committee forcefield cannot be initialized from an empty list"
            )
        self.fflist = fflist
        self.ff_requests = {}
        self.baseline_uncertainty = baseline_uncertainty
        self.baseline_name = baseline_name
        if len(ffweights) == 0 and self.baseline_uncertainty < 0:
            ffweights = np.ones(len(fflist))
        elif len(ffweights) == 0 and self.baseline_uncertainty > 0:
            ffweights = np.ones(len(fflist) - 1)
        if len(ffweights) != len(fflist) and self.baseline_uncertainty < 0:
            raise ValueError("List of weights does not match length of committee model")
        elif len(ffweights) != len(fflist) - 1 and self.baseline_uncertainty > 0:
            raise ValueError("List of weights does not match length of committee model")
        if (self.baseline_name == "") != (self.baseline_uncertainty < 0):
            raise ValueError(
                "Name and the uncertainty of the baseline are not simultaneously defined"
            )
        self.ffweights = ffweights
        self.alpha = alpha
        self.active_thresh = active_thresh
        self.active_out = active_out
        self.parse_json = parse_json

    def bind(self, output_maker):
        super(FFCommittee, self).bind(output_maker)
        if self.active_thresh > 0:
            if self.active_out is None:
                raise ValueError(
                    "Must specify an output file if you want to save structures for active learning"
                )
            else:
                self.active_file = self.output_maker.get_output(self.active_out, "w")

    def start(self):
        for ff in self.fflist:
            ff.start()
        super(FFCommittee, self).start()

    def queue(self, atoms, cell, reqid=-1):
        # launches requests for all of the committee FF objects
        ffh = []
        for ff in self.fflist:
            ffh.append(ff.queue(atoms, cell, reqid))

        # creates the request with the help of the base class,
        # making sure it already contains a handle to the list of FF
        # requests
        req = super(FFCommittee, self).queue(
            atoms, cell, reqid, template=dict(ff_handles=ffh)
        )
        req["t_dispatched"] = time.time()
        return req

    def check_finish(self, r):
        """Checks if all sub-requests associated with a given
        request are finished"""
        for ff_r in r["ff_handles"]:
            if ff_r["status"] != "Done":
                return False
        return True

    def gather(self, r):
        """Collects results from all sub-requests, and assemble the committee of models."""

        r["result"] = [
            0.0,
            np.zeros(len(r["pos"]), float),
            np.zeros((3, 3), float),
            "",
        ]

        # list of pointers to the forcefield requests. shallow copy so we can remove stuff
        com_handles = r["ff_handles"].copy()
        if self.baseline_name != "":
            # looks for the baseline potential, store its value and drops it from the list
            names = [ff.name for ff in self.fflist]

            for i, ff_r in enumerate(com_handles):
                if names[i] == self.baseline_name:
                    baseline_pot = ff_r["result"][0]
                    baseline_frc = ff_r["result"][1]
                    baseline_vir = ff_r["result"][2]
                    baseline_xtr = ff_r["result"][3]
                    com_handles.pop(i)
                    break

        # Gathers the forcefield energetics and extras
        pots = []
        frcs = []
        virs = []
        xtrs = []

        all_have_frc = True
        all_have_vir = True

        for ff_r in com_handles:
            # if required, tries to extract multiple committe members from the extras JSON string
            if "committee_pot" in ff_r["result"][3] and self.parse_json:
                pots += ff_r["result"][3]["committee_pot"]
                if "committee_force" in ff_r["result"][3]:
                    frcs += ff_r["result"][3]["committee_force"]
                    ff_r["result"][3].pop("committee_force")
                else:
                    # if the commitee doesn't have forces, just add the mean force from this model
                    frcs.append(ff_r["result"][1])
                    warning("JSON committee doesn't have forces", verbosity.medium)
                    all_have_frc = False

                if "committee_virial" in ff_r["result"][3]:
                    virs += ff_r["result"][3]["committee_virial"]
                    ff_r["result"][3].pop("committee_virial")
                else:
                    # if the commitee doesn't have virials, just add the mean virial from this model
                    virs.append(ff_r["result"][2])
                    warning("JSON committee doesn't have virials", verbosity.medium)
                    all_have_vir = False

            else:
                pots.append(ff_r["result"][0])
                frcs.append(ff_r["result"][1])
                virs.append(ff_r["result"][2])

        pots = np.array(pots)
        if len(pots) != len(frcs) and len(frcs) > 1:
            raise ValueError(
                "If the committee returns forces, we need *all* components"
            )
        frcs = np.array(frcs).reshape(len(frcs), -1)

        if len(pots) != len(virs) and len(virs) > 1:
            raise ValueError(
                "If the committee returns virials, we need *all* components"
            )
        virs = np.array(virs).reshape(-1, 3, 3)

        xtrs.append(ff_r["result"][3])

        # Computes the mean energetics
        mean_pot = np.mean(pots, axis=0)
        mean_frc = np.mean(frcs, axis=0)
        mean_vir = np.mean(virs, axis=0)

        # Rescales the committee energetics so that their standard deviation corresponds to the error
        rescaled_pots = np.asarray(
            [mean_pot + self.alpha * (pot - mean_pot) for pot in pots]
        )
        rescaled_frcs = np.asarray(
            [mean_frc + self.alpha * (frc - mean_frc) for frc in frcs]
        )
        rescaled_virs = np.asarray(
            [mean_vir + self.alpha * (vir - mean_vir) for vir in virs]
        )

        # Calculates the error associated with the committee
        var_pot = np.var(rescaled_pots, ddof=1)
        std_pot = np.sqrt(var_pot)

        if self.baseline_name != "":
            if not (all_have_frc and all_have_vir):
                raise ValueError(
                    "Cannot use weighted baseline without a force ensemble"
                )

            # Computes the additional component of the energetics due to a position
            # dependent weight. This is based on the assumption that V_committee is
            # a correction over the baseline, that V = V_baseline + V_committe, that
            # V_baseline has an uncertainty given by baseline_uncertainty,
            # and V_committee the committee error. Then
            # V = V_baseline + s_b^2/(s_c^2+s_b^2) V_committe

            s_b2 = self.baseline_uncertainty**2

            nmodels = len(pots)
            uncertain_frc = (
                self.alpha**2
                * np.sum(
                    [
                        (pot - mean_pot) * (frc - mean_frc)
                        for pot, frc in zip(pots, frcs)
                    ],
                    axis=0,
                )
                / (nmodels - 1)
            )
            uncertain_vir = (
                self.alpha**2
                * np.sum(
                    [
                        (pot - mean_pot) * (vir - mean_vir)
                        for pot, vir in zip(pots, virs)
                    ],
                    axis=0,
                )
                / (nmodels - 1)
            )

            # Computes the final average energetics
            final_pot = baseline_pot + mean_pot * s_b2 / (s_b2 + var_pot)
            final_frc = (
                baseline_frc
                + mean_frc * s_b2 / (s_b2 + var_pot)
                - 2.0 * mean_pot * s_b2 / (s_b2 + var_pot) ** 2 * uncertain_frc
            )
            final_vir = (
                baseline_vir
                + mean_vir * s_b2 / (s_b2 + var_pot)
                - 2.0 * mean_pot * s_b2 / (s_b2 + var_pot) ** 2 * uncertain_vir
            )

            # Sets the output of the committee model.
            r["result"][0] = final_pot
            r["result"][1] = final_frc
            r["result"][2] = final_vir
        else:
            # Sets the output of the committee model.
            r["result"][0] = mean_pot
            r["result"][1] = mean_frc
            r["result"][2] = mean_vir

        r["result"][3] = {
            "committee_pot": rescaled_pots,
            "committee_uncertainty": std_pot,
        }

        if all_have_frc:
            r["result"][3]["committee_force"] = rescaled_frcs.reshape(
                len(rescaled_pots), -1
            )
        if all_have_vir:
            r["result"][3]["committee_virial"] = rescaled_virs.reshape(
                len(rescaled_pots), -1
            )

        if self.baseline_name != "":
            r["result"][3]["baseline_pot"] = (baseline_pot,)
            r["result"][3]["baseline_force"] = (baseline_frc,)
            r["result"][3]["baseline_virial"] = ((baseline_vir.flatten()),)
            r["result"][3]["baseline_extras"] = (baseline_xtr,)
            r["result"][3]["wb_mixing"] = (s_b2 / (s_b2 + var_pot),)

        # "dissolve" the extras dictionaries into a list
        for k in xtrs[0].keys():
            if ("committee_" + k) in r["result"][3].keys():
                raise ValueError(
                    "Name clash between extras key "
                    + k
                    + " and default committee extras"
                )
            r["result"][3][("committee_" + k)] = []
            for x in xtrs:
                r["result"][3][("committee_" + k)].append(x[k])

        if self.active_thresh > 0.0 and std_pot > self.active_thresh:
            dumps = json.dumps(
                {
                    "position": list(r["pos"]),
                    "cell": list(r["cell"][0].flatten()),
                    "uncertainty": std_pot,
                }
            )
            self.active_file.write(dumps)

        # releases the requests from the committee FF
        for ff, ff_r in zip(self.fflist, r["ff_handles"]):
            ff.release(ff_r)

    def poll(self):
        """Polls the forcefield object to check if it has finished."""

        with self._threadlock:
            for r in self.requests:
                if r["status"] != "Done" and self.check_finish(r):
                    r["t_finished"] = time.time()
                    self.gather(r)
                    r["result"][0] -= self.offset
                    r["status"] = "Done"


class PhotonDriver:
    """
    Photon driver for a single cavity mode
    """

    def __init__(self, apply_photon=True, E0=1e-4, omega_c=0.01, ph_rep="loose"):
        """
        Initialise PhotonDriver

        In this implementation, the photonic masses are set as 1 a.u.

        Args:
            apply_photon: Determine if applying light-matter interactions
            E0: varepsilon in the paper (doi.org/10.1073/pnas.2009272117), light-matter coupling strength
            omega_c: cavity frequency at normal incidence
            ph_rep: 'loose' or 'dense'. In the current implementation, two energy-degenerate photon modes polarized along x and y directions
                are coupled to the molecular system. If 'loose', the cavity photons polarized along the x, y directions are represented by two 'L' atoms;
                the x dimension of the first 'L' atom is coupled to the molecules, and the y dimension of the second 'L' atom is coupled to the molecules.
                If 'dense', the cavity photons polarized along the x, y directions are represented by one 'L' atom;
                the x and y dimensions of this 'L' atom are coupled to the molecules.
        """
        self.apply_photon = apply_photon
        self.E0 = E0
        self.omega_c = omega_c

        if self.apply_photon == False:
            self.n_mode = 0
        elif self.apply_photon == True:
            self.n_mode = 1
            self.n_grid = 1
            self.ph_rep = ph_rep
            self.init_photon()

    def init_photon(self):
        """
        Initialize the photon environment parameters
        """

        if self.ph_rep == "loose":
            self.n_photon = 2 * self.n_mode
        elif self.ph_rep == "dense":
            self.n_photon = self.n_mode
        self.n_photon_3 = self.n_photon * 3
        self.pos_ph = np.zeros(self.n_photon_3)

        # construct cavity mode frequency array for all photons
        self.omega_k = np.array([self.omega_c])
        if self.ph_rep == "loose":
            self.omega_klambda = np.concatenate((self.omega_k, self.omega_k))
        elif self.ph_rep == "dense":
            self.omega_klambda = self.omega_k
        self.omega_klambda3 = np.reshape(
            np.array([[x, x, x] for x in self.omega_klambda]), -1
        )

        # construct varepsilon array for all photons
        self.varepsilon_k = self.E0
        self.varepsilon_klambda = (
            self.E0 * self.omega_klambda / np.min(self.omega_klambda)
        )
        self.varepsilon_klambda3 = (
            self.E0 * self.omega_klambda3 / np.min(self.omega_klambda3)
        )

        # cavity mode function acting on the dipole moment
        self.ftilde_kx = np.array([1.0])
        self.ftilde_ky = np.array([1.0])
        self.ftilde_kx3 = np.reshape(np.array([[x, x, x] for x in self.ftilde_kx]), -1)
        self.ftilde_ky3 = np.reshape(np.array([[x, x, x] for x in self.ftilde_ky]), -1)

    def split_atom_ph_coord(self, pos):
        """
        Split atomic and photonic coordinates and update our photonic coordinates

        Args:
            pos: A 3*N position numpy array, [1x, 1y, 1z, 2x, ...]

        Returns:
            Atomic coordinates, Photonic coordinates
        """
        if self.apply_photon:
            pos_at = pos[: -self.n_photon_3]
            pos_ph = pos[-self.n_photon_3 :]
            self.pos_ph = pos_ph
        else:
            pos_at = pos
            pos_ph = pos[0:0]
        return pos_at, pos_ph

    def get_ph_energy(self, dx_array, dy_array):
        """
        Calculate the total photonic potential energy, including the light-matter
        interaction and dipole self energy

        Args:
            dx_array: x-direction dipole array of molecular subsystems
            dy_array: y-direction dipole array of molecular subsystems

        Returns:
            total energy of photonic system
        """
        # calculate the photonic potential energy
        e_ph = np.sum(0.5 * self.omega_klambda3**2 * self.pos_ph**2)

        # calculate the dot products between mode functions and dipole array
        d_dot_f_x = np.dot(self.ftilde_kx, dx_array)
        d_dot_f_y = np.dot(self.ftilde_ky, dy_array)

        # calculate the light-matter interaction
        if self.ph_rep == "loose":
            e_int_x = np.sum(
                self.varepsilon_k * d_dot_f_x * self.pos_ph[: self.n_mode * 3 : 3]
            )
            e_int_y = np.sum(
                self.varepsilon_k * d_dot_f_y * self.pos_ph[1 + self.n_mode * 3 :: 3]
            )
        elif self.ph_rep == "dense":
            e_int_x = np.sum(self.varepsilon_k * d_dot_f_x * self.pos_ph[::3])
            e_int_y = np.sum(self.varepsilon_k * d_dot_f_y * self.pos_ph[1::3])

        # calculate the dipole self-energy term
        dse = np.sum(
            (self.varepsilon_k**2 / 2.0 / self.omega_k**2)
            * (d_dot_f_x**2 + d_dot_f_y**2)
        )

        e_tot = e_ph + e_int_x + e_int_y + dse

        return e_tot

    def get_ph_forces(self, dx_array, dy_array):
        """
        Calculate the photonic forces

        Args:
            dx_array: x-direction dipole array of molecular subsystems
            dy_array: y-direction dipole array of molecular subsystems

        Returns:
            force array of all photonic dimensions (3*nphoton) [1x, 1y, 1z, 2x..]
        """
        # calculat the bare photonic contribution of the force
        f_ph = -self.omega_klambda3**2 * self.pos_ph

        # calculate the dot products between mode functions and dipole array
        d_dot_f_x = np.dot(self.ftilde_kx, dx_array)
        d_dot_f_y = np.dot(self.ftilde_ky, dy_array)

        # calculate the force due to light-matter interactions
        if self.ph_rep == "loose":
            f_ph[: self.n_mode * 3 : 3] -= self.varepsilon_k * d_dot_f_x
            f_ph[self.n_mode * 3 + 1 :: 3] -= self.varepsilon_k * d_dot_f_y
        elif self.ph_rep == "dense":
            f_ph[::3] -= self.varepsilon_k * d_dot_f_x
            f_ph[1::3] -= self.varepsilon_k * d_dot_f_y
        return f_ph

    def get_nuc_cav_forces(self, dx_array, dy_array, charge_array_bath):
        """
        Calculate the photonic forces on nuclei from MM partial charges

        Args:
            dx_array: x-direction dipole array of molecular subsystems
            dy_array: y-direction dipole array of molecular subsystems
            charge_array_bath: partial charges of all atoms in a single bath

        Returns:
            force array of all nuclear dimensions (3*natoms) [1x, 1y, 1z, 2x..]
        """

        # calculate the dot products between mode functions and dipole array
        d_dot_f_x = np.dot(self.ftilde_kx, dx_array)
        d_dot_f_y = np.dot(self.ftilde_ky, dy_array)

        # cavity force on x direction
        if self.ph_rep == "loose":
            Ekx = self.varepsilon_k * self.pos_ph[: self.n_mode * 3 : 3]
            Eky = self.varepsilon_k * self.pos_ph[self.n_mode * 3 + 1 :: 3]
        elif self.ph_rep == "dense":
            Ekx = self.varepsilon_k * self.pos_ph[::3]
            Eky = self.varepsilon_k * self.pos_ph[1::3]
        Ekx += self.varepsilon_k**2 / self.omega_k**2 * d_dot_f_x
        Eky += self.varepsilon_k**2 / self.omega_k**2 * d_dot_f_y

        # dimension of independent baths (xy grid points)
        coeff_x = np.dot(np.transpose(Ekx), self.ftilde_kx)
        coeff_y = np.dot(np.transpose(Eky), self.ftilde_ky)
        fx = -np.kron(coeff_x, charge_array_bath)
        fy = -np.kron(coeff_y, charge_array_bath)
        return fx, fy


class FFCavPhSocket(FFSocket):
    """
    Socket for dealing with cavity photons interacting with molecules by
    Tao E. Li @ 2023-02-25
    Check https://doi.org/10.1073/pnas.2009272117 for details

    Independent bath approximation will be made to communicate with many sockets
    """

    def __init__(
        self,
        latency=1.0,
        offset=0.0,
        name="",
        pars=None,
        dopbc=False,
        active=np.array([-1]),
        threaded=True,
        interface=None,
        charge_array=None,
        apply_photon=True,
        E0=1e-4,
        omega_c=0.01,
        ph_rep="loose",
    ):
        """Initialises FFCavPhFPSocket.

        Args:
           latency: The number of seconds the socket will wait before updating
              the client list.
           name: The name of the forcefield.
           pars: A dictionary used to initialize the forcefield, if required.
              Of the form {'name1': value1, 'name2': value2, ... }.
           dopbc: Decides whether or not to apply the periodic boundary conditions
              before sending the positions to the client code.
           interface: The object used to create the socket used to interact
              with the client codes.
           charge_array: An N-dimensional numpy array for fixed point charges of all atoms
           apply_photon: If add photonic degrees of freedom in the dynamics
           E0: Effective light-matter coupling strength
           omega_c: Cavity mode frequency
           ph_rep: A string to control how to represent the photonic coordinates: 'loose' or 'dense'.
                In the current implementation, two energy-degenerate photon modes polarized along x and
                y directions are coupled to the molecular system. If 'loose', the cavity photons polarized
                along the x, y directions are represented by two 'L' atoms; the x dimension of the first
                'L' atom is coupled to the molecules, and the y dimension of the second 'L' atom is coupled
                to the molecules. If 'dense', the cavity photons polarized along the x, y directions are
                represented by one 'L' atom; the x and y dimensions of this 'L' atom are coupled to the molecules.
        """

        # a socket to the communication library is created or linked
        super(FFCavPhSocket, self).__init__(
            latency, offset, name, pars, dopbc, active, threaded, interface
        )

        # definition of independent baths
        self.n_independent_bath = 1
        self.charge_array = charge_array

        # store photonic variables
        self.apply_photon = apply_photon
        self.E0 = E0
        self.omega_c = omega_c
        self.ph_rep = ph_rep
        # define the photon environment
        self.ph = PhotonDriver(
            apply_photon=apply_photon, E0=E0, omega_c=omega_c, ph_rep=ph_rep
        )

        self._getallcount = 0

    def calc_dipole_xyz_mm(self, pos, n_bath, charge_array_bath):
        """
        Calculate the x, y, and z components of total dipole moment for a single molecular subsystem (bath)

        Args:
            pos: position of all atoms (3*n) in all subsystems
            n_bath: total number of molecular subsystems (baths)
            charge_array_bath: charge_array of all atoms (n) in a single subsystem (bath)

        Returns:
            dx_array, dy_array, dz_array: total dipole moment array along x, y, and z directions
        """
        ndim_tot = np.size(pos)
        ndim_local = int(ndim_tot // n_bath)

        dx_array, dy_array, dz_array = [], [], []
        for idx in range(n_bath):
            pos_bath = pos[ndim_local * idx : ndim_local * (idx + 1)]
            # check the dimension of charge array
            if np.size(pos_bath[::3]) != np.size(charge_array_bath):
                softexit.trigger(
                    "The size of charge array = {}  does not match the size of atoms = {} ".format(
                        np.size(charge_array_bath), np.size(pos_bath[::3])
                    )
                )
            dx = np.sum(pos_bath[::3] * charge_array_bath)
            dy = np.sum(pos_bath[1::3] * charge_array_bath)
            dz = np.sum(pos_bath[2::3] * charge_array_bath)
            dx_array.append(dx)
            dy_array.append(dy)
            dz_array.append(dz)
        dx_array = np.array(dx_array)
        dy_array = np.array(dy_array)
        dz_array = np.array(dz_array)
        return dx_array, dy_array, dz_array

    def queue(self, atoms, cell, reqid=-1):
        """Adds a request.

        Note that the pars dictionary need to be sent as a string of a
        standard format so that the initialisation of the driver can be done.

        Args:
            atoms: An Atoms object giving the atom positions.
            cell: A Cell object giving the system box.
            pars: An optional dictionary giving the parameters to be sent to the
                driver for initialisation. Defaults to {}.
            reqid: An optional integer that identifies requests of the same type,
               e.g. the bead index

        Returns:
            A list giving the status of the request of the form {'pos': An array
            giving the atom positions folded back into the unit cell,
            'cell': Cell object giving the system box, 'pars': parameter string,
            'result': holds the result as a list once the computation is done,
            'status': a string labelling the status of the calculation,
            'id': the id of the request, usually the bead number, 'start':
            the starting time for the calculation, used to check for timeouts.}.
        """

        par_str = " "

        if not self.pars is None:
            for k, v in list(self.pars.items()):
                par_str += k + " : " + str(v) + " , "
        else:
            par_str = " "

        pbcpos = dstrip(atoms.q).copy()

        # Indexes come from input in a per atom basis and we need to make a per atom-coordinate basis
        # Reformat indexes for full system (default) or piece of system
        # active atoms do not change but we only know how to build this array once we get the positions once
        if self.iactive is None:
            if self.active[0] == -1:
                activehere = np.arange(len(pbcpos))
            else:
                activehere = np.array(
                    [[3 * n, 3 * n + 1, 3 * n + 2] for n in self.active]
                )

            # Reassign active indexes in order to use them
            activehere = activehere.flatten()

            # Perform sanity check for active atoms
            if len(activehere) > len(pbcpos) or activehere[-1] > (len(pbcpos) - 1):
                raise ValueError("There are more active atoms than atoms!")

            self.iactive = activehere

        newreq_lst = []

        # 1. split coordinates to atoms and photons
        pbcpos_atoms, pbcpos_phs = self.ph.split_atom_ph_coord(pbcpos)
        ndim_tot = np.size(pbcpos_atoms)
        ndim_local = int(ndim_tot // self.n_independent_bath)

        # 2. for atomic coordinates, we now evaluate their atomic forces
        for idx in range(self.n_independent_bath):
            pbcpos_local = pbcpos_atoms[
                ndim_local * idx : ndim_local * (idx + 1)
            ].copy()
            iactive_local = self.iactive[0:ndim_local]
            # Let's try to do PBC for the small regions
            if self.dopbc:
                cell.array_pbc(pbcpos_local)
            newreq_local = ForceRequest(
                {
                    "id": int(reqid * self.n_independent_bath) + idx,
                    "pos": pbcpos_local,
                    "active": iactive_local,
                    "cell": (dstrip(cell.h).copy(), dstrip(cell.ih).copy()),
                    "pars": par_str,
                    "result": None,
                    "status": "Queued",
                    "start": -1,
                    "t_queued": time.time(),
                    "t_dispatched": 0,
                    "t_finished": 0,
                }
            )
            newreq_lst.append(newreq_local)

        with self._threadlock:
            for newreq in newreq_lst:
                self.requests.append(newreq)
                self._getallcount += 1

        if not self.threaded:
            self.poll()

        # sleeps until all the new requests have been evaluated
        for self.request in newreq_lst:
            while self.request["status"] != "Done":
                if self.request["status"] == "Exit" or softexit.triggered:
                    # now, this is tricky. we are stuck here and we cannot return meaningful results.
                    # if we return, we may as well output wrong numbers, or mess up things.
                    # so we can only call soft-exit and wait until that is done. then kill the thread
                    # we are in.
                    softexit.trigger(" @ FORCES : cannot return so will die off here")
                    while softexit.exiting:
                        time.sleep(self.latency)
                    sys.exit()
                time.sleep(self.latency)

            """
            with self._threadlock:
                self._getallcount -= 1

            # releases just once, but wait for all requests to be complete
            if self._getallcount == 0:
                self.release(self.request)
                self.request = None
            else:
                while self._getallcount > 0:
                    time.sleep(self.latency)
            """
            self.release(self.request)
            self.request = None

        # ...atomic forces have been calculated at this point

        # 3. At this moment, we combine the small requests to a big mega request (updated results)
        result_tot = [0.0, np.zeros(len(pbcpos), float), np.zeros((3, 3), float), {}]
        for idx, newreq in enumerate(newreq_lst):
            u, f, vir, extra = newreq["result"]
            result_tot[0] += u
            result_tot[1][ndim_local * idx : ndim_local * (idx + 1)] = f
            result_tot[2] += vir
            result_tot[3][idx] = extra

        if self.ph.apply_photon:
            # 4. calculate total dipole moment array for N baths
            dx_array, dy_array, dz_array = self.calc_dipole_xyz_mm(
                pos=pbcpos_atoms,
                n_bath=self.n_independent_bath,
                charge_array_bath=self.charge_array,
            )
            # check the size of photon modes + molecules to match the total number of particles
            if (
                self.ph.n_photon + self.n_independent_bath * self.charge_array.size
                != int(len(pbcpos) // 3)
            ):
                softexit.trigger(
                    "Total number of photons + molecules does not match total number of particles"
                )
            # info("mux = %.6f muy = %.6f muz = %.6f [units of a.u.]" %(dipole_x_tot, dipole_y_tot, dipole_z_tot), verbosity.medium)
            # 5. calculate photonic contribution of total energy
            e_ph = self.ph.get_ph_energy(dx_array=dx_array, dy_array=dy_array)
            # 6. calculate photonic forces
            f_ph = self.ph.get_ph_forces(dx_array=dx_array, dy_array=dy_array)
            # 7. calculate cavity forces on nuclei
            fx_cav, fy_cav = self.ph.get_nuc_cav_forces(
                dx_array=dx_array,
                dy_array=dy_array,
                charge_array_bath=self.charge_array,
            )
            # 8. add cavity effects to our output
            result_tot[0] += e_ph
            result_tot[1][:ndim_tot:3] += fx_cav
            result_tot[1][1:ndim_tot:3] += fy_cav
            result_tot[1][ndim_tot:] = f_ph

        result_tot[0] -= self.offset

        # At this moment, we have sucessfully gathered the CavMD forces
        newreq = ForceRequest(
            {
                "id": reqid,
                "pos": pbcpos,
                "active": self.iactive,
                "cell": (dstrip(cell.h).copy(), dstrip(cell.ih).copy()),
                "pars": par_str,
                "result": result_tot,
                "status": newreq_lst[-1]["status"],
                "start": newreq_lst[0]["start"],
                "t_queued": newreq_lst[0]["t_queued"],
                "t_dispatched": newreq_lst[0]["t_dispatched"],
                "t_finished": newreq_lst[-1]["t_finished"],
            }
        )

        return newreq

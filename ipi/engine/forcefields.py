"""Contains the classes that connect the driver to the python code.

ForceField objects are force providers, i.e. they are the abstraction
layer for a driver that gets positions and returns forces (and energy).
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import time
import threading

import numpy as np

from ipi.utils.softexit import softexit
from ipi.utils.messages import verbosity
from ipi.utils.messages import info
from ipi.interfaces.sockets import InterfaceSocket
from ipi.utils.depend import dobject
from ipi.utils.depend import dstrip
from ipi.utils.io import read_file
from ipi.utils.units import unit_to_internal


__all__ = ['ForceField', 'FFSocket', 'FFLennardJones', 'FFDebye', 'FFPlumed', 'FFYaff']


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


class ForceField(dobject):

    """Base forcefield class.

    Gives the standard methods and quantities needed in all the forcefield
    classes.

    Attributes:
        pars: A dictionary of the parameters needed to initialize the forcefield.
            Of the form {'name1': value1, 'name2': value2, ... }.
        name: The name of the forcefield.
        latency: A float giving the number of seconds the socket will wait
            before updating the client list.
        requests: A list of all the jobs to be given to the client codes.
        dopbc: A boolean giving whether or not to apply the periodic boundary
            conditions before sending the positions to the client code.
        _thread: The thread on which the socket polling loop is being run.
        _doloop: A list of booleans. Used to decide when to stop running the
            polling loop.
        _threadlock: Python handle used to lock the thread held in _thread.
    """

    def __init__(self, latency=1.0, name="", pars=None, dopbc=True, active=np.array([-1])):
        """Initialises ForceField.

        Args:
            latency: The number of seconds the socket will wait before updating
                the client list.
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
        self.requests = []
        self.dopbc = dopbc
        self.active = active
        self._thread = None
        self._doloop = [False]
        self._threadlock = threading.Lock()

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
            for k, v in self.pars.items():
                par_str += k + " : " + str(v) + " , "
        else:
            par_str = " "

        pbcpos = dstrip(atoms.q).copy()

        # Indexes come from input in a per atom basis and we need to make a per atom-coordinate basis
        # Reformat indexes for full system (default) or piece of system
        if self.active[0] == -1:
            activehere = np.array([i for i in range(len(pbcpos))])
        else:
            activehere = np.array([[3 * n, 3 * n + 1, 3 * n + 2] for n in self.active])

        # Reassign active indexes in order to use them
        activehere = activehere.flatten()

        # Perform sanity check for active atoms
        if (len(activehere) > len(pbcpos) or activehere[-1] > (len(pbcpos) - 1)):
            raise ValueError("There are more active atoms than atoms!")

        if self.dopbc:
            cell.array_pbc(pbcpos)

        newreq = ForceRequest({
            "id": reqid,
            "pos": pbcpos,
            "active": activehere,
            "cell": (dstrip(cell.h).copy(), dstrip(cell.ih).copy()),
            "pars": par_str,
            "result": None,
            "status": "Queued",
            "start": -1,
            "t_queued": time.time(),
            "t_dispatched": 0,
            "t_finished": 0
        })

        with self._threadlock:
            self.requests.append(newreq)

        return newreq

    def poll(self):
        """Polls the forcefield object to check if it has finished."""

        with self._threadlock:
            for r in self.requests:
                if r["status"] == "Queued":
                    r["t_dispatched"] = time.time()
                    r["result"] = [0.0, np.zeros(len(r["pos"]), float), np.zeros((3, 3), float), ""]
                    r["status"] = "Done"
                    r["t_finished"] = time.time()

    def _poll_loop(self):
        """Polling loop.

        Loops over the different requests, checking to see when they have
        finished.
        """

        info(" @ForceField: Starting the polling thread main loop.", verbosity.low)
        while self._doloop[0]:
            if len(self.requests) == 0 :
                time.sleep(self.latency)
            else:
                self.poll()

    def release(self, request):
        """Shuts down the client code interface thread.

        Args:
            request: The id of the job to release.
        """

        """Frees up a request."""

        with self._threadlock:
            if request in self.requests:
                try:
                    self.requests.remove(request)
                except ValueError:
                    print "failed removing request", id(request), ' ',
                    print [id(r) for r in self.requests], "@", threading.currentThread()
                    raise

    def stop(self):
        """Dummy stop method."""

        self._doloop[0] = False
        for r in self.requests:
            r["status"] = "Exit"

    def run(self):
        """Spawns a new thread.

        Splits the main program into two threads, one that runs the polling loop
        which updates the client list, and one which gets the data.

        Raises:
            NameError: Raised if the polling thread already exists.
        """

        if not self._thread is None:
            raise NameError("Polling thread already started")

        self._doloop[0] = True
        self._thread = threading.Thread(target=self._poll_loop, name="poll_" + self.name)
        self._thread.daemon = True
        self._thread.start()
        softexit.register_function(self.softexit)
        softexit.register_thread(self._thread, self._doloop)

    def softexit(self):
        """ Takes care of cleaning up upon softexit """

        self.stop()

    def update(self):
        """ Makes updates to the potential that only need to be triggered
        upon completion of a time step. """

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

    def __init__(self, latency=1.0, name="", pars=None, dopbc=True, active=np.array([-1]), interface=None):
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
        super(FFSocket, self).__init__(latency, name, pars, dopbc, active)
        if interface is None:
            self.socket = InterfaceSocket()
        else:
            self.socket = interface
        self.socket.requests = self.requests

    def poll(self):
        """Function to check the status of the client calculations."""

        self.socket.poll()

    def run(self):
        """Spawns a new thread."""

        self.socket.open()
        super(FFSocket, self).run()

    def stop(self):
        """Closes the socket and the thread."""

        super(FFSocket, self).stop()
        if self._thread is not None:
            # must wait until loop has ended before closing the socket
            self._thread.join()
        self.socket.close()


class FFLennardJones(ForceField):

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

    def __init__(self, latency=1.0e-3, name="", pars=None, dopbc=False):
        """Initialises FFLennardJones.

        Args:
           pars: Optional dictionary, giving the parameters needed by the driver.
        """

        # check input - PBCs are not implemented here
        if dopbc:
            raise ValueError("Periodic boundary conditions are not supported by FFLennardJones.")

        # a socket to the communication library is created or linked
        super(FFLennardJones, self).__init__(latency, name, pars, dopbc=False)
        self.epsfour = float(self.pars["eps"]) * 4
        self.sixepsfour = 6 * self.epsfour
        self.sigma2 = float(self.pars["sigma"]) * float(self.pars["sigma"])

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

            x6 = (self.sigma2 / rij2)**3
            x12 = x6**2

            v += (x12 - x6).sum()
            dij *= (self.sixepsfour * (2.0 * x12 - x6) / rij2)[:, np.newaxis]
            f[i] += dij.sum(axis=0)
            f[:i] -= dij

        v *= self.epsfour

        r["result"] = [v, f.reshape(nat * 3), np.zeros((3, 3), float), ""]
        r["status"] = "Done"


try:
    import quippy
except Exception as e:
    quippy = None
    quippy_exc = e


class FFQUIP(ForceField):

    """Basic fully pythonic force provider.

    Computes an arbitrary interaction potential implemented in QUIP. 
    Parallel evaluation with threads.

    Attributes:
        parameters: A dictionary of the parameters used by QUIP. Of the
            form {'name': value}.
        requests: During the force calculation step this holds a dictionary
            containing the relevant data for determining the progress of the step.
            Of the form {'atoms': atoms, 'cell': cell, 'pars': parameters,
                         'status': status, 'result': result, 'id': bead id,
                         'start': starting time}.
    """

    def __init__(self, init_file, args_str, param_file, latency=1.0e-3, name="", pars=None, dopbc=True):
        """Initialises QUIP.

        Args:
        pars: Mandatory dictionaru, giving the parameters needed by QUIP.
        """
        if quippy is None:
            raise ImportError("QUIPPY import failed due to exception : " + str(e))

        # a socket to the communication library is created or linked
        super(FFQUIP, self).__init__(latency, name, pars, dopbc)
        self.init_file = init_file
        self.args_str = args_str
        self.param_file = param_file

        # Initializes an atoms object and the interaction potential
        self.atoms = quippy.Atoms(self.init_file)
        self.pot = quippy.Potential(self.args_str, param_filename=self.param_file)
       
        # Initializes the conversion factors from i-pi to QUIP
        self.len_conv = 0.529177 
        self.energy_conv = 27.211386

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

    def evaluate(self, r):
        """ The function that evaluates the
        QUIP interaction potential."""

        # Obtains the positions and the cell.
        q = r["pos"].reshape((-1, 3)) 
        h, ih = r["cell"]

        nat = len(q)

        # Performs conversion of units.
        q *= self.len_conv
        h *= self.len_conv

        # Updates the QUIP atoms object.
        self.atoms.set_positions(q)
        self.atoms.set_cell(h)

        # Calculates the energies, forces and the virial.
        self.pot.calc(self.atoms, energy=True, force=True, virial=True) 

        # Obtains the energetics and converts to i-pi units.
        u = self.atoms.energy  / self.energy_conv
        f = self.atoms.force.T.flatten()   / self.energy_conv * self.len_conv 
        v = np.triu(self.atoms.virial)

        r["result"] = [u, f.reshape(nat * 3), v, ""]
        r["status"] = "Done"


class FFDebye(ForceField):

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

    def __init__(self, latency=1.0, name="", H=None, xref=None, vref=0.0, pars=None, dopbc=False, threaded=True):
        """Initialises FFDebye.

        Args:
           pars: Optional dictionary, giving the parameters needed by the driver.
        """

        # a socket to the communication library is created or linked
        # NEVER DO PBC -- forces here are computed without.
        super(FFDebye, self).__init__(latency, name, pars, dopbc=False)

        if H is None:
            raise ValueError("Must provide the Hessian for the Debye crystal.")
        if xref is None:
            raise ValueError("Must provide a reference configuration for the Debye crystal.")

        self.H = H
        self.xref = xref
        self.vref = vref

        eigsys = np.linalg.eigh(self.H)
        info(" @ForceField: Hamiltonian eigenvalues: " + ' '.join(map(str, eigsys[0])), verbosity.medium)

    def poll(self):
        """ Polls the forcefield checking if there are requests that should
        be answered, and if necessary evaluates the associated forces and energy. """

        # we have to be thread-safe, as in multi-system mode this might get called by many threads at once
        with self._threadlock:
            for r in self.requests:
                if r["status"] == "Queued":
                    r["status"] = "Running"
                    self.evaluate(r)

    def evaluate(self, r):
        """ A simple evaluator for a harmonic Debye crystal potential. """

        q = r["pos"]
        n3 = len(q)
        if self.H.shape != (n3, n3):
            raise ValueError("Hessian size mismatch")
        if self.xref.shape != (n3,):
            raise ValueError("Reference structure size mismatch")

        d = q - self.xref
        mf = np.dot(self.H, d)

        r["result"] = [self.vref + 0.5 * np.dot(d, mf), -mf, np.zeros((3, 3), float), ""]
        r["status"] = "Done"
        r["t_finished"] = time.time()


try:
    import plumed
except:
    plumed = None


class FFPlumed(ForceField):
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

    def __init__(self, latency=1.0e-3, name="", pars=None, dopbc=False, init_file="", plumeddat="", precision=8, plumedstep=0):
        """Initialises FFPlumed.

        Args:
           pars: Optional dictionary, giving the parameters needed by the driver.
        """

        # a socket to the communication library is created or linked
        if plumed is None:
            raise ImportError("Cannot find plumed libraries to link to a FFPlumed object/")
        super(FFPlumed, self).__init__(latency, name, pars, dopbc=False)
        self.plumed = plumed.Plumed(precision)
        self.precision = precision
        self.plumeddat = plumeddat
        self.plumedstep = plumedstep
        self.init_file = init_file

        if self.init_file.mode == "xyz":
            infile = open(self.init_file.value, "r")
            myframe = read_file(self.init_file.mode, infile)
            myatoms = myframe['atoms']
            mycell = myframe['cell']
            myatoms.q *= unit_to_internal("length", self.init_file.units, 1.0)
            mycell.h *= unit_to_internal("length", self.init_file.units, 1.0)

        self.natoms = myatoms.natoms
        self.plumed.cmd("setNatoms", self.natoms)
        self.plumed.cmd("setPlumedDat", self.plumeddat)
        self.plumed.cmd("setTimestep", 1.)
        self.plumed.cmd("setMDEnergyUnits", 2625.4996)        # Pass a pointer to the conversion factor between the energy unit used in your code and kJ mol-1
        self.plumed.cmd("setMDLengthUnits", 0.052917721)        # Pass a pointer to the conversion factor between the length unit used in your code and nm
        self.plumed.cmd("setMDTimeUnits", 2.4188843e-05)
        self.plumedrestart = False
        if self.plumedstep > 0:
            # we are restarting, signal that PLUMED should continue
            self.plumedrestart = True
            self.plumed.cmd("setRestart", 1)
        self.plumed.cmd("init")
        self.charges = dstrip(myatoms.q) * 0.0
        self.masses = dstrip(myatoms.m)
        self.lastq = np.zeros(3 * self.natoms)

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
                    r["t_finished"] = time.time()

    def evaluate(self, r):
        """A wrapper function to call the PLUMED evaluation routines
        and return forces."""

        if self.natoms != len(r["pos"]) / 3:
            raise ValueError("Size of atom array changed after initialization of FFPlumed")

        v = 0.0
        f = np.zeros(3 * self.natoms)
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
        self.plumed.cmd("setBox", r["cell"][0])
        self.plumed.cmd("setPositions", r["pos"])
        self.plumed.cmd("setForces", f)
        self.plumed.cmd("setVirial", vir)
        self.plumed.cmd("prepareCalc");
        self.plumed.cmd("performCalcNoUpdate");

        bias = np.zeros(1, float)
        self.plumed.cmd("getBias", bias)
        v = bias[0]

        r["result"] = [v, f, vir, ""]
        r["status"] = "Done"

    def mtd_update(self, pos, cell):
        """ Makes updates to the potential that only need to be triggered
        upon completion of a time step. """

        self.plumedstep += 1
        f = np.zeros(3 * self.natoms)
        vir = np.zeros((3, 3))

        self.plumed.cmd("setStep", self.plumedstep)
        self.plumed.cmd("setCharges", self.charges)
        self.plumed.cmd("setMasses", self.masses)
        self.plumed.cmd("setPositions", pos)
        self.plumed.cmd("setBox", cell)
        self.plumed.cmd("setForces", f)
        self.plumed.cmd("setVirial", vir)
        self.plumed.cmd("prepareCalc");
        self.plumed.cmd("performCalcNoUpdate");
        self.plumed.cmd("update")

        return True


class FFYaff(ForceField):

    """ Use Yaff as a library to construct a force field """

    def __init__(self, latency=1.0, name="", yaffpara=None, yaffsys=None, yafflog='yaff.log', rcut=18.89726133921252, alpha_scale=3.5, gcut_scale=1.1, skin=0, smooth_ei=False, reci_ei='ewald', pars=None, dopbc=False, threaded=True):
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
        super(FFYaff, self).__init__(latency, name, pars, dopbc)

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
        logf = open(yafflog, 'w')
        # Tell Python to close the file when the script exits
        atexit.register(logf.close)

        # Redirect Yaff log to file
        log._file = codecs.getwriter(locale.getpreferredencoding())(logf)

        self.system = System.from_file(self.yaffsys)
        self.ff = ForceField.generate(self.system, self.yaffpara, rcut=self.rcut, alpha_scale=self.alpha_scale, gcut_scale=self.gcut_scale, skin=self.skin, smooth_ei=self.smooth_ei, reci_ei=self.reci_ei)

        log._active = False

    def poll(self):
        """ Polls the forcefield checking if there are requests that should
        be answered, and if necessary evaluates the associated forces and energy. """

        # we have to be thread-safe, as in multi-system mode this might get called by many threads at once
        with self._threadlock:
            for r in self.requests:
                if r["status"] == "Queued":
                    r["status"] = "Running"
                    self.evaluate(r)

    def evaluate(self, r):
        """ Evaluate the energy and forces with the Yaff force field. """

        q = r["pos"]
        nat = len(q) / 3
        rvecs = r["cell"][0]

        self.ff.update_rvecs(np.ascontiguousarray(rvecs.T, dtype=np.float64))
        self.ff.update_pos(q.reshape((nat, 3)))
        gpos = np.zeros((nat, 3))
        vtens = np.zeros((3, 3))
        e = self.ff.compute(gpos, vtens)

        r["result"] = [e, -gpos.ravel(), -vtens, ""]
        r["status"] = "Done"
        r["t_finished"] = time.time()

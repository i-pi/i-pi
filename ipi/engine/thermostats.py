"""Classes that deal with constant temperature simulations.

Contains the algorithms which propagate the thermostatting steps in the constant
temperature ensembles. Includes the new GLE thermostat, which can be used to
run PI+GLE dynamics, reducing the number of path integral beads required.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np

from ipi.utils.depend import *
from ipi.utils.units import *
from ipi.utils.mathtools import matrix_exp, stab_cholesky, root_herm
from ipi.utils.prng import Random
from ipi.utils.messages import verbosity, warning, info
from ipi.engine.normalmodes import NormalModes


__all__ = [
    "Thermostat",
    "ThermoLangevin",
    "ThermoFFL",
    "ThermoPILE_L",
    "ThermoPILE_G",
    "ThermoSVR",
    "ThermoGLE",
    "ThermoNMGLE",
    "ThermoNMGLEG",
    "ThermoCL",
    "MultiThermo",
]


class Thermostat:
    """Base thermostat class.

    Gives the standard methods and attributes needed in all the thermostat
    classes.

    Attributes:
       prng: A pseudo random number generator object.
       ndof: The number of degrees of freedom that the thermostat will be
          attached to.

    Depend objects:
       dt: The time step used in the algorithms. Depends on the simulation dt.
       temp: The simulation temperature. Higher than the system temperature by
          a factor of the number of beads. Depends on the simulation temp.
       ethermo: The total energy exchanged with the bath due to the thermostat.
       p: The momentum vector that the thermostat is coupled to. Depends on the
          beads p object.
       m: The mass vector associated with p. Depends on the beads m object.
       sm: The square root of the mass vector.
    """

    def __init__(self, temp=1.0, dt=1.0, ethermo=0.0):
        """Initialises Thermostat.

        Args:
           temp: The simulation temperature. Defaults to 1.0.
           dt: The simulation time step. Defaults to 1.0.
           ethermo: The initial heat energy transferred to the bath.
              Defaults to 0.0. Will be non-zero if the thermostat is
              initialised from a checkpoint file.
        """

        self._temp = depend_value(name="temp", value=temp)
        self._dt = depend_value(name="dt", value=dt)
        self._ethermo = depend_value(name="ethermo", value=ethermo)

    def bind(self, beads=None, atoms=None, pm=None, nm=None, prng=None, fixdof=None):
        """Binds the appropriate degrees of freedom to the thermostat.

        This takes an object with degrees of freedom, and makes their momentum
        and mass vectors members of the thermostat. It also then creates the
        objects that will hold the data needed in the thermostat algorithms
        and the dependency network.

        Args:
           beads: An optional beads object to take the mass and momentum vectors
              from.
           atoms: An optional atoms object to take the mass and momentum vectors
              from.
           pm: An optional tuple containing a single momentum value and its
              conjugate mass.
           prng: An optional pseudo random number generator object. Defaults to
              Random().
           fixdof: An optional integer which can specify the number of constraints
              applied to the system. Defaults to zero.

        Raises:
           TypeError: Raised if no appropriate degree of freedom or object
              containing a momentum vector is specified for
              the thermostat to couple to.
        """

        if prng is None:
            warning(
                "Initializing thermostat from standard random PRNG", verbosity.medium
            )
            self.prng = Random()
        else:
            self.prng = prng

        if beads is not None:
            self._p = beads.p.flatten()
            self._m = beads.m3.flatten()
        elif atoms is not None:
            self._p = atoms._p
            self._m = atoms._m3
        elif pm is not None:
            self._p = pm[
                0
            ].flatten()  # MR this should allow to simply pass the cell momenta in the anisotropic barostat
            self._m = pm[1].flatten()
        else:
            raise TypeError(
                "Thermostat.bind expects either Beads, Atoms, NormalModes, or a (p,m) tuple to bind to"
            )

        if fixdof is None:
            self.ndof = len(self.p)
        else:
            self.ndof = float(len(self.p) - fixdof)

        self._sm = depend_array(
            name="sm",
            value=np.zeros(len(self._m)),
            func=self.get_sm,
            dependencies=[self._m],
        )

    def get_sm(self):
        """Retrieves the square root of the mass matrix.

        Returns:
           A vector of the square root of the mass matrix with one value for
           each degree of freedom.
        """

        return np.sqrt(self.m)

    def step(self):
        """Dummy thermostat step."""

        pass


dproperties(Thermostat, ["temp", "dt", "ethermo", "p", "m", "sm", "dt"])


class ThermoLangevin(Thermostat):
    """Represents a langevin thermostat.

    Depend objects:
       tau: Thermostat damping time scale. Larger values give a less strongly
          coupled thermostat.
       T: Coefficient of the diffusive contribution of the thermostat, i.e. the
          drift back towards equilibrium. Depends on tau and the time step.
       S: Coefficient of the stochastic contribution of the thermostat, i.e.
          the uncorrelated Gaussian noise. Depends on T and the temperature.
    """

    def get_T(self):
        """Calculates the coefficient of the overall drift of the velocities."""

        return np.exp(-self.dt / self.tau)

    def get_S(self):
        """Calculates the coefficient of the white noise."""
        return np.sqrt(Constants.kb * self.temp * (1 - self.T**2))

    def get_T_on_sm(self):
        """Calculates the combined mass scaling and thermostat damping."""
        return self.T / dstrip(self.sm)

    def __init__(self, temp=1.0, dt=1.0, tau=1.0, ethermo=0.0):
        """Initialises ThermoLangevin.

        Args:
           temp: The simulation temperature. Defaults to 1.0.
           dt: The simulation time step. Defaults to 1.0.
           tau: The thermostat damping timescale. Defaults to 1.0.
           ethermo: The initial heat energy transferred to the bath.
              Defaults to 0.0. Will be non-zero if the thermostat is
              initialised from a checkpoint file.
        """

        super(ThermoLangevin, self).__init__(temp, dt, ethermo)

        self._tau = depend_value(value=tau, name="tau")
        self._T = depend_value(
            name="T", func=self.get_T, dependencies=[self._tau, self._dt]
        )
        self._S = depend_value(
            name="S", func=self.get_S, dependencies=[self._temp, self._T]
        )

    def bind(self, beads=None, atoms=None, pm=None, nm=None, prng=None, fixdof=None):
        """Binds the appropriate degrees of freedom to the thermostat."""

        super(ThermoLangevin, self).bind(beads, atoms, pm, nm, prng, fixdof)

        self._T_on_sm = depend_value(
            name="T_on_sm", func=self.get_T_on_sm, dependencies=[self._T, self._sm]
        )

    def step(self):
        """Updates the bound momentum vector with a langevin thermostat."""

        # This performs the following steps
        # p <- p*exp(-dt/tau)+xi sqrt(C(1-exp-2dt/tau))
        # where dt is the thermostat time step (half of the MD step)
        # and C is the equilibrium fluctuations. to have a single
        # coefficient (and facilitate computing the kinetic energy)
        # the change in kinetic energy is computed to accumulate the work made
        # by the thermostat

        # goes in a single step to mass scaled coordinates and applies damping
        p = dstrip(self.p) * dstrip(self.T_on_sm)

        deltah = noddot(p, p) / (self.T**2)  # must correct for the "pre-damping"
        p += self.S * self.prng.gvec(len(p))  # random part (in ms coordinates)
        deltah -= noddot(p, p)  # new energy

        self.p[:] = p * dstrip(
            self.sm
        )  # back to physical momentum and updates actual p
        self.ethermo += deltah * 0.5


dproperties(ThermoLangevin, ["tau", "T", "S", "T_on_sm"])


class ThermoPILE_L(Thermostat):
    """Represents a PILE thermostat with a local centroid thermostat.

    Attributes:
       _thermos: The list of the different thermostats for all the ring polymer
          normal modes.
       nm: A normal modes object to attach the thermostat to.
       prng: Random number generator used in the stochastic integration
          algorithms.

    Depend objects:
       tau: Centroid thermostat damping time scale. Larger values give a
          less strongly coupled centroid thermostat.
       tauk: Thermostat damping time scale for the non-centroid normal modes.
          Depends on the ring polymer spring constant, and thus the simulation
          temperature.
       pilescale: A float used to reduce the intensity of the PILE thermostat if
          required.
    """

    def __init__(self, temp=1.0, dt=1.0, tau=1.0, ethermo=0.0, scale=1.0, pilect=0.0):
        """Initialises ThermoPILE_L.

        Args:
           temp: The simulation temperature. Defaults to 1.0.
           dt: The simulation time step. Defaults to 1.0.
           tau: The centroid thermostat damping timescale. Defaults to 1.0.
           ethermo: The initial conserved energy quantity. Defaults to 0.0. Will
              be non-zero if the thermostat is initialised from a checkpoint file.
           scale: A float used to reduce the intensity of the PILE thermostat if
              required.
           pilect: centroid mode temperature (if different from ensemble and set by input)

        Raises:
           TypeError: Raised if the thermostat is used with any object other than
              a beads object, so that we make sure that the objects needed for the
              normal mode transformation exist.
        """

        super(ThermoPILE_L, self).__init__(temp, dt, ethermo)
        self._tau = depend_value(value=tau, name="tau")
        self._pilescale = depend_value(value=scale, name="pilescale")
        self._pilect = depend_value(value=pilect, name="pilect")
        self._npilect = depend_value(
            func=self.get_npilect, name="npilect", dependencies=[self._pilect]
        )

    def bind(
        self,
        beads=None,
        atoms=None,
        pm=None,
        nm=None,
        prng=None,
        bindcentroid=True,
        fixdof=None,
    ):
        """Binds the appropriate degrees of freedom to the thermostat.

        This takes a beads object with degrees of freedom, and makes its momentum
        and mass vectors members of the thermostat. It also then creates the
        objects that will hold the data needed in the thermostat algorithms
        and the dependency network.

        Gives the interface for both the PILE_L and PILE_G thermostats, which
        only differ in their treatment of the centroid coordinate momenta.

        Args:
           nm: An optional normal mode object to take the mass and momentum
              vectors from.
           prng: An optional pseudo random number generator object. Defaults to
              Random().
           bindcentroid: An optional boolean which decides whether a Langevin
              thermostat is attached to the centroid mode of each atom
              separately, or the total kinetic energy. Defaults to True, which
              gives a thermostat bound to each centroid momentum.
           fixdof: An optional integer which can specify the number of constraints
              applied to the system. Defaults to zero.

        Raises:
           TypeError: Raised if no appropriate degree of freedom or object
              containing a momentum vector is specified for
              the thermostat to couple to.
        """

        if nm is None or not type(nm) is NormalModes:
            raise TypeError(
                "ThermoPILE_L.bind expects a NormalModes argument to bind to"
            )
        if prng is None:
            self.prng = Random()
        else:
            self.prng = prng

        prev_ethermo = self.ethermo

        # creates a set of thermostats to be applied to individual normal modes
        self._thermos = [ThermoLangevin(temp=1, dt=1, tau=1) for b in range(nm.nbeads)]
        # optionally does not bind the centroid, so we can re-use all of this
        # in the PILE_G case
        if not bindcentroid:
            self._thermos[0] = None

        self.nm = nm

        self._tauk = depend_array(
            name="tauk",
            value=np.zeros(nm.nbeads - 1, float),
            func=self.get_tauk,
            dependencies=[self._pilescale, nm._dynomegak],
        )

        # must pipe all the dependencies in such a way that values for the nm thermostats
        # are automatically updated based on the "main" thermostat
        def make_taugetter(k):
            return lambda: self.tauk[k - 1]

        it = 0
        for t in self._thermos:
            if t is None:
                it += 1
                continue
            if it > 0:
                fixdof = None  # only the centroid thermostat may have constraints

            # bind thermostat t to the it-th bead

            t.bind(pm=(nm.pnm[it, :], nm.dynm3[it, :]), prng=self.prng, fixdof=fixdof)
            if it == 0:
                # the following lines pipe a different temperature to the centroid, if requested
                if self.pilect > 0.0:
                    dpipe(self._npilect, t._temp)
                else:
                    dpipe(self._temp, t._temp)
            else:
                # pipes temp
                dpipe(self._temp, t._temp)
            # pipes dt
            dpipe(self._dt, t._dt)

            # for tau it is slightly more complex
            if it == 0:
                dpipe(self._tau, t._tau)
            else:
                # Here we manually connect _thermos[i].tau to tauk[i].
                # Simple and clear.
                t._tau.add_dependency(self._tauk)
                t._tau._func = make_taugetter(it)
            self._ethermo.add_dependency(t._ethermo)
            it += 1

        # since the ethermo will be "delegated" to the normal modes thermostats,
        # one has to split
        # any previously-stored value between the sub-thermostats
        if bindcentroid:
            for t in self._thermos:
                t.ethermo = prev_ethermo / nm.nbeads
            self._ethermo._func = self.get_ethermo
            # if we are not binding the centroid just yet, this bit of the piping
            # is delegated to the function which is actually calling this

    def get_tauk(self):
        """Computes the thermostat damping time scale for the non-centroid
        normal modes.

        Returns:
           An array with the damping time scales for the non-centroid modes.
        """

        # Also include an optional scaling factor to reduce the intensity of NM thermostats
        return np.array(
            [
                1.0 / (2 * self.pilescale * self.nm.dynomegak[k])
                for k in range(1, len(self._thermos))
            ]
        )

    def get_ethermo(self):
        """Computes the total energy transferred to the heat bath for all the
        thermostats.
        """

        et = 0.0
        for t in self._thermos:
            et += t.ethermo
        return et

    def get_npilect(self):
        """Multiplies centroid temperature by nbeads"""
        return self.nm.nbeads * self.pilect

    def step(self):
        """Updates the bound momentum vector with a PILE thermostat."""

        # super-cool! just loop over the thermostats! it's as easy as that!
        for t in self._thermos:
            t.step()


dproperties(ThermoPILE_L, ["tau", "pilescale", "pilect", "npilect", "tauk"])


class ThermoSVR(Thermostat):
    """Represents a stochastic velocity rescaling thermostat.

    Depend objects:
       tau: Centroid thermostat damping time scale. Larger values give a
          less strongly coupled centroid thermostat.
       K: Scaling factor for the total kinetic energy. Depends on the
          temperature.
       et: Parameter determining the strength of the thermostat coupling.
          Depends on tau and the time step.
    """

    def get_et(self):
        """Calculates the damping term in the propagator."""

        return np.exp(-self.dt / self.tau)

    def get_K(self):
        """Calculates the average kinetic energy per degree of freedom."""
        return Constants.kb * self.temp * 0.5

    def __init__(self, temp=1.0, dt=1.0, tau=1.0, ethermo=0.0):
        """Initialises ThermoSVR.

        Args:
           temp: The simulation temperature. Defaults to 1.0.
           dt: The simulation time step. Defaults to 1.0.
           tau: The thermostat damping timescale. Defaults to 1.0.
           ethermo: The initial conserved energy quantity. Defaults to 0.0. Will
              be non-zero if the thermostat is initialised from a checkpoint file.
        """

        super(ThermoSVR, self).__init__(temp, dt, ethermo)

        self._tau = depend_value(value=tau, name="tau")
        self._et = depend_value(
            name="et", func=self.get_et, dependencies=[self._tau, self._dt]
        )
        self._K = depend_value(name="K", func=self.get_K, dependencies=[self._temp])

    def step(self):
        """Updates the bound momentum vector with a stochastic velocity rescaling
        thermostat. See G Bussi, D Donadio, M Parrinello,
        Journal of Chemical Physics 126, 014101 (2007)
        """

        K = np.dot(dstrip(self.p), dstrip(self.p) / dstrip(self.m)) * 0.5

        # rescaling is un-defined if the KE is zero
        if K == 0.0:
            return

        # gets the stochastic term (basically a Gamma distribution for the kinetic energy)
        r1 = self.prng.g
        if (self.ndof - 1) % 2 == 0:
            rg = 2.0 * self.prng.gamma((self.ndof - 1) / 2)
        else:
            rg = 2.0 * self.prng.gamma((self.ndof - 2) / 2) + self.prng.g**2

        alpha2 = (
            self.et
            + self.K / K * (1 - self.et) * (r1**2 + rg)
            + 2.0 * r1 * np.sqrt(self.K / K * self.et * (1 - self.et))
        )
        alpha = np.sqrt(alpha2)
        if (r1 + np.sqrt(2 * K / self.K * self.et / (1 - self.et))) < 0:
            alpha *= -1

        self.ethermo += K * (1 - alpha2)
        self.p *= alpha


dproperties(ThermoSVR, ["tau", "et", "K"])


class ThermoPILE_G(ThermoPILE_L):
    """Represents a PILE thermostat with a global centroid thermostat.

    Simply replaces the Langevin thermostat for the centroid normal mode with
    a global velocity rescaling thermostat.
    """

    def __init__(self, temp=1.0, dt=1.0, tau=1.0, ethermo=0.0, scale=1.0, pilect=0.0):
        """Initialises ThermoPILE_G.

        Args:
           temp: The simulation temperature. Defaults to 1.0.
           dt: The simulation time step. Defaults to 1.0.
           tau: The centroid thermostat damping timescale. Defaults to 1.0.
           ethermo: The initial conserved energy quantity. Defaults to 0.0. Will
              be non-zero if the thermostat is initialised from a checkpoint file.
           scale: A float used to reduce the intensity of the PILE thermostat if
              required.
           pilect: centroid mode temperature (if different from ensemble and set by input).
              Default of 0.0 means it is not used and all modes have equal temperatures.
        """

        super(ThermoPILE_G, self).__init__(temp, dt, tau, ethermo)
        self._pilescale = depend_value(value=scale, name="pilescale")
        self._pilect = depend_value(value=pilect, name="pilect")
        self._npilect = depend_value(
            func=self.get_npilect, name="npilect", dependencies=[self._pilect]
        )

    def bind(self, beads=None, atoms=None, pm=None, nm=None, prng=None, fixdof=None):
        """Binds the appropriate degrees of freedom to the thermostat.

        This takes a beads object with degrees of freedom, and makes its momentum
        and mass vectors members of the thermostat. It also then creates the
        objects that will hold the data needed in the thermostat algorithms
        and the dependency network.

        Uses the PILE_L bind interface, with bindcentroid set to false so we can
        specify that thermostat separately, by binding a global
        thermostat to the centroid mode.

        Args:
           beads: An optional beads object to take the mass and momentum vectors
              from.
           prng: An optional pseudo random number generator object. Defaults to
              Random().
           fixdof: An optional integer which can specify the number of constraints
              applied to the system. Defaults to zero.

        """

        # first binds as a local PILE, then substitutes the thermostat on the centroid
        prev_ethermo = self.ethermo
        super(ThermoPILE_G, self).bind(
            nm=nm, prng=prng, bindcentroid=False, fixdof=fixdof
        )

        # centroid thermostat
        self._thermos[0] = ThermoSVR(temp=1, dt=1, tau=1)

        t = self._thermos[0]
        t.bind(pm=(nm.pnm[0, :], nm.dynm3[0, :]), prng=self.prng, fixdof=fixdof)
        # the next lines pipe a different temperatures to the centroid modes, if requested.
        if self.pilect > 0.0:
            dpipe(self._npilect, t._temp)
        else:
            dpipe(self._temp, t._temp)

        dpipe(self._dt, t._dt)
        dpipe(self._tau, t._tau)
        self._ethermo.add_dependency(t._ethermo)

        # splits any previous ethermo between the thermostats, and finishes to bind ethermo to the sum function
        for t in self._thermos:
            t.ethermo = prev_ethermo / nm.nbeads
        self._ethermo._func = self.get_ethermo


class ThermoGLE(Thermostat):
    """Represents a generalized Langevin equation thermostat.

    This is similar to a langevin thermostat, in that it uses Gaussian random
    numbers to simulate a heat bath acting on the system, but simulates a
    non-Markovian system by using a Markovian formulation in an extended phase
    space. This allows for a much greater degree of flexibility, and this
    thermostat, properly fitted, can give an approximation to the correct
    quantum ensemble even for a classical, 1-bead simulation. More reasonably,
    using this thermostat allows for a far smaller number of replicas of the
    system to be used, as the convergence of the properties
    of the system is accelerated with respect to number of beads when PI+GLE
    are used in combination. (See M. Ceriotti, D. E. Manolopoulos, M. Parinello,
    J. Chem. Phys. 134, 084104 (2011)).

    Attributes:
       ns: The number of auxilliary degrees of freedom.
       s: An array holding all the momenta, including the ones for the
          auxilliary degrees of freedom.

    Depend objects:
       A: Drift matrix giving the damping time scales for all the different
          degrees of freedom.
       C: Static covariance matrix.
          Satisfies A.C + C.transpose(A) = B.transpose(B), where B is the
          diffusion matrix, giving the strength of the coupling of the system
          with the heat bath, and thus the size of the stochastic
          contribution of the thermostat.
       T: Matrix for the diffusive contribution of the thermostat, i.e. the
          drift back towards equilibrium. Depends on A and the time step.
       S: Matrix for the stochastic contribution of the thermostat, i.e.
          the uncorrelated Gaussian noise. Depends on C and T.
    """

    def get_T(self):
        """Calculates the matrix for the overall drift of the velocities."""

        return matrix_exp(-self.dt * self.A)

    def get_S(self):
        """Calculates the matrix for the coloured noise."""

        SST = Constants.kb * (self.C - np.dot(self.T, np.dot(self.C, self.T.T)))
        # Uses a symetric decomposition rather than Cholesky, since it is more stable
        return root_herm(SST)
        # return stab_cholesky(SST)

    def get_C(self):
        """Calculates C from temp (if C is not set explicitly)"""

        rC = np.identity(self.ns + 1, float) * self.temp
        return rC[:]

    def __init__(self, temp=1.0, dt=1.0, A=None, C=None, ethermo=0.0):
        """Initialises ThermoGLE.

        Args:
           temp: The simulation temperature. Defaults to 1.0.
           dt: The simulation time step. Defaults to 1.0.
           A: An optional matrix giving the drift matrix. Defaults to a single
              value of 1.0.
           C: An optional matrix giving the covariance matrix. Defaults to an
              identity matrix times temperature with the same dimensions as the
              total number of degrees of freedom in the system.
           ethermo: The initial heat energy transferred to the bath.
              Defaults to 0.0. Will be non-zero if the thermostat is
              initialised from a checkpoint file.
        """

        super(ThermoGLE, self).__init__(temp, dt, ethermo)

        if A is None:
            A = np.identity(1, float)
        self._A = depend_value(value=A.copy(), name="A")

        self.ns = len(self.A) - 1

        # now, this is tricky. if C is taken from temp, then we want it to be updated
        # as a depend of temp. Otherwise, we want it to be an independent beast.
        if C is None:
            C = np.identity(self.ns + 1, float) * self.temp
            self._C = depend_value(name="C", func=self.get_C, dependencies=[self._temp])
        else:
            self._C = depend_value(value=C.copy(), name="C")

        self._T = depend_value(
            name="T", func=self.get_T, dependencies=[self._A, self._dt]
        )
        self._S = depend_value(
            name="S", func=self.get_S, dependencies=[self._C, self._T]
        )

        self.s = np.zeros(0)

    def bind(self, beads=None, atoms=None, pm=None, nm=None, prng=None, fixdof=None):
        """Binds the appropriate degrees of freedom to the thermostat.

        This takes an object with degrees of freedom, and makes their momentum
        and mass vectors members of the thermostat. It also then creates the
        objects that will hold the data needed in the thermostat algorithms
        and the dependency network.

        Args:
           beads: An optional beads object to take the mass and momentum vectors
              from.
           atoms: An optional atoms object to take the mass and momentum vectors
              from.
           pm: An optional tuple containing a single momentum value and its
              conjugate mass.
           prng: An optional pseudo random number generator object. Defaults to
              Random().
           fixdof: An optional integer which can specify the number of constraints
              applied to the system. Defaults to zero.

        Raises:
           TypeError: Raised if no appropriate degree of freedom or object
              containing a momentum vector is specified for
              the thermostat to couple to.
        """

        super(ThermoGLE, self).bind(
            beads=beads, atoms=atoms, pm=pm, prng=prng, fixdof=fixdof
        )

        # allocates, initializes or restarts an array of s's
        if self.s.shape != (self.ns + 1, len(self._m)):
            if len(self.s) > 0:
                warning(
                    "Mismatch in GLE s array size on restart, will reinitialise to free particle.",
                    verbosity.low,
                )
            self.s = np.zeros((self.ns + 1, len(self._m)))

            # Initializes the s vector in the free-particle limit
            info(
                " GLE additional DOFs initialised to the free-particle limit.",
                verbosity.low,
            )
            SC = stab_cholesky(self.C * Constants.kb)
            self.s[:] = np.dot(SC, self.prng.gvec(self.s.shape))
        else:
            info("GLE additional DOFs initialised from input.", verbosity.medium)

    def step(self):
        """Updates the bound momentum vector with a GLE thermostat"""

        self.s[0, :] = self.p / self.sm

        self.ethermo += np.dot(self.s[0], self.s[0]) * 0.5
        self.s[:] = np.dot(self.T, self.s) + np.dot(
            self.S, self.prng.gvec(self.s.shape)
        )
        self.ethermo -= np.dot(self.s[0], self.s[0]) * 0.5

        self.p = self.s[0] * self.sm


dproperties(ThermoGLE, ["A", "C", "T", "S"])


class ThermoNMGLE(Thermostat):
    """Represents a 'normal-modes' generalized Langevin equation thermostat.

    An extension to the GLE thermostat which is applied in the
    normal modes representation, and which allows to use a different
    GLE for each normal mode

    Attributes:
       ns: The number of auxilliary degrees of freedom.
       nb: The number of beads.
       s: An array holding all the momenta, including the ones for the
          auxilliary degrees of freedom.

    Depend objects:
       A: Drift matrix giving the damping time scales for all the different
          degrees of freedom (must contain nb terms).
       C: Static covariance matrix.
          Satisfies A.C + C.transpose(A) = B.transpose(B), where B is the
          diffusion matrix, giving the strength of the coupling of the system
          with the heat bath, and thus the size of the stochastic
          contribution of the thermostat.
    """

    def get_C(self):
        """Calculates C from temp (if C is not set explicitely)."""

        rv = np.ndarray((self.nb, self.ns + 1, self.ns + 1), float)
        for b in range(0, self.nb):
            rv[b] = np.identity(self.ns + 1, float) * self.temp
        return rv[:]

    def __init__(self, temp=1.0, dt=1.0, A=None, C=None, ethermo=0.0):
        """Initialises ThermoGLE.

        Args:
           temp: The simulation temperature. Defaults to 1.0.
           dt: The simulation time step. Defaults to 1.0.
           A: An optional matrix giving the drift matrix. Defaults to a single
              value of 1.0.
           C: An optional matrix giving the covariance matrix. Defaults to an
              identity matrix times temperature with the same dimensions as the
              total number of degrees of freedom in the system.
           ethermo: The initial heat energy transferred to the bath.
              Defaults to 0.0. Will be non-zero if the thermostat is
              initialised from a checkpoint file.
        """

        super(ThermoNMGLE, self).__init__(temp, dt, ethermo)

        if A is None:
            A = np.identity(1, float)
        self._A = depend_value(value=A.copy(), name="A")

        self.nb = len(self.A)
        self.ns = len(self.A[0]) - 1

        # now, this is tricky. if C is taken from temp, then we want it to be
        # updated as a depend of temp.
        # Otherwise, we want it to be an independent beast.
        if C is None:
            self._C = depend_value(name="C", func=self.get_C, dependencies=[self._temp])
        else:
            self._C = depend_value(value=C.copy(), name="C")

    def bind(self, beads=None, atoms=None, pm=None, nm=None, prng=None, fixdof=None):
        """Binds the appropriate degrees of freedom to the thermostat.

        This takes an object with degrees of freedom, and makes their momentum
        and mass vectors members of the thermostat. It also then creates the
        objects that will hold the data needed in the thermostat algorithms
        and the dependency network. Actually, this specific thermostat requires
        being called on a beads object.

        Args:
           nm: An optional normal modes object to take the mass and momentum
              vectors from.
           prng: An optional pseudo random number generator object. Defaults to
              Random().
           fixdof: An optional integer which can specify the number of constraints
              applied to the system. Defaults to zero.

        Raises:
           TypeError: Raised if no beads object is specified for
              the thermostat to couple to.
        """

        if nm is None or not type(nm) is NormalModes:
            raise TypeError(
                "ThermoNMGLE.bind expects a NormalModes argument to bind to"
            )

        if prng is None:
            self.prng = Random()
        else:
            self.prng = prng

        if nm.nbeads != self.nb:
            raise IndexError(
                "The parameters in nm_gle options correspond to a bead number "
                + str(self.nb)
                + " which does not match the number of beads in the path"
                + str(nm.nbeads)
            )

        # allocates, initializes or restarts an array of s's
        if self.s.shape != (self.nb, self.ns + 1, nm.natoms * 3):
            if len(self.s) > 0:
                warning(
                    "Mismatch in GLE s array size on restart, will reinitialise to free particle.",
                    verbosity.low,
                )
            self.s = np.zeros((self.nb, self.ns + 1, nm.natoms * 3))

            # Initializes the s vector in the free-particle limit
            info(
                " GLE additional DOFs initialised to the free-particle limit.",
                verbosity.low,
            )
            for b in range(self.nb):
                SC = stab_cholesky(self.C[b] * Constants.kb)
                self.s[b] = np.dot(SC, self.prng.gvec(self.s[b].shape))
        else:
            info("GLE additional DOFs initialised from input.", verbosity.medium)

        prev_ethermo = self.ethermo

        # creates a set of thermostats to be applied to individual normal modes
        self._thermos = [
            ThermoGLE(temp=1, dt=1, A=self.A[b], C=self.C[b]) for b in range(self.nb)
        ]

        # must pipe all the dependencies in such a way that values for the nm
        # thermostats are automatically updated based on the "main" thermostat
        def make_Agetter(k):
            return lambda: self.A[k]

        def make_Cgetter(k):
            return lambda: self.C[k]

        it = 0
        for t in self._thermos:
            t.s = self.s[it]  # gets the s's as a slice of self.s
            t.bind(
                pm=(nm.pnm[it, :], nm.dynm3[it, :]), prng=self.prng
            )  # bind thermostat t to the it-th normal mode

            # pipes temp and dt
            dpipe(self._temp, t._temp)
            dpipe(self._dt, t._dt)

            # here we pipe the A and C of individual NM to the "main" arrays
            t._A.add_dependency(self._A)
            t._A._func = make_Agetter(it)
            t._C.add_dependency(self._C)
            t._C._func = make_Cgetter(it)
            self._ethermo.add_dependency(t._ethermo)
            it += 1

        # since the ethermo will be "delegated" to the normal modes thermostats,
        # one has to split
        # any previously-stored value between the sub-thermostats
        for t in self._thermos:
            t.ethermo = prev_ethermo / self.nb

        self._ethermo._func = self.get_ethermo

    def step(self):
        """Updates the thermostat in NM representation by looping over the
        individual DOFs.
        """

        for t in self._thermos:
            t.step()

    def get_ethermo(self):
        """Computes the total energy transferred to the heat bath for all the nm
        thermostats.
        """

        et = 0.0
        for t in self._thermos:
            et += t.ethermo
        return et


dproperties(ThermoNMGLE, ["A", "C", "T", "S"])


class ThermoNMGLEG(ThermoNMGLE):
    """Represents a 'normal-modes' generalized Langevin equation thermostat + SVR.

    An extension to the above NMGLE thermostat which also adds a stochastic
    velocity rescaling to the centroid. Allows kinetic energy as well as
    potential energy sampling optimization.

    Depend objects:
       tau: Thermostat damping time scale. Larger values give a less strongly
          coupled thermostat.
    """

    def __init__(self, temp=1.0, dt=1.0, A=None, C=None, tau=1.0, ethermo=0.0):
        super(ThermoNMGLEG, self).__init__(temp, dt, A, C, ethermo)
        self._tau = depend_value(value=tau, name="tau")

    def bind(self, beads=None, atoms=None, pm=None, nm=None, prng=None, fixdof=None):
        """Binds the appropriate degrees of freedom to the thermostat.

        This takes an object with degrees of freedom, and makes their momentum
        and mass vectors members of the thermostat. It also then creates the
        objects that will hold the data needed in the thermostat algorithms
        and the dependency network. Actually, this specific thermostat requires
        being called on a beads object.

        Args:
           nm: An optional normal modes object to take the mass and momentum
              vectors from.
           prng: An optional pseudo random number generator object. Defaults to
              Random().
           fixdof: An optional integer which can specify the number of constraints
              applied to the system. Defaults to zero.
        """

        super(ThermoNMGLEG, self).bind(nm=nm, prng=prng, fixdof=fixdof)

        t = ThermoSVR(self.temp, self.dt, self.tau)

        t.bind(
            pm=(nm.pnm[0, :], nm.dynm3[0, :]), prng=self.prng
        )  # bind global thermostat to centroid

        # pipes temp and dt
        dpipe(self._temp, t._temp)
        dpipe(self._dt, t._dt)
        dpipe(self._tau, t._tau)

        self._ethermo.add_dependency(t._ethermo)
        self._thermos.append(t)


dproperties(ThermoNMGLE, ["tau"])


class ThermoCL(Thermostat):
    """Represents a Langevin thermostat for systems driven by forces which are statistical
       in nature, i.e. they contain noise, and/or have an unwanted dissipative effect.
       This thermostat's dissipative term is modified to compensate for the inherent noise.
       Similarly, the noise term is modified to compensate for the inherent dissipation.
       The time scales of inherent diffusion and drift (intau and idtau) must be set
       to adequate values in order to make the system reach the target temperature.
       Alternatively, if the automatic parameter adjustment time scale apat is set > 0
       and either intau OR idtau are initially set to a value > 0, the respective parameter
       will be automatically adjusted over time, driven by current versus target temperature.

       Author: Jan Kessler <jakessle@t-online.de>

    Depend objects:
       tau: Thermostat damping time scale. Larger values give a less strongly
          coupled thermostat.
       intau: Inherent noise time scale. Smaller intau results in more compensating drift.
       idtau: Inherent dissipation time scale. Smaller idtau results in more compensating diffusion.
       apat: Automatic parameter adjustment time scale.
       lgT: Langevin coefficient of drift.
       inT: Coefficient of drift compensating for inherent noise.
       idT: Coefficient of drift of inherent dissipation.
       T: Coefficient of the drift contribution of the thermostat, i.e. the
          dissipation back towards equilibrium. Depends on lgT and inT.
       S: Coefficient of the diffusive contribution of the thermostat, i.e.
          the uncorrelated Gaussian noise. Depends on lgT, idT and the temperature.
    """

    def get_lgT(self):
        """Calculates the Langevin coefficient of drift."""

        if self.tau > 0:
            return np.exp(-0.5 * self.dt / self.tau)
        else:
            return 1.0

    def get_inT(self):
        """Calculates the coefficient of drift to compensate for the inherent noise."""

        if self.intau > 0:
            return np.exp(-0.5 * self.dt / self.intau)
        else:
            return 1.0

    def get_idT(self):
        """Calculates the coefficient of drift for the inherent dissipation."""

        if self.idtau > 0:
            return np.exp(-0.5 * self.dt / self.idtau)
        else:
            return 1.0

    def get_T(self):
        """Calculates the coefficient of the thermosta'st drift of the velocities."""

        return self.lgT * self.inT

    def get_S(self):
        """Calculates the coefficient of the white noise."""

        return np.sqrt(Constants.kb * self.temp * (1.0 - (self.lgT * self.idT) ** 2))

    def __init__(self, temp=1.0, dt=1.0, tau=0, intau=0, idtau=0, apat=0, ethermo=0.0):
        """Initialises ThermoCL.

        Args:
           temp: The simulation temperature. Defaults to 1.0.
           dt: The simulation time step. Defaults to 1.0.
           tau: The thermostat damping timescale. Defaults to 0 (off).
           intau: Estimated inherent noise time coefficient. Defaults to 0 (off).
           idtau: Estimated inherent dissipation time coefficient. Defaults to 0 (off).
           apat: Automatic parameter adjustment time scale. Defaults to 0 (off).
           ethermo: The initial heat energy transferred to the bath.
              Defaults to 0.0. Will be non-zero if the thermostat is
              initialised from a checkpoint file.
        """

        super(ThermoCL, self).__init__(temp, dt, ethermo)

        self.idstep = False
        self._tau = depend_value(value=tau, name="tau")
        self._intau = depend_value(value=intau, name="intau")
        self._idtau = depend_value(value=idtau, name="idtau")
        self._apat = depend_value(value=apat, name="apat")
        self._lgT = depend_value(
            name="lgT", func=self.get_lgT, dependencies=[self._tau, self._dt]
        )
        self._inT = depend_value(
            name="inT", func=self.get_inT, dependencies=[self._intau, self._dt]
        )
        self._idT = depend_value(
            name="idT", func=self.get_idT, dependencies=[self._idtau, self._dt]
        )
        self._T = depend_value(
            name="T", func=self.get_T, dependencies=[self._lgT, self._inT]
        )
        self._S = depend_value(
            name="S", func=self.get_S, dependencies=[self._temp, self._lgT, self._idT]
        )

    def step(self):
        """Updates the bound momentum vector with a langevin thermostat."""

        et = self.ethermo
        p = dstrip(self.p).copy()
        sm = dstrip(self.sm)

        p /= sm

        et += np.dot(p, p) * 0.5
        p *= self.T
        p += self.S * self.prng.gvec(len(p))
        et -= np.dot(p, p) * 0.5

        p *= sm

        self.p = p
        self.ethermo = et

        if self.apat > 0 and self.idstep and ((self.intau != 0) ^ (self.idtau != 0)):
            ekin = np.dot(dstrip(self.p), dstrip(self.p) / dstrip(self.m)) * 0.5
            mytemp = ekin / Constants.kb / self.ndof * 2

            if self.intau != 0:
                if mytemp != 0:
                    self.intau /= (mytemp / self.temp) ** (self.dt / self.apat)
                print(("ThermoCL inherent noise time scale: " + str(self.intau)))
            else:
                self.idtau *= (mytemp / self.temp) ** (self.dt / self.apat)
                print(("ThermoCL inherent dissipation time scale: " + str(self.idtau)))

        self.idstep = not self.idstep


dproperties(ThermoCL, ["tau", "T", "S", "idtau", "intau", "apat", "lgT", "inT", "idT"])


class ThermoFFL(ThermoLangevin):
    """Represents a fast-forward langevin thermostat.

    Depend objects:
       tau: Thermostat damping time scale. Larger values give a less strongly
          coupled thermostat.
       T: Coefficient of the diffusive contribution of the thermostat, i.e. the
          drift back towards equilibrium. Depends on tau and the time step.
       S: Coefficient of the stochastic contribution of the thermostat, i.e.
          the uncorrelated Gaussian noise. Depends on T and the temperature.
       flip: Type of flip to use ('soft', 'hard', 'rescale').
    """

    def get_T(self):
        """Calculates the coefficient of the overall drift of the velocities."""

        return np.exp(-self.dt / self.tau)

    def get_S(self):
        """Calculates the coefficient of the white noise."""

        return np.sqrt(Constants.kb * self.temp * (1 - self.T**2))

    def __init__(self, temp=1.0, dt=1.0, tau=1.0, ethermo=0.0, flip="rescale"):
        """Initialises ThermoFFL.

        Args:
           temp: The simulation temperature. Defaults to 1.0.
           dt: The simulation time step. Defaults to 1.0.
           tau: The thermostat damping timescale. Defaults to 1.0.
           ethermo: The initial heat energy transferred to the bath.
              Defaults to 0.0. Will be non-zero if the thermostat is
              initialised from a checkpoint file.
           flip: The flipping type. Defaults to 'rescale'.
              Allowed values are 'soft', 'hard', 'rescale', 'none'.
        """

        super(ThermoFFL, self).__init__(temp, dt, ethermo)

        self._dt = depend_value(value=dt, name="dt")
        self._tau = depend_value(value=tau, name="tau")
        self._T = depend_value(
            name="T", func=self.get_T, dependencies=[self._tau, self._dt]
        )
        self._S = depend_value(
            name="S", func=self.get_S, dependencies=[self._temp, self._T]
        )
        self._flip = depend_value(value=flip, name="flip")

        allowed_flips = ["soft", "hard", "rescale", "none"]
        if not (self.flip in allowed_flips):
            raise KeyError(
                "Invalid flip type "
                + self.flip
                + '; allowed flip types: "'
                + '", "'.join(allowed_flips)
                + '"'
            )

    def step(self):
        """Updates the bound momentum vector with a fast-forward langevin thermostat."""

        et = self.ethermo
        p = dstrip(self.p).copy()
        sm = dstrip(self.sm)

        p /= sm

        # Store momentum before langevin step
        p_old = list(p)

        # Accumulate conserved quantity
        et += np.dot(p, p) * 0.5

        # Do standard langevin thermostatting
        p *= self.T
        p += self.S * self.prng.gvec(len(p))

        # Check whether to flip momenta back
        if self.flip == "soft":
            # Soft flip
            p_old = np.reshape(p_old, (len(p) / 3, 3))
            p_new = np.reshape(p, (len(p) / 3, 3))
            dotpr = hfunc(
                np.sum(np.multiply(p_old, p_new), axis=1)
                / np.sum(np.multiply(p_old, p_old), axis=1)
            )
            p += np.reshape(np.multiply(dotpr, p_old.T).T, len(p))
        elif self.flip == "hard":
            # Hard flip
            p = np.multiply(p, np.sign(np.multiply(p, p_old)))
        elif self.flip == "rescale":
            # Rescale flip
            p_old = np.reshape(p_old, (len(p) / 3, 3))
            p_new = np.reshape(p, (len(p) / 3, 3))
            scalfac = np.linalg.norm(p_new, axis=1) / np.linalg.norm(p_old, axis=1)
            p = np.reshape(np.multiply(scalfac, p_old.T).T, len(p))
        # Otherwise we have chosen 'none', and we just don't do anything here

        # Accumulate conserved quantity
        et -= np.dot(p, p) * 0.5

        p *= sm

        self.p = p
        self.ethermo = et


dproperties(ThermoFFL, "flip")


class MultiThermo(Thermostat):
    def __init__(self, temp=1.0, dt=1.0, ethermo=0.0, thermolist=None):
        """Initialises Thermostat.

        Args:
           temp: The simulation temperature. Defaults to 1.0.
           dt: The simulation time step. Defaults to 1.0.
           ethermo: The initial heat energy transferred to the bath.
              Defaults to 0.0. Will be non-zero if the thermostat is
              initialised from a checkpoint file.
        """

        self.tlist = thermolist
        self._temp = depend_value(name="temp", value=temp)
        self._dt = depend_value(name="dt", value=dt)
        self._ethermo = depend_value(name="ethermo", value=ethermo)
        for t in self.tlist:
            dpipe(self._dt, t._dt)
            dpipe(self._temp, t._temp)

    def get_ethermo(self):
        et = 0.0
        for t in self.tlist:
            et += t.ethermo
        return et

    def bind(self, beads=None, atoms=None, pm=None, nm=None, prng=None, fixdof=None):
        """Binds the appropriate degrees of freedom to the thermostat."""

        # just binds all the sub-thermostats
        for t in self.tlist:
            t.bind(beads=beads, atoms=atoms, pm=pm, nm=nm, prng=prng, fixdof=fixdof)
            self._ethermo.add_dependency(t._ethermo)

        self._ethermo._func = self.get_ethermo

    def step(self):
        """Steps through all sub-thermostats."""

        for t in self.tlist:
            t.step()
        pass


def hfunc(x):
    return (np.sign(x) - 1.0) * x

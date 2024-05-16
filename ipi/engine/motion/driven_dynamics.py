"""Creates objects to deal with dynamics driven by external fields."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

import numpy as np

from ipi.utils.depend import *
from ipi.utils.units import Constants
from ipi.engine.motion.dynamics import NVEIntegrator, DummyIntegrator, Dynamics
from ipi.utils.units import UnitMap
import re


class DrivenDynamics(Dynamics):
    """self (path integral) molecular dynamics class.

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

    driven_integrators = ["eda-nve"]

    def __init__(
        self,
        efield=None,
        bec=None,
        *argc,
        **argv,
    ):
        """Initialises a "dynamics" motion object.

        Args:
            dt: The timestep of the simulation algorithms.
            fixcom: An optional boolean which decides whether the centre of mass
                motion will be constrained or not. Defaults to False.
        """

        super().__init__(*argc, **argv)

        if self.enstype == "eda-nve":
            # NVE integrator with an external time-dependent electric field
            self.integrator = EDANVEIntegrator()
        else:
            self.integrator = DummyIntegrator()

        self.Electric_Field = efield
        self.Born_Charges = bec
        self.Electric_Dipole = ElectricDipole()
        self._time = depend_value(name="time", value=0)

    def bind(self, ens, beads, nm, cell, bforce, prng, omaker):
        dpipe(dfrom=ens._time, dto=self._time)
        super().bind(ens, beads, nm, cell, bforce, prng, omaker)
        self.Electric_Dipole.bind(ens)
        self.Electric_Field.bind(self, self.enstype)
        self.Born_Charges.bind(ens, self.enstype)


dproperties(DrivenDynamics, ["dt", "nmts", "splitting", "ntemp", "time", "efield_time"])


class EDAIntegrator(DummyIntegrator):
    """Integrator object for simulations using the Electric Dipole Approximation (EDA)
    when an external electric field is applied.
    """

    def __init__(self):
        super().__init__()
        self._time = depend_value(name="time", value=0)

    def bind(self, motion: DrivenDynamics):
        """bind variables"""
        super().bind(motion)
        dpipe(dfrom=motion._time, dto=self._time)
        self.Born_Charges: BEC = motion.Born_Charges
        self.Electric_Field: ElectricField = motion.Electric_Field

    def td_pstep(self, time, level=0):
        """Velocity Verlet momentum propagator with a time dependent force."""
        # This method adds the time-dependent force contribution to the momenta.
        # and it is called only twice per MD step, as implemented in `EDANVEIntegrator.step`.

        if level != 0:
            # A time-dependent force integrator is ill-defined in a MTS algorithm framework.
            # For this reason, the time-dependent contribution can added in the outermost layer only,
            # such that to avoid any inconsistency in the definition of the time used to evaluate the forces.
            raise ValueError(
                "EDAIntegrator can be called only in the outermost layer of a Multiple Time Step algorithm."
            )

        edaforces = self.EDAforces(time)
        self.beads.p += edaforces * self.pdt[level]

        pass

    def EDAforces(self, time):
        """Compute the EDA contribution to the forces, i.e. `q_e Z^* @ E(t)`"""
        Z = dstrip(self.Born_Charges.bec)  # tensor of shape (nbeads,3xNatoms,3)
        E = self.Electric_Field.Efield(time)  # vector of shape (3)
        forces = Constants.e * Z @ E  # array of shape (nbeads,3xNatoms)
        return forces


dproperties(EDAIntegrator, ["time"])


class EDANVEIntegrator(EDAIntegrator, NVEIntegrator):
    """Integrator object for simulations with constant Number of particles, Volume, and Energy (NVE)
    using the Electric Dipole Approximation (EDA) when an external electric field is applied.
    """

    def step(self, step):
        time = self.time
        EDAIntegrator.td_pstep(self, time, 0)
        NVEIntegrator.step(self, step)
        time = self.time + self.dt
        EDAIntegrator.td_pstep(self, time, 0)


class BEC:
    """Class to handle the Born Effective Charge tensors when performing driven dynamics (with 'eda-nve')"""

    # The BEC tensors Z^* are defined as the derivative of the electric dipole of the system w.r.t. nuclear positions
    # in units of the elementary charge, i.e. Z^*_{ij} = 1/q_e \frac{\partial d_i }{\partial R_j}
    # The dipole is a vector of shape (3), while the nuclear positions have shape (3xNatoms)
    # The BEC tensors have then shape (3xNatoms,3).
    # A tensor of this shape is stored for each bead --> self._bec.shape = (self.nbeads, 3 * self.natoms, 3)
    #
    # If an external time-dependent electric field E(t) is applied, this couples to the dipole of the system,
    # and the resulting additional term to the forces is given by q_e Z^* @ E(t) --> have a look at EDAIntegrator._eda_forces
    #
    # The BEC tensors Z^* can be given to i-PI by an external driver through the etxras strings in the forces
    # of they can be kept fixed during the dynamics: in this case you can provide them through a txt file.

    def __init__(self, cbec=None, bec=None):
        self.cbec = cbec
        if bec is None:
            bec = np.full((0, 3), np.nan)
        self._bec = depend_array(name="bec", value=bec)
        pass

    def bind(self, ensemble, enstype):
        self.enstype = enstype
        self.nbeads = ensemble.beads.nbeads
        self.natoms = ensemble.beads.natoms
        self.forces = ensemble.forces

        if self.enstype in DrivenDynamics.driven_integrators and self.cbec:
            self._bec = depend_array(
                name="bec",
                value=np.full((self.nbeads, 3 * self.natoms, 3), np.nan),
                func=self._get_driver_BEC,
                dependencies=[ensemble.beads._q],
            )
        elif self.enstype in DrivenDynamics.driven_integrators:
            temp = self._get_fixed_BEC()  # reshape the BEC once and for all
            self._bec = depend_array(name="bec", value=temp)
        else:
            self._bec = depend_array(
                name="bec", value=np.full((self.nbeads, 3 * self.natoms, 3), np.nan)
            )

        pass

    def store(self, bec):
        super(BEC, self).store(bec)
        self.cbec.store(bec.cbec)

    def _get_driver_BEC(self, bead=None):
        """Return the BEC tensors (in cartesian coordinates), when computed by the driver"""

        msg = "Error in '_get_driver_BEC'"

        if bead is not None:
            if bead < 0:
                raise ValueError("Error in '_get_driver_BEC': 'bead' is negative")
            if bead >= self.nbeads:
                raise ValueError(
                    "Error in '_get_driver_BEC': 'bead' is greater than the number of beads"
                )
        else:
            if self.nbeads != 1:
                raise ValueError(
                    "Error in '_get_driver_BEC': EDA integration has not implemented yet for 'nbeads' > 1"
                )

        if self.cbec:
            if "BEC" not in self.forces.extras:
                raise ValueError(
                    msg
                    + ": BEC tensors are not returned to i-PI (or at least not accessible in '_get_driver_BEC')."
                )
        else:
            raise ValueError(
                msg + ": you should not get into this functon if 'cbec' is False."
            )

        BEC = np.full((self.nbeads, 3 * self.natoms, 3), np.nan)
        for n in range(self.nbeads):
            bec = np.asarray(self.forces.extras["BEC"][n])

            if bec.shape[0] != 3 * self.natoms:
                raise ValueError(
                    msg
                    + ": number of BEC tensors is not equal to the number fo atoms x 3."
                )
            if bec.shape[1] != 3:
                raise ValueError(
                    msg
                    + ": BEC tensors with wrong shape. They should have 3 components."
                )

            BEC[n, :, :] = np.copy(bec)

        return BEC

    def _get_fixed_BEC(self):
        """Return the BEC tensors (in cartesian coordinates).
        The BEC tensor are stored in a compact form.
        This method trasform the BEC tensors into another data structure, suitable for computation.
        A lambda function is also returned to perform fast matrix multiplication.
        """
        try:
            return self.bec.reshape((self.nbeads, 3 * self.natoms, 3))
        except:
            line = (
                "Error in '_get_fixed_BEC': i-PI is going to stop.\n"
                + "The BEC tensor is: "
                + str(self.bec)
                + "\nPay attention that in 'input.xml' you should have the following line:\n"
                + "\t'<bec mode=\"file\"> filepath </bec>'\n"
                + "The default mode could be 'none'. Please change it."
            )
            raise ValueError(line)


dproperties(BEC, ["bec"])


class ElectricDipole:
    """Class to handle the electric dipole of the system when performing driven dynamics (with 'eda-nve')"""

    def __init__(self):
        pass

    def bind(self, ensemble):
        self._nbeads = depend_value(name="nbeads", value=ensemble.beads.nbeads)

        self.ens = ensemble

        val = np.full((self.nbeads, 3), np.nan)
        self._dipole = depend_array(
            name="dipole",
            func=lambda: self._get_dipole(),
            value=val,
            dependencies=[ensemble.beads._q],
        )

        pass

    def store(self, dipole):
        super().store(dipole)
        pass

    def _get_dipole(self, bead=None):
        """Return the electric dipole of all the beads as a list of np.array"""

        # check that 'bead' is correct
        if bead is not None:
            if bead < 0:
                raise ValueError("Error in '_get_dipole': 'beads' is negative.")
            if bead >= self.nbeads:
                raise ValueError(
                    "Error in '_get_dipole': 'beads' is greater than the number of beads."
                )

        dipole = np.full((self.nbeads, 3), np.nan)
        try:
            if "dipole" in self.ens.forces.extras:
                raws = [self.ens.forces.extras["dipole"][i] for i in range(self.nbeads)]
                for n, raw in enumerate(raws):
                    if len(raw) != 3:
                        raise ValueError("'dipole' has not length 3")
                    dipole[n] = np.asarray(raw)
                return dipole

            elif "raw" not in self.ens.forces.extras:
                raise ValueError("'raw' has to be in 'forces.extras'")

            elif np.all(
                ["Total dipole moment" in s for s in self.ens.forces.extras["raw"]]
            ):
                raws = [self.ens.forces.extras["raw"][i] for i in range(self.nbeads)]
                for n, raw in enumerate(raws):
                    factor = 1.0
                    if "[eAng]" in raw:
                        factor = UnitMap["length"]["angstrom"]

                    pattern = r"[-+]?\d*\.\d+(?:[eE][-+]?\d+)?|\b[-+]?\d+\b"
                    matches = list(re.findall(pattern, raw))
                    if len(matches) != 3:
                        raise ValueError(
                            "wrong number of extracted values from the extra string: they should be 3."
                        )
                    else:
                        dipole[n] = float(factor) * np.asarray(matches)
                    return dipole
            else:
                raise ValueError(
                    "Error in '_get_dipole': can not extract dipole from the extra string."
                )
        except:
            return np.full((self.nbeads, 3), np.nan)


dproperties(ElectricDipole, ["dipole", "nbeads", "forces"])


class ElectricField:
    """Class to handle the time dependent electric field when performing driven dynamics (with 'eda-nve')"""

    def __init__(self, amp=None, freq=None, phase=None, peak=None, sigma=None):
        self._amp = depend_array(
            name="amp", value=amp if amp is not None else np.zeros(3)
        )
        self._freq = depend_value(name="freq", value=freq if freq is not None else 0.0)
        self._phase = depend_value(
            name="phase", value=phase if phase is not None else 0.0
        )
        self._peak = depend_value(name="peak", value=peak if peak is not None else 0.0)
        self._sigma = depend_value(
            name="sigma", value=sigma if sigma is not None else np.inf
        )
        self.enabled = True

    def bind(self, driven_dyn, enstype: str):
        self.enstype = enstype
        self.enabled = enstype in driven_dyn.driven_integrators
        pass

    def store(self, ef):
        super(ElectricField, self).store(ef)
        self.amp.store(ef.amp)
        self.freq.store(ef.freq)
        self.phase.store(ef.phase)
        self.peak.store(ef.peak)
        self.sigma.store(ef.sigma)
        pass

    def Efield(self, time):
        """Get the value of the external electric field (cartesian axes)"""
        if not self.enabled:  # return a zero field with the correct shape
            if hasattr(time, "__len__"):
                return np.zeros((len(time), 3), float)
            else:
                return np.zeros(3, float)
        else:
            Eenv = self.Eenvelope(time)  # evaluate the envelope function
            Ecos = self._get_Ecos(time)  # evaluate the cos function
            if hasattr(time, "__len__"):
                return np.outer(Ecos * Eenv, self.amp)
            else:
                return Ecos * Eenv * self.amp

    def _Eenvelope_is_on(self):
        return self.peak > 0.0 and self.sigma != np.inf

    def Eenvelope(self, time):
        """Get the gaussian envelope function of the external electric field"""
        if self._Eenvelope_is_on():
            x = time  # indipendent variable
            u = self.peak  # mean value
            s = self.sigma  # standard deviation
            return np.exp(
                -0.5 * ((x - u) / s) ** 2
            )  # the returned maximum value is 1, when x = u
        else:
            return 1.0

    def _get_Ecos(self, time):
        """Get the sinusoidal part of the external electric field"""
        return np.cos(self.freq * time + self.phase)


dproperties(
    ElectricField,
    ["amp", "phase", "peak", "sigma", "freq"],
)

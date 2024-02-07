"""TODO
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

import numpy as np

from ipi.utils.depend import depend_value, dproperties
from ipi.utils.depend import dstrip
from ipi.utils.depend import *
from ipi.utils.units import UnitMap
import re

__all__ = ["Motion","BEC", "ElectricField", "EDA"]


class Motion:
    """Base motion calculation class.

    Gives the standard methods and attributes needed in all the
    motion calculation classes.

    Attributes:
        beads: A beads object giving the atoms positions.
        cell: A cell object giving the system box.
        forces: A forces object giving the virial and the forces acting on
            each bead.
        fixcom: A boolean which decides whether the centre of mass
            motion will be constrained or not.
        fixatoms: A list of atoms that should be held fixed to their
            initial positions.

    Depend objects:
        none
    """

    def __init__(self, fixcom=False, fixatoms=None):
        """Initialises Motion object.

        Args:
           fixcom: An optional boolean which decides whether the centre of mass
              motion will be constrained or not. Defaults to False.
           fixatoms: A list of atoms that should be held fixed to their
             initial positions.
        """

        self._dt = depend_value(name="dt", value=0.0)
        self.fixcom = fixcom
        if fixatoms is None:
            self.fixatoms = np.zeros(0, int)
        else:
            self.fixatoms = fixatoms

        self.beads = self.cell = self.forces = self.prng = self.nm = self.enstype = None

    def bind(self, ens, beads, nm, cell, bforce, prng, omaker):
        """Binds beads, cell, bforce, and prng to the calculator.

        This takes a beads object, a cell object, a forcefield object and a
        random number generator object and makes them members of the atom motion calculator.
        It also then creates the objects that will hold the data needed in the
        ensemble algorithms and the dependency network. Note that the conserved
        quantity is defined in the init, but as each ensemble has a different
        conserved quantity the dependencies are defined in bind.

        Args:
            beads: The beads object from which the bead positions are taken.
            nm: A normal modes object used to do the normal modes transformation.
            cell: The cell object from which the system box is taken.
            bforce: The forcefield object from which the force and virial are taken.
            prng: The random number generator object which controls random number
                generation.
        """

        # store local references to the different bits of the simulation
        self.beads = beads
        self.cell = cell
        self.forces = bforce
        self.prng = prng
        self.nm = nm
        self.ensemble = ens
        self.output_maker = omaker

    def step(self, step=None):
        """Dummy simulation time step which does nothing."""

        pass


dproperties(Motion, "dt")

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
    # The BEC tensors Z^* can be given to i-PI by an external driven through the etxra strings in the forces
    # of they can be kept fixed during the dynamics: in this case you can provide them through a txt file.

    def __init__(self, cbec=None, bec=None):
        self.cbec = cbec
        if bec is None:
            bec = np.full((0, 3), np.nan)
        self._bec = depend_array(name="bec", value=bec)
        pass

    def bind(self, eda, ensemble, enstype):
        self.enstype = enstype
        self.nbeads = ensemble.beads.nbeads
        self.natoms = ensemble.beads.natoms
        self.forces = ensemble.forces

        if self.enstype in EDA.integrators and self.cbec:
            self._bec = depend_array(
                name="bec",
                value=np.full((self.nbeads, 3 * self.natoms, 3), np.nan),
                func=self._get_driver_BEC,
                dependencies=[ensemble.beads._q],
            )
        elif self.enstype in EDA.integrators:
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

    def bind(self, eda, ensemble):
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
        self._Efield = depend_array(
            name="Efield",
            value=np.zeros(3, float),
            func=lambda: np.zeros(3, float),
        )
        self._Eenvelope = depend_value(
            name="Eenvelope",
            value=1.0,
            func=lambda: np.zeros(3, float),
        )

    def bind(self, eda, enstype):
        self.enstype = enstype
        self._mts_time = eda._mts_time

        dep = [self._mts_time, self._peak, self._sigma]
        self._Eenvelope = depend_value(
            name="Eenvelope", value=1.0, func=self._get_Eenvelope, dependencies=dep
        )

        if enstype in EDA.integrators:
            # dynamics is not driven --> add dependencies to the electric field
            dep = [self._mts_time, self._amp, self._freq, self._phase, self._Eenvelope]
            self._Efield = depend_array(
                name="Efield",
                value=np.zeros(3, float),
                func=self._get_Efield,
                dependencies=dep,
            )
        else:
            # dynamics is not driven --> no dependencies for the electric field
            self._Efield = depend_array(
                name="Efield",
                value=np.zeros(3, float),
                func=lambda: np.zeros(3, float),
            )
        pass

    def store(self, ef):
        super(ElectricField, self).store(ef)
        self.amp.store(ef.amp)
        self.freq.store(ef.freq)
        self.phase.store(ef.phase)
        self.peak.store(ef.peak)
        self.sigma.store(ef.sigma)
        pass

    def _get_Efield(self):
        """Get the value of the external electric field (cartesian axes)"""
        time = dstrip(self.mts_time)
        if hasattr(time, "__len__"):
            return np.outer(self._get_Ecos(time) * self.Eenvelope, self.amp)
        else:
            return self._get_Ecos(time) * self.Eenvelope * self.amp

    def _Eenvelope_is_on(self):
        return self.peak > 0.0 and self.sigma != np.inf

    def _get_Eenvelope(self):
        time = dstrip(self.mts_time)
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
    ["amp", "phase", "peak", "sigma", "freq", "Eenvelope", "Efield", "mts_time"],
)


class EDA:
    """Class to handle in a compact way 'BEC', 'ElectricDipole', and 'ElectricField' objects when performing driven dynamics (with 'eda-nve')"""

    integrators = ["eda-nve"]

    def __init__(self, efield: ElectricField, bec: BEC, **kwargv):
        super(EDA, self).__init__(**kwargv)
        self.Electric_Field = efield
        self.Electric_Dipole = ElectricDipole()
        self.Born_Charges = bec
        self._time = depend_value(name="time", value=0)
        self._mts_time = depend_value(name="mts_time", value=0)
        pass

    def bind(self, ensemble, enstype):
        self.enstype = enstype
        self._mts_time = depend_value(name="mts_time", value=dstrip(ensemble).time)
        self._time = ensemble._time
        self.Electric_Field.bind(self, enstype)
        self.Electric_Dipole.bind(self, ensemble)
        self.Born_Charges.bind(self, ensemble, enstype)
        pass

    def store(self, eda):
        super(EDA, self).store(eda)
        self.Electric_Field.store(eda.Electric_Field)
        self.Electric_Dipole.store(eda.Electric_Dipole)
        self.Born_Charges.store(eda.bec)
        pass


dproperties(EDA, ["mts_time", "time"])

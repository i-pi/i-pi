import numpy as np
from ipi.utils.depend import dd, dpipe, dstrip
from ipi.utils.depend import dobject, depend_array, depend_value
from ipi.utils.units import UnitMap
import re

__all__ = ["BEC", "ElectricField", "EDA"]


class BEC(dobject):
    def __init__(self, cbec=None, bec=None):
        self.cbec = cbec
        if bec is None:
            bec = np.full((0, 3), np.nan)
        dd(self).bec = depend_array(name="bec", value=bec)
        pass

    def bind(self, eda, ensemble, enstype):
        self.enstype = enstype
        self.nbeads = ensemble.beads.nbeads
        self.natoms = ensemble.beads.natoms
        self.forces = ensemble.forces
        # my_nbeads = 1  # self.nbeads

        dself = dd(self)
        if self.enstype in EDA.integrators and self.cbec:
            dself.bec = depend_array(
                name="bec",
                value=np.full(
                    (self.nbeads, 3 * self.natoms, 3), np.nan
                ),  # value=np.full((self.natoms,3,3),np.nan),\
                func=self._get_driver_BEC,
                dependencies=[dd(ensemble.beads).q],
            )
        elif self.enstype in EDA.integrators:
            temp = self._get_fixed_BEC()  # reshape the BEC once and for all
            dself.bec = depend_array(name="bec", value=temp)
        else:
            # dself.bec = depend_array(name="bec",value=np.full((self.natoms,3,3),np.nan))
            dself.bec = depend_array(
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


class ElectricDipole(dobject):
    # def __init__(self, cdip):
    #     self.cdip = cdip
    #     pass

    def bind(self, eda, ensemble):
        dself = dd(self)
        dself.nbeads = depend_value(name="nbeads", value=ensemble.beads.nbeads)

        self.forces = ensemble.forces  # is this a weakref??

        # val = np.zeros(
        #     3, dtype=float
        # )
        val = np.full(
            (self.nbeads, 3), np.nan
        )  # if self.nbeads > 1 else np.zeros(3,dtype=float)
        dself.dipole = depend_array(
            name="dipole",
            func=lambda: self._get_dipole(),
            value=val,
            dependencies=[dd(ensemble.beads).q],
        )

        pass

    def store(self, dipole):
        super().store(dipole)
        # self.cdip.store(dipole.cdip)
        pass

    def _get_dipole(self, bead=None):
        """Return the electric dipole of all the beads as a list of np.array"""

        # check that bead is a correct value
        if bead is not None:
            if bead < 0:
                raise ValueError("Error in '_get_dipole': 'beads' is negative.")
            if bead >= self.nbeads:
                raise ValueError(
                    "Error in '_get_dipole': 'beads' is greater than the number of beads."
                )

        dipole = np.full((self.nbeads, 3), np.nan)
        try:
            if "dipole" in self.forces.extras:
                raws = [self.forces.extras["dipole"][i] for i in range(self.nbeads)]
                for n, raw in enumerate(raws):
                    if len(raw) != 3:
                        raise ValueError("'dipole' has not length 3")
                    dipole[n] = np.asarray(raw)
                return dipole

            elif "raw" not in self.forces.extras:
                raise ValueError("'raw' has to be in 'forces.extras'")

            elif np.all(
                ["Total dipole moment" in s for s in self.forces.extras["raw"]]
            ):
                raws = [self.forces.extras["raw"][i] for i in range(self.nbeads)]
                for n, raw in enumerate(raws):
                    factor = 1.0
                    if "[eAng]" in raw:
                        factor = UnitMap["length"]["angstrom"]

                    pattern = "[-+]?\d*\.\d+(?:[eE][-+]?\d+)?|\b[-+]?\d+\b"
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
            # if bead is not None:
            #     return np.full(3,np.nan)
            # else :
            return np.full((self.nbeads, 3), np.nan)


class ElectricField(dobject):
    def __init__(self, amp=None, freq=None, phase=None, peak=None, sigma=None):
        dself = dd(self)
        dself.amp = depend_array(
            name="amp", value=amp if amp is not None else np.zeros(3)
        )
        dself.freq = depend_value(name="freq", value=freq if freq is not None else 0.0)
        dself.phase = depend_value(
            name="phase", value=phase if phase is not None else 0.0
        )
        dself.peak = depend_value(name="peak", value=peak if peak is not None else 0.0)
        dself.sigma = depend_value(
            name="sigma", value=sigma if sigma is not None else np.inf
        )

    def bind(self, eda, enstype):
        self.enstype = enstype
        dself = dd(self)
        dself.mts_time = dd(eda).mts_time

        # same dependencies for Eenvelope and its time derivative
        dep = [dself.mts_time, dself.peak, dself.sigma]
        dself.Eenvelope = depend_value(
            name="Eenvelope", value=1.0, func=self._get_Eenvelope, dependencies=dep
        )

        if enstype in EDA.integrators:
            # with dependencies
            dep = [dself.mts_time, dself.amp, dself.freq, dself.phase, dself.Eenvelope]
            dself.Efield = depend_array(
                name="Efield",
                value=np.zeros(3, float),
                func=self._get_Efield,
                dependencies=dep,
            )
        else:
            # no dependencies
            dself.Efield = depend_array(
                name="Efield",
                value=np.zeros(3, float),
                func=lambda time=None: np.zeros(3, float),
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
        # https://en.wikipedia.org/wiki/Normal_distribution
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
        # it's easier to define a function and compute this 'cos'
        # again everytime instead of define a 'depend_value'
        return np.cos(self.freq * time + self.phase)


class EDA(dobject):
    integrators = ["eda-nve"]

    def __init__(self, efield: ElectricField, bec: BEC, **kwargv):
        super(EDA, self).__init__(**kwargv)
        self.Electric_Field = efield
        self.Electric_Dipole = ElectricDipole()
        self.Born_Charges = bec
        pass

    def bind(self, ensemble, enstype):
        self.enstype = enstype
        dself = dd(self)

        dself.mts_time = depend_value(name="mts_time", value=0.0)
        dself.time = dd(ensemble).time

        self.Electric_Field.bind(self, enstype)
        self.Electric_Dipole.bind(self, ensemble)
        self.Born_Charges.bind(self, ensemble, enstype)

        dself.Efield = dd(self.Electric_Field).Efield
        dself.Eenvelope = dd(self.Electric_Field).Eenvelope
        dself.bec = dd(self.Born_Charges).bec
        dself.dipole = dd(self.Electric_Dipole).dipole

        pass

    def store(self, eda):
        super(EDA, self).store(eda)
        self.Electric_Field.store(eda.Electric_Field)
        self.Electric_Dipole.store(eda.Electric_Dipole)
        self.Born_Charges.store(eda.bec)
        pass

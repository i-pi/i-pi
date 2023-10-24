import numpy as np
from ipi.utils.depend import dd, dpipe
from ipi.utils.depend import dobject, depend_array, depend_value
from ipi.utils.units import UnitMap
import re


class BEC(dobject):
    def __init__(self, cbec, bec):
        self.cbec = cbec
        dd(self).bec = depend_array(name="bec", value=bec)
        pass

    def bind(self, eda, ensemble, enstype):
        self.enstype = enstype
        self.nbeads = ensemble.beads.nbeads
        self.natoms = ensemble.beads.natoms
        self.forces = ensemble.forces
        my_nbeads = 1  # self.nbeads

        dself = dd(self)
        if self.cbec:
            dself.bec = depend_array(
                name="bec",
                value=np.full(
                    (my_nbeads, 3 * self.natoms, 3), np.nan
                ),  # value=np.full((self.natoms,3,3),np.nan),\
                func=self._get_otf_BEC,
                dependencies=[dd(eda).time, dd(ensemble.beads).q],
            )
        elif self.enstype in EDA.integrators:
            temp = self._get_static_BEC()  # reshape the BEC once and for all
            dself.bec = depend_array(name="bec", value=temp)
        else:
            # dself.bec = depend_array(name="bec",value=np.full((self.natoms,3,3),np.nan))
            dself.bec = depend_array(
                name="bec", value=np.full((my_nbeads, 3 * self.natoms, 3), np.nan)
            )

        pass

    def store(self, bec):
        super(BEC, self).store(bec)
        self.cbec.store(bec.cbec)

    def _get_otf_BEC(self, bead=None):
        """Return the BEC tensors (in cartesian coordinates), when computed on the fly (otf) by the driver"""

        msg = "Error in '_get_otf_BEC'"

        if bead is not None:
            if bead < 0:
                raise ValueError("Error in '_get_otf_BEC': 'bead' is negative")
            if bead >= self.nbeads:
                raise ValueError(
                    "Error in '_get_otf_BEC': 'bead' is greater than the number of beads"
                )
        else:
            if self.nbeads != 1:
                raise ValueError(
                    "Error in '_get_otf_BEC': EDA integration has not implemented yet for 'nbeads' > 1"
                )

        if self.cbec:
            if "BEC" not in self.forces.extras:
                raise ValueError(
                    msg
                    + ": BEC tensors are not returned to i-PI (or at least not accessible in '_get_otf_BEC')."
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

    def _get_static_BEC(self):
        """Return the BEC tensors (in cartesian coordinates).
        The BEC tensor are stored in a compact form.
        This method trasform the BEC tensors into another data structure, suitable for computation.
        A lambda function is also returned to perform fast matrix multiplication.
        """

        return self.bec.reshape((self.nbeads, 3 * self.natoms, 3))


class Dipole(dobject):
    # def __init__(self, cdip):
    #     self.cdip = cdip
    #     pass

    def bind(self, eda, ensemble):
        dself = dd(self)
        dself.nbeads = depend_value(name="nbeads", value=ensemble.beads.nbeads)

        self.forces = ensemble.forces  # is this a weakref??

        val = np.zeros(
            3, dtype=float
        )  # np.full(self.nbeads,np.zeros(3,dtype=float)) if self.nbeads > 1 else np.zeros(3,dtype=float)
        dself._dipole_ = depend_array(
            name="_dipole_",
            func=lambda: self._get_dipole(bead=0),
            value=val,
            dependencies=[dd(eda).time, dd(ensemble.beads).q],
        )

        pass

    def store(self, dipole):
        super(Dipole, self).store(dipole)
        # self.cdip.store(dipole.cdip)
        pass

    def _get_dipole(self, bead=None):
        """Return the electric dipole of all the beads as a list of np.array"""
        # self._check_dipole()

        # check that bead is a correct value
        # N = self.beads.nbeads
        if bead is not None:
            if bead < 0:
                raise ValueError("Error in '_get_dipole': 'beads' is negative.")
            if bead >= self.nbeads:
                raise ValueError(
                    "Error in '_get_dipole': 'beads' is greater than the number of beads."
                )
            if bead > 1:
                raise ValueError(
                    "The case with 'beads' != 0 has not been implemeted yet"
                )

        # if not self.cdip:
        #     return np.asarray([0, 0, 0])
        # else:
        try :
            if "dipole" in self.forces.extras:
                dipole = np.asarray(self.forces.extras["dipole"]).flatten()
                if len(dipole) != 3:
                    print("dipole:", dipole)
                    raise ValueError("'dipole' has not length 3")
                return dipole
                # dipole = [ self.forces.extras["dipole"][i] for i in range(self.nbeads)]
                # return dipole[0] if bead is None else dipole[bead]

            elif "raw" not in self.forces.extras:
                raise ValueError("'raw' has to be in 'forces.extras'")

            elif np.all(
                ["Total dipole moment" in s for s in self.forces.extras["raw"]]
            ):
                raw = [self.forces.extras["raw"][i] for i in range(self.nbeads)]
                raw = raw[0] if bead is None else raw[bead]
                factor = 1.0
                if "[eAng]" in raw:
                    factor = UnitMap["length"]["angstrom"]

                dipole = np.full(3, np.nan)
                pattern = "[-+]?\d*\.\d+(?:[eE][-+]?\d+)?|\b[-+]?\d+\b"
                matches = re.findall(pattern, raw)
                if len(matches) != 3:
                    raise ValueError(
                        "wrong number of extracted values from the extra string: they should be 3."
                    )
                else:
                    for n in range(3):
                        dipole[n] = float(matches[n])
                return float(factor) * dipole
            else:
                raise ValueError(
                    "Error in '_get_dipole': can not extract dipole from the extra string."
                )
        except:
            if bead is not None:
                return np.full(3,np.nan)
            else :
                return np.full((self.nbeads,3),np.nan) 


class ElectricField(dobject):
    def __init__(self, Eamp, Efreq, Ephase, Epeak, Esigma):
        dself = dd(self)
        dself.Eamp = depend_array(
            name="Eamp", value=Eamp if Eamp is not None else np.zeros(3)
        )
        dself.Efreq = depend_value(
            name="Efreq", value=Efreq if Efreq is not None else 0.0
        )
        dself.Ephase = depend_value(
            name="Ephase", value=Ephase if Ephase is not None else 0.0
        )
        dself.Epeak = depend_value(
            name="Epeak", value=Epeak if Epeak is not None else 0.0
        )
        dself.Esigma = depend_value(
            name="Esigma", value=Esigma if Esigma is not None else np.inf
        )

        # these will be overwritten
        dself.Eenvelope = depend_value(
            name="Eenvelope", value=1.0, func=self._get_Eenvelope
        )
        dself.Efield = depend_array(
            name="Efield", value=np.zeros(3, float), func=self._get_Efield
        )

    def bind(self, eda, enstype):
        self.enstype = enstype
        dself = dd(self)
        dself.cptime = depend_value(name="cptime", value=0.0)
        dpipe(dfrom=dd(eda).cptime, dto=dd(self).cptime)

        # same dependencies for Eenvelope and its time derivative
        dep = [dself.cptime, dself.Epeak, dself.Esigma]
        dself.Eenvelope = depend_value(
            name="Eenvelope", value=1.0, func=self._get_Eenvelope, dependencies=dep
        )

        if enstype in EDA.integrators:
            # with dependencies
            dep = [dself.cptime, dself.Eamp, dself.Efreq, dself.Ephase, dself.Eenvelope]
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
        self.Eamp.store(ef.Eamp)
        self.Efreq.store(ef.Efreq)
        self.Ephase.store(ef.Ephase)
        self.Epeak.store(ef.Epeak)
        self.Esigma.store(ef.Esigma)
        pass

    def _get_Efield(self, time=None):
        """Get the value of the external electric field (cartesian axes)"""
        if time is None:
            raise ValueError(
                "Hey man! Don't you think it's better to specify the time you want to evaluate the electric field?"
            )
        if hasattr(time, "__len__"):
            return np.outer(self._get_Ecos(time) * dd(self).Eenvelope(time), self.Eamp)
        else:
            return self._get_Ecos(time) * dd(self).Eenvelope(time) * self.Eamp

    def _Eenvelope_is_on(self):
        return self.Epeak > 0.0 and self.Esigma != np.inf

    def _get_Eenvelope(self, time=None):
        """Get the gaussian envelope function of the external electric field"""
        # https://en.wikipedia.org/wiki/Normal_distribution
        if self._Eenvelope_is_on():
            x = self.cptime if time is None else time  # indipendent variable
            u = self.Epeak  # mean value
            s = self.Esigma  # standard deviation
            return np.exp(
                -0.5 * ((x - u) / s) ** 2
            )  # the returned maximum value is 1, when x = u
        else:
            return 1.0

    def _get_Ecos(self, time=None):
        """Get the sinusoidal part of the external electric field"""
        t = self.cptime if time is None else time
        return np.cos(self.Efreq * t + self.Ephase)


class EDA(dobject):
    integrators = ["eda-nve"]

    def __init__(self, Eamp, Efreq, Ephase, Epeak, Esigma, cbec, bec, **kwargv): # cdip
        super(EDA, self).__init__(**kwargv)
        self.Electric_Field = ElectricField(Eamp, Efreq, Ephase, Epeak, Esigma)
        self.Dipole = Dipole() # (cdip)
        self.Born_Charges = BEC(cbec, bec)
        pass

    def bind(self, ensemble, enstype):
        self.enstype = enstype
        dself = dd(self)

        dself.econs = depend_value(name="econs", value=0.0)
        dself.time = depend_value(name="time", value=0.0)
        dself.cptime = depend_value(name="cptime", value=0.0)

        dpipe(dfrom=dd(ensemble).econs, dto=dself.econs)
        dpipe(dfrom=dd(ensemble).time, dto=dself.time)

        self.Electric_Field.bind(self, enstype)
        self.Dipole.bind(self, ensemble)
        self.Born_Charges.bind(self, ensemble, enstype)

        # for easier access
        dself.Efield = depend_array(
            name="Efield",
            value=np.full(dd(self.Electric_Field).Efield.shape, np.nan),
            func=lambda time=None: dd(self.Electric_Field).Efield(time),
        )
        dself.bec = depend_array(
            name="bec", value=np.full(dd(self.Born_Charges).bec.shape, np.nan)
        )
        dself.dipole = depend_array(
            name="dipole", value=np.full(dd(self.Dipole)._dipole_.shape, np.nan)
        )

        dpipe(dfrom=dd(self.Born_Charges).bec, dto=dself.bec)
        dpipe(dfrom=dd(self.Dipole)._dipole_, dto=dself.dipole)
        dpipe(dfrom=dself.time, dto=dself.cptime)

        pass

    def store(self, eda):
        super(EDA, self).store(eda)
        self.Electric_Field.store(eda.Electric_Field)
        self.Dipole.store(eda.Dipole)
        self.Born_Charges.store(eda.bec)
        pass

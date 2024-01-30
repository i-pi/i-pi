"""Creates objects that deal with constant temperature simulations.

Chooses between the different possible thermostat options and creates the
appropriate thermostat object, with suitable parameters.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


from copy import copy

import numpy as np

import ipi.engine.thermostats as ethermostats
from ipi.utils.depend import *
from ipi.utils.inputvalue import *


__all__ = ["InputThermo"]


class InputThermoBase(Input):
    """Thermostat input class.

    Handles generating the appropriate thermostat class from the xml input file,
    and generating the xml checkpoint tags and data from an instance of the
    object.

    Attributes:
       mode: An optional string giving the type of the thermostat used. Defaults
          to 'langevin'.

    Fields:
       ethermo: An optional float giving the amount of heat energy transferred
          to the bath. Defaults to 0.0.
       tau: An optional float giving the damping time scale. Defaults to 1.0.
       pile_lambda: Scaling for the PILE damping relative to the critical damping.
       pile_centroid_t: Option to set a different temperature to the centroid in PILE thermostats.
       A: An optional array of floats giving the drift matrix. Defaults to 0.0.
       C: An optional array of floats giving the static covariance matrix.
          Defaults to 0.0.
       s: An optional array of floats giving the additional momentum-scaled
          momenta in GLE. Defaults to 0.0.
       intau: Sets the estimated time scale of the inherent force noise term. Defaults to 0.0.
       idtau: Sets the estimated time scale of the forces' inherent dissipation. Defaults to 0.0.
       apat: Time scale for automatic parameter adjustment of intau or idtau. Defaults to 0.0.
       flip: Flipping type for FFL thermostat ('soft', 'hard', 'rescale', 'none'). Defaults to 'rescale'

    """

    attribs = {
        "mode": (
            InputAttribute,
            {
                "dtype": str,
                "options": [
                    "",
                    "langevin",
                    "svr",
                    "pile_l",
                    "pile_g",
                    "gle",
                    "nm_gle",
                    "nm_gle_g",
                    "cl",
                    "ffl",
                ],
                "help": "The style of thermostatting. 'langevin' specifies a white noise langevin equation to be attached to the cartesian representation of the momenta. 'svr' attaches a velocity rescaling thermostat to the cartesian representation of the momenta. Both 'pile_l' and 'pile_g' attaches a white noise langevin thermostat to the normal mode representation, with 'pile_l' attaching a local langevin thermostat to the centroid mode and 'pile_g' instead attaching a global velocity rescaling thermostat. 'gle' attaches a coloured noise langevin thermostat to the cartesian representation of the momenta, 'nm_gle' attaches a coloured noise langevin thermostat to the normal mode representation of the momenta and a langevin thermostat to the centroid and 'nm_gle_g' attaches a gle thermostat to the normal modes and a svr thermostat to the centroid. 'cl' represents a modified langevin thermostat which compensates for additional white noise from noisy forces or for dissipative effects. 'ffl' is the fast-forward langevin thermostat, in which momenta are flipped back whenever the action of the thermostat changes its direction. 'multiple' is a special thermostat mode, in which one can define multiple thermostats _inside_ the thermostat tag.",
            },
        )
    }
    fields = {
        "ethermo": (
            InputValue,
            {
                "dtype": float,
                "default": 0.0,
                "help": "The initial value of the thermostat energy. Used when the simulation is restarted to guarantee continuity of the conserved quantity.",
                "dimension": "energy",
            },
        ),
        "tau": (
            InputValue,
            {
                "dtype": float,
                "default": 0.0,
                "help": "The friction coefficient for white noise thermostats.",
                "dimension": "time",
            },
        ),
        "pile_lambda": (
            InputValue,
            {
                "dtype": float,
                "default": 1.0,
                "help": "Scaling for the PILE damping relative to the critical damping. (gamma_k=2*lambda*omega_k",
            },
        ),
        "pile_centroid_t": (
            InputValue,
            {
                "dtype": float,
                "default": 0.0,
                "help": "Option to set a different centroid temperature wrt. that of the ensemble. Only used if value other than 0.0.",
                "dimension": "temperature",
            },
        ),
        "A": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "The friction matrix for GLE thermostats.",
                "dimension": "frequency",
            },
        ),
        "C": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "The covariance matrix for GLE thermostats.",
                "dimension": "temperature",
            },
        ),
        "s": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "Input values for the additional momenta in GLE.",
                "dimension": "ms-momentum",
            },
        ),
        "intau": (
            InputValue,
            {
                "dtype": float,
                "default": 0.0,
                "help": "The inherent noise time scale for compensating langevin thermostats.",
                "dimension": "time",
            },
        ),
        "idtau": (
            InputValue,
            {
                "dtype": float,
                "default": 0.0,
                "help": "The inherent dissipation time scale for compensating langevin thermostats.",
                "dimension": "time",
            },
        ),
        "apat": (
            InputValue,
            {
                "dtype": float,
                "default": 0.0,
                "help": "The time scale for automatic adjustment of CL thermostat's parameters.",
                "dimension": "time",
            },
        ),
        "flip": (
            InputValue,
            {
                "dtype": str,
                "default": "rescale",
                "help": "Flipping type for ffl thermostat ('soft', 'hard', 'rescale', 'none')",
            },
        ),
    }

    dynamic = {}

    default_help = "Simulates an external heat bath to keep the velocity distribution at the correct temperature."
    default_label = "THERMOSTATS"

    def store(self, thermo):
        """Takes a thermostat instance and stores a minimal representation of it.

        Args:
           thermo: A thermostat object.

        Raises:
           TypeError: Raised if the thermostat is not a recognized type.
        """

        super(InputThermoBase, self).store(thermo)
        if type(thermo) is ethermostats.ThermoLangevin:
            self.mode.store("langevin")
            self.tau.store(thermo.tau)
        elif type(thermo) is ethermostats.ThermoSVR:
            self.mode.store("svr")
            self.tau.store(thermo.tau)
        elif type(thermo) is ethermostats.ThermoPILE_L:
            self.mode.store("pile_l")
            self.tau.store(thermo.tau)
            self.pile_lambda.store(thermo.pilescale)
            self.pile_centroid_t.store(thermo.pilect)
        elif type(thermo) is ethermostats.ThermoPILE_G:
            self.mode.store("pile_g")
            self.tau.store(thermo.tau)
            self.pile_lambda.store(thermo.pilescale)
            self.pile_centroid_t.store(thermo.pilect)
        elif type(thermo) is ethermostats.ThermoGLE:
            self.mode.store("gle")
            self.A.store(thermo.A)
            if thermo._C._func is None:
                self.C.store(thermo.C)
            self.s.store(thermo.s)
        elif type(thermo) is ethermostats.ThermoNMGLE:
            self.mode.store("nm_gle")
            self.A.store(thermo.A)
            if thermo._C._func is None:
                self.C.store(thermo.C)
            self.s.store(thermo.s)
        elif type(thermo) is ethermostats.ThermoNMGLEG:
            self.mode.store("nm_gle_g")
            self.A.store(thermo.A)
            self.tau.store(thermo.tau)
            if thermo._C._func is None:
                self.C.store(thermo.C)
            self.s.store(thermo.s)
        elif type(thermo) is ethermostats.ThermoCL:
            self.mode.store("cl")
            self.tau.store(thermo.tau)
            self.intau.store(thermo.intau)
            self.idtau.store(thermo.idtau)
            self.apat.store(thermo.apat)
        elif type(thermo) is ethermostats.ThermoFFL:
            self.mode.store("ffl")
            self.tau.store(thermo.tau)
            self.flip.store(thermo.flip)
        elif type(thermo) is ethermostats.Thermostat:
            self.mode.store("")
        else:
            raise TypeError("Unknown thermostat mode " + type(thermo).__name__)
        self.ethermo.store(thermo.ethermo)

    def fetch(self):
        """Creates a thermostat object.

        Returns:
           A thermostat object of the appropriate type and with the appropriate
           parameters given the attributes of the InputThermo object.

        Raises:
           TypeError: Raised if the thermostat type is not a recognized option.
        """

        super(InputThermoBase, self).fetch()
        if self.mode.fetch() == "langevin":
            thermo = ethermostats.ThermoLangevin(tau=self.tau.fetch())
        elif self.mode.fetch() == "svr":
            thermo = ethermostats.ThermoSVR(tau=self.tau.fetch())
        elif self.mode.fetch() == "pile_l":
            thermo = ethermostats.ThermoPILE_L(
                tau=self.tau.fetch(),
                scale=self.pile_lambda.fetch(),
                pilect=self.pile_centroid_t.fetch(),
            )
        elif self.mode.fetch() == "pile_g":
            thermo = ethermostats.ThermoPILE_G(
                tau=self.tau.fetch(),
                scale=self.pile_lambda.fetch(),
                pilect=self.pile_centroid_t.fetch(),
            )
        elif self.mode.fetch() == "gle":
            rC = self.C.fetch()
            if len(rC) == 0:
                rC = None
            thermo = ethermostats.ThermoGLE(A=self.A.fetch(), C=rC)
            thermo.s = self.s.fetch()
        elif self.mode.fetch() == "nm_gle":
            rC = self.C.fetch()
            if len(rC) == 0:
                rC = None
            thermo = ethermostats.ThermoNMGLE(A=self.A.fetch(), C=rC)
            thermo.s = self.s.fetch()
        elif self.mode.fetch() == "nm_gle_g":
            rC = self.C.fetch()
            if len(rC) == 0:
                rC = None
            thermo = ethermostats.ThermoNMGLEG(
                A=self.A.fetch(), C=rC, tau=self.tau.fetch()
            )
            thermo.s = self.s.fetch()
        elif self.mode.fetch() == "cl":
            thermo = ethermostats.ThermoCL(
                tau=self.tau.fetch(),
                intau=self.intau.fetch(),
                idtau=self.idtau.fetch(),
                apat=self.apat.fetch(),
            )
        elif self.mode.fetch() == "ffl":
            thermo = ethermostats.ThermoFFL(
                tau=self.tau.fetch(), flip=self.flip.fetch()
            )
        elif self.mode.fetch() == "":
            thermo = ethermostats.Thermostat()
        else:
            raise TypeError("Invalid thermostat mode " + self.mode.fetch())

        thermo.ethermo = self.ethermo.fetch()

        return thermo

    def check(self):
        """Checks that the parameter arrays represents a valid thermostat."""

        super(InputThermoBase, self).check()
        mode = self.mode.fetch()

        if mode in ["langevin", "svr", "pile_l", "pile_g", "nm_gle_g", "ffl"]:
            if self.tau.fetch() <= 0:
                raise ValueError(
                    "The thermostat friction coefficient must be set to a positive value"
                )
        if mode == "cl":
            if self.tau.fetch() < 0:
                raise ValueError(
                    "The thermostat friction coefficient must be set to a non-negative value"
                )
            if self.intau.fetch() < 0:
                raise ValueError(
                    "The inherent noise time scale must be set to a non-negative value"
                )
            if self.idtau.fetch() < 0:
                raise ValueError(
                    "The inherent dissipation time scale must be set to a non-negative value"
                )
            if self.apat.fetch() < 0:
                raise ValueError(
                    "The automatic parameter adjustment time scale must be set to a non-negative value"
                )
        if mode in ["gle", "nm_gle", "nm_gle_g"]:
            pass  # PERHAPS DO CHECKS THAT MATRICES SATISFY REASONABLE CONDITIONS (POSITIVE-DEFINITENESS, ETC)
        # MR Check that pilect is not less than 0.0


class InputThermo(InputThermoBase):
    """Extends InputThermoBase to allow the definition of a multithermo"""

    attribs = copy(InputThermoBase.attribs)

    attribs["mode"][1]["options"].append("multi")

    dynamic = {
        "thermostat": (
            InputThermoBase,
            {
                "default": input_default(factory=ethermostats.Thermostat),
                "help": "The thermostat for the atoms, keeps the atom velocity distribution at the correct temperature.",
            },
        )
    }

    def store(self, thermo):
        if type(thermo) is ethermostats.MultiThermo:
            self.mode.store("multi")

            if len(self.extra) != len(thermo.tlist):
                self.extra = [0] * len(thermo.tlist)
            for ii, t in enumerate(thermo.tlist):
                if self.extra[ii] == 0:
                    it = InputThermoBase()
                    it.store(t)
                    self.extra[ii] = ("thermostat", it)
                else:
                    self.extra[ii][1].store(t)

            self.ethermo.store(thermo.ethermo)
        else:
            super(InputThermo, self).store(thermo)

    def fetch(self):
        if self.mode.fetch() == "multi":
            tlist = []
            for k, t in self.extra:
                tlist.append(t.fetch())
            thermo = ethermostats.MultiThermo(thermolist=tlist)
            thermo.ethermo = self.ethermo.fetch()
        else:
            thermo = super(InputThermo, self).fetch()

        return thermo

"""Creates objects that deal with constant pressure and stress simulations."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np

import ipi.engine.thermostats
from ipi.engine.barostats import *
from ipi.utils.inputvalue import *
from ipi.inputs.thermostats import *
from ipi.engine.cell import Cell
from ipi.inputs.cell import *


__all__ = ["InputBaro"]


class InputBaro(Input):

    """Barostat input class.

    Handles generating the appropriate barostat class from the xml input file,
    and generating the xml checkpoint tags and data from an
    instance of the object.

    Attributes:
       mode: An optional string giving the type of barostat used. Defaults to
          'rigid'.

    Fields:
       thermostat: A thermostat object giving the cell thermostat.
       tau: The time constant associated with the dynamics of the piston.
       p: The conjugate momentum to the volume degree of freedom.
    """

    attribs = {
        "mode": (
            InputAttribute,
            {
                "dtype": str,
                "default": "dummy",
                "help": """The type of barostat.  Currently, only a 'isotropic' barostat is implemented, that combines
                                    ideas from the Bussi-Zykova-Parrinello barostat for classical MD with ideas from the
                                    Martyna-Hughes-Tuckerman centroid barostat for PIMD; see Ceriotti, More, Manolopoulos, Comp. Phys. Comm. 2013 for
                                    implementation details.""",
                "options": [
                    "dummy",
                    "isotropic",
                    "flexible",
                    "anisotropic",
                    "sc-isotropic",
                ],
            },
        )
    }
    fields = {
        "thermostat": (
            InputThermo,
            {
                "default": input_default(factory=ipi.engine.thermostats.Thermostat),
                "help": "The thermostat for the cell. Keeps the cell velocity distribution at the correct temperature. Note that the 'pile_l', 'pile_g', 'nm_gle' and 'nm_gle_g' options will not work for this thermostat.",
            },
        ),
        "tau": (
            InputValue,
            {
                "default": 1.0,
                "dtype": float,
                "dimension": "time",
                "help": "The time constant associated with the dynamics of the piston.",
            },
        ),
        "p": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "Momentum (or momenta) of the piston.",
                "dimension": "momentum",
            },
        ),
        "h0": (
            InputCell,
            {
                "dtype": float,
                "default": input_default(factory=Cell),
                "help": "Reference cell for Parrinello-Rahman-like barostats.",
                "dimension": "length",
            },
        ),
        "hfix": (
            InputArray,
            {
                "default": input_default(factory=np.zeros, args=(0, str)),
                "dtype": str,
                "help": "A list of the cell entries that should be held fixed (xx, yy, zz, xy, xz, yz). 'offdiagonal' is an alias for xy, xz, yz.",
            },
        ),
    }

    default_help = "Simulates an external pressure bath."
    default_label = "BAROSTAT"

    def store(self, baro):
        """Takes a barostat instance and stores a minimal representation of it.

        Args:
           baro: A barostat object.
        """

        super(InputBaro, self).store(baro)
        self.thermostat.store(baro.thermostat)
        self.tau.store(baro.tau)
        if type(baro) is BaroBZP:
            self.mode.store("isotropic")
            self.p.store(baro.p)
        elif type(baro) is BaroSCBZP:
            self.mode.store("sc-isotropic")
            self.p.store(baro.p)
        elif type(baro) is BaroMTK:
            self.mode.store("flexible")
            self.p.store(baro.p)
            self.hfix.store(baro.hfix)
        elif type(baro) is BaroRGB:
            self.mode.store("anisotropic")
            self.p.store(baro.p)
            self.h0.store(baro.h0)
            self.hfix.store(baro.hfix)
        elif type(baro) is Barostat:
            self.mode.store("dummy")
        else:
            raise TypeError(
                "The type " + type(baro).__name__ + " is not a valid barostat type"
            )

    def fetch(self):
        """Creates a barostat object.

        Returns:
           A barostat object of the appropriate type and with the appropriate
           thermostat given the attributes of the InputBaro object.
        """

        super(InputBaro, self).fetch()
        if self.mode.fetch() == "isotropic":
            baro = BaroBZP(thermostat=self.thermostat.fetch(), tau=self.tau.fetch())
            if self.p._explicit:
                baro.p = self.p.fetch()
        elif self.mode.fetch() == "sc-isotropic":
            baro = BaroSCBZP(thermostat=self.thermostat.fetch(), tau=self.tau.fetch())
            if self.p._explicit:
                baro.p = self.p.fetch()
        elif self.mode.fetch() == "flexible":
            baro = BaroMTK(
                thermostat=self.thermostat.fetch(),
                tau=self.tau.fetch(),
                hfix=self.hfix.fetch(),
            )
            if self.p._explicit:
                baro.p = self.p.fetch()
        elif self.mode.fetch() == "anisotropic":
            baro = BaroRGB(
                thermostat=self.thermostat.fetch(),
                tau=self.tau.fetch(),
                hfix=self.hfix.fetch(),
            )
            if self.p._explicit:
                baro.p = self.p.fetch()
            if self.h0._explicit:
                baro.h0 = self.h0.fetch()
            else:
                raise ValueError(
                    "Reference cell MUST be specified for an anisotropic barostat"
                )
        elif self.mode.fetch() == "dummy":
            baro = Barostat(thermostat=self.thermostat.fetch(), tau=self.tau.fetch())
        else:
            raise ValueError(self.mode.fetch() + " is not a valid mode of barostat")

        return baro

"""Harmonic potential"""

import sys

try:
    from .dummy import Dummy_driver
except:
    from dummy import Dummy_driver

import numpy as np
from ipi.utils import units


__DRIVER_NAME__ = "DW"
__DRIVER_CLASS__ = "DoubleWell_driver"


invcm2au = units.unit_to_internal("frequency", "inversecm", 1.0)
A2au = units.unit_to_internal("length", "angstrom", 1.0)

# ------------DOUBLE WELL POTENTIAL-----------------------------
#
#                                                 m^2*w_b^4
# V(x,w_b,v0) =  - 0.5 *m*w_b^2*(x-delta)^2 +    ---------- x^4
#                                                  16V0
#
# ---------------------------------------------------------


class DoubleWell_driver(Dummy_driver):
    """
    DW driver accepts 0 or 4 arguments.\nExample: python driver.py -m DoubleWell -o omega_b (cm^-1) V0 (cm^-1) mass(a.u) delta(angs) \n
        python driver.py -m DoubleWell -o 500,2085,1837,0.00 \n"""

    def __init__(self, w_b=None, v0=None, m=None, delta=None, *args, **kwargs):

        if w_b == None or v0 == None or m == None or delta == None:
            print("using default values from Craig-JCP-2005")
            # We used Craig's values (J. Chem. Phys. 122, 084106, 2005)
            w_b = 500 * invcm2au  # Tc = 115K
            v0 = 2085 * invcm2au
            m = 1837.36223469
            self.delta = 00
        else:
            try:
                w_b = w_b * invcm2au
                v0 = v0 * invcm2au
                self.delta = delta * A2au
            except:
                sys.exit(self.__doc__)

        self.k = 1837.36223469 * (3800.0 / 219323.0) ** 2
        self.A = -0.5 * m * (w_b) ** 2
        self.B = ((m**2) * (w_b) ** 4) / (16 * v0)
        super().__init__(*args, **kwargs)

    def __call__(self, cell, pos):
        """DoubleWell potential l"""
        pot = 0
        pos3 = pos.reshape(-1, 3)
        force3 = np.zeros(pos.shape)

        # DW
        pot += self.A * (pos3[:, 0] - self.delta) ** 2 + self.B * (pos3[:, 0] ** 4)
        force3[:, 0] = -2.0 * self.A * (pos3[:, 0] - self.delta) - 4.0 * self.B * (
            pos3[:, 0] ** 3
        )

        # Harmonic
        pot += 0.5 * self.k * (pos3[:, 1] ** 2).sum()
        pot += 0.5 * self.k * (pos3[:, 2] ** 2).sum()
        force3[:, 1] = -self.k * pos3[:, 1]
        force3[:, 2] = -self.k * pos3[:, 2]

        vir = cell * 0.0  # makes a zero virial with same shape as cell
        extras = "empty"
        pos = pos3.reshape(pos.shape)
        force = force3.reshape(pos.shape)

        return pot, force, vir, extras

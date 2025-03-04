"""Harmonic potential"""

import sys

try:
    from .dummy import Dummy_driver
    from .bath import Harmonic_Bath_explicit
except:
    from dummy import Dummy_driver
    from .bath import Harmonic_Bath_explicit

import numpy as np
from ipi.utils import units


__DRIVER_NAME__ = "DW_bath"
__DRIVER_CLASS__ = "DoubleWell_with_explicit_bath_driver"


invcm2au = units.unit_to_internal("frequency", "inversecm", 1.0)
A2au = units.unit_to_internal("length", "angstrom", 1.0)

# ------------DOUBLE WELL POTENTIAL-----------------------------
#
#                                                 m^2*w_b^4
# V(x,w_b,v0) =  - 0.5 *m*w_b^2*(x-delta)^2 +    ---------- x^4
#                                                  16V0
#
# ---------------------------------------------------------


class DoubleWell_with_explicit_bath_driver(Dummy_driver):
    r"""Adds to the double well potential an explicit harmonic bath. First dof correpond to the DW, the rest to the bath discretization

    !  V(q,x1,..,xn) = DW(q) +
    !      sum_2^(3*N) [ 0.5 m w_i^2(q - (c_i s(q))/(m w_i)^2)^2 ]
    !      s(q) = q *sd(q)
    !      sd(q) = [1+eps1*exp( q^2 / (2 deltaQ^2) ) ] + eps2 tanh(q/deltaQ)
    !      If eps1=eps2=0 then sd(q) =1 and s(q) = q --->Spatially independent bath

    DW+explicit_bath driver expects 9 arguments.
        Example: python driver.py -m DoubleWell_with_explicit_bath -o omega_b (cm^-1) V0 (cm^-1) mass delta(\AA) eta0  eps1 eps2  deltaQ omega_c(cm^-1)     \n
        python driver.py -m DoubleWell -o 500,2085,1837,0.00,1,0,0,1,500
    """

    def __init__(
        self,
        w_b=None,
        v0=None,
        m=None,
        delta=None,
        eta0=None,
        eps1=None,
        eps2=None,
        deltaQ=None,
        w_c=None,
        *args,
        **kwargs
    ):

        self.init = False

        try:
            w_b = w_b * invcm2au
            v0 = v0 * invcm2au
            self.m = m
            self.delta = delta * A2au

            self.bath_parameters = {}
            self.bath_parameters["m"] = self.m

            self.bath_parameters["eta0"] = eta0
            self.bath_parameters["eps1"] = eps1
            self.bath_parameters["eps2"] = eps2
            self.bath_parameters["deltaQ"] = deltaQ
            self.bath_parameters["w_c"] = w_c * invcm2au

        except:
            sys.exit(self.__doc__)

        self.A = -0.5 * self.m * (w_b) ** 2
        self.B = ((self.m**2) * (w_b) ** 4) / (16 * v0)
        super().__init__(*args, **kwargs)

    def __call__(self, cell, pos):
        """DoubleWell potential"""
        pot = 0
        q = pos.reshape(-1, 1)[0]
        x = pos.reshape(-1, 1)[1:]
        fq = np.zeros(q.shape)
        fx = np.zeros(x.shape)

        if not self.init:
            self.nbath = np.size(x)
            self.bath = Harmonic_Bath_explicit(self.nbath, **self.bath_parameters)

        # Harmonic bath
        pot, fq, fx = self.bath(q, x)

        # DW
        pot += self.A * (q - self.delta) ** 2 + self.B * (q**4)
        fq += -2.0 * self.A * (q - self.delta) - 4.0 * self.B * (q**3)

        force = np.concatenate((fq, fx.flatten())).reshape(pos.shape)
        vir = cell * 0.0  # makes a zero virial with same shape as cell
        extras = "empty"

        return pot, force, vir, extras

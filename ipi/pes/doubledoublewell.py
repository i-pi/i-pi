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

__DRIVER_NAME__ = "double_double_well"
__DRIVER_CLASS__ = "DDW_with_explicit_bath_driver"


invcm2au = units.unit_to_internal("frequency", "inversecm", 1.0)
A2au = units.unit_to_internal("length", "angstrom", 1.0)

# ------------DOUBLE WELL POTENTIAL-----------------------------
#
#                                                 m^2*w_b^4
# V(x,w_b,v0) =  - 0.5 *m*w_b^2*(x-delta)^2 +    ---------- x^4
#                                                  16V0
#
# ---------------------------------------------------------


class DDW_with_explicit_bath_driver(Dummy_driver):
    r"""Adds to a double-double well (DDW) potential coupled to two (explicit) harmonic baths.
           pos[0:2] = DDW
           pos[2:n//2+1] = bath1
           pos[n//2+1:] = bath2
           where n = 3*natoms = # dofs

    !  V(q1,q2,x1,..,xm,y1,...,ym) =
    !      DDW(q1,q2) +
    !      sum_i^(m) [ 0.5 m w_i^2(xi - (c_i s(q1))/(m w_i)^2)^2 ] +
    !      sum_i^(m) [ 0.5 m w_i^2(yi - (c_i s(q2))/(m w_i)^2)^2 ] +
    !
    !      with
    !      s(q) = q *sd(q)
    !      sd(q) = [1+eps1*exp( q^2 / (2 deltaQ^2) ) ] + eps2 tanh(q/deltaQ)
    !      If eps1=eps2=0 then sd(q) =1 and s(q) = q --->Spatially independent bath
    !
    !      and
    !
    !      DDW(q1,q2) = DW(q1) + DW(q2) + C(q1q2)^2

    DW+explicit_bath driver expects 11 arguments.\n
    Example:
        i-pi-py_driver.py -m DoubleWell_with_explicit_bath -o wb1 (cm^-1) V1 (cm^-1) wb2 (cm^-1) V2 (cm^-1) coupling(au) mass delta(\AA) eta0  eps1 eps2  deltaQ omega_c(cm^-1)     \n
        i-pi-py_driver.py -m DoubleWell -o 500,2085,500,2085,0.1,1837,0.00,1,0,0,1,500
    """

    def __init__(
        self,
        wb1=None,
        v1=None,
        wb2=None,
        v2=None,
        C=None,
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
            wb1 = wb1 * invcm2au
            v1 = v1 * invcm2au
            wb2 = wb2 * invcm2au
            v2 = v2 * invcm2au
            self.C = C
            self.m = m
            self.delta = delta * A2au

            self.bath_parameters = {}
            self.bath_parameters["m"] = self.m
            self.bath_parameters["eta0"] = eta0
            self.bath_parameters["eps1"] = eps1
            self.bath_parameters["eps2"] = eps2
            self.bath_parameters["deltaQ"] = deltaQ
            self.bath_parameters["w_c"] = w_c * invcm2au
            self.bath_parameters["delta"] = self.delta
        # def __init__(self, nbath, m, eta0, eps1, eps2, deltaQ, w_c, *args, **kwargs):

        except:
            print("Received arguments:")
            sys.exit(self.__doc__)

        self.A1 = -0.5 * self.m * (wb1) ** 2
        self.B1 = ((self.m**2) * (wb1) ** 4) / (16 * v1)

        self.A2 = -0.5 * self.m * (wb2) ** 2
        self.B2 = ((self.m**2) * (wb2) ** 4) / (16 * v2)
        super().__init__(*args, **kwargs)

    def __call__(self, cell, pos):
        """DoubleWell potential"""
        if not self.init:
            self.ndof = np.size(pos)
            assert self.ndof % 2 == 0, "Sorry we need an even number of ndof"
            self.nbath = (self.ndof - 2) // 2
            self.bath1 = Harmonic_Bath_explicit(self.nbath, self.bath_parameters)
            self.bath2 = Harmonic_Bath_explicit(self.nbath, self.bath_parameters)
            self.init = True

        pot = 0
        q1 = pos.reshape(-1, 1)[0:1]
        q2 = pos.reshape(-1, 1)[1:2]
        x = pos.reshape(-1, 1)[2 : self.ndof // 2 + 1]
        y = pos.reshape(-1, 1)[self.ndof // 2 + 1 :]
        assert self.ndof == q1.size + q2.size + x.size + y.size

        fx = np.zeros(x.shape)
        fy = np.zeros(x.shape)

        # Harmonic bath
        pot_x, fq1, fx = self.bath1(q1, x)
        pot_y, fq2, fy = self.bath2(q2, y)

        pot = pot_x + pot_y

        # DW
        pot += self.A1 * (q1 - self.delta) ** 2 + self.B1 * (q1**4)
        fq1 += -2.0 * self.A1 * (q1 - self.delta) - 4.0 * self.B1 * (q1**3)

        pot += self.A2 * (q2 - self.delta) ** 2 + self.B2 * (q2**4)
        fq2 += -2.0 * self.A2 * (q2 - self.delta) - 4.0 * self.B2 * (q2**3)

        # Coupling
        pot += self.C * q1 * q2
        fq1 += -self.C * q2
        fq2 += -self.C * q1

        force = np.concatenate((fq1, fq2, fx, fy), axis=None).reshape(pos.shape)
        vir = cell * 0.0  # makes a zero virial with same shape as cell
        extras = "empty"

        return pot, force, vir, extras

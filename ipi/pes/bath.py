""" Harmonic Bath """

import numpy as np

__DRIVER_NAME__ = None
__DRIVER_CLASS__ = "driver_tools"


class Harmonic_Bath_explicit(object):
    """Explicit description of an Harmonic bath"""

    def __init__(self, nbath, parameters):
        self.nbath = nbath
        self.m = parameters["m"]
        # self.delta = parameters["delta"]
        self.eta0 = parameters["eta0"]
        self.eps1 = parameters["eps1"]
        self.eps2 = parameters["eps2"]
        self.deltaQ = parameters["deltaQ"]
        self.w_c = parameters["w_c"]

        self.set_ci_and_wi()

    def set_ci_and_wi(self):
        """Computes the the parameters to represent the harmonic bath"""

        omega = np.zeros(self.nbath)
        self.omega2 = np.zeros(self.nbath)
        self.coef = np.zeros(self.nbath)

        if self.w_c > 0:
            for i in range(self.nbath):
                omega[i] = -self.w_c * np.log((i + 0.5) / self.nbath)
                self.omega2[i] = omega[i] ** 2
                self.coef[i] = (
                    omega[i]
                    * ((2 * self.eta0 * self.m * self.w_c) / (self.nbath * np.pi))
                    ** 0.5
                )
        self.init = True

    def S(self, q):
        """S coupling function"""
        return q * self.SD(q)

    def SD(self, q):
        """Auxiliary function to compute friction tensor"""
        dx = q / self.deltaQ
        SD = 1.0 + self.eps1 * np.exp(-0.5 * (dx**2)) + self.eps2 * np.tanh(dx)
        return SD

    def dSD_dq(self, q):
        """Auxiliary function to compute friction tensor"""
        dx = q / self.deltaQ
        dsddq1 = self.eps1 * np.exp(-0.5 * (dx**2)) * (-dx / self.deltaQ)
        dsddq2 = self.eps2 * (1 - np.tanh(dx) ** 2) / self.deltaQ
        dSD_dq = q * (dsddq1 + dsddq2) + self.SD(q)

        return dSD_dq

    def __call__(self, q, x):
        pot_bath = 0
        fq_bath = np.zeros(q.shape)
        fx = np.zeros(x.shape)
        assert x.size == self.nbath, "x.size {} , nbath {}".format(x.size, self.nbath)

        if self.w_c > 0:
            for i in range(self.nbath):
                aux = x[i] - (self.coef[i] * self.S(q) / (self.m * self.omega2[i]))
                pot_bath += 0.5 * self.m * self.omega2[i] * aux**2
                fq_bath += aux * self.coef[i] * self.dSD_dq(q)
                fx[i] -= (self.m * self.omega2[i] * aux).flatten()

        return pot_bath, fq_bath, fx

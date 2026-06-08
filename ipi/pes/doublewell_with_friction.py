"""Harmonic potential"""

try:
    from .doublewell import DoubleWell_driver
except:
    from doublewell import DoubleWell_driver

import numpy as np
from ipi.utils import units
import json
import sys

__DRIVER_NAME__ = "DW_friction"
__DRIVER_CLASS__ = "DoubleWell_with_friction_driver"


invcm2au = units.unit_to_internal("frequency", "inversecm", 1.0)
A2au = units.unit_to_internal("length", "angstrom", 1.0)

# ------------DOUBLE WELL POTENTIAL-----------------------------
#
#                                                 m^2*w_b^4
# V(x,w_b,v0) =  - 0.5 *m*w_b^2*(x-delta)^2 +    ---------- x^4
#                                                  16V0
#
# ---------------------------------------------------------


class DoubleWell_with_friction_driver(DoubleWell_driver):
    r"""Adds to the double well potential the calculation of the friction tensor.

    friction(q) = eta0 [\partial sd(q) \partial q ]^2
    with
    q = position, and
    sd(q) = [1+eps1 exp( (q-0)^2 / (2deltaQ^2) ) ] + eps2 tanh(q/deltaQ)

    DW+fric driver expects 8 arguments.
        Example: python driver.py -m DoubleWell_with_fric -o omega_b (cm^-1) V0 (cm^-1) mass eta0 eps1 eps2 delta(\AA) deltaQ
        python driver.py -m DoubleWell -o 500,2085,1837,0.00,1,0,0,1
    """

    def __init__(
        self,
        w_b=None,
        v0=None,
        m=None,
        eta0=None,
        eps1=None,
        eps2=None,
        delta=None,
        deltaQ=None,
        *args,
        **kwargs
    ):

        try:
            self.eta0 = float(eta0)
            self.root_eta0 = np.sqrt(self.eta0)
            self.eps1 = float(eps1)
            self.eps2 = float(eps2)
            self.deltaQ = float(deltaQ)

        except:
            sys.exit(self.__doc__)

        super().__init__(w_b=w_b, v0=v0, m=m, delta=delta, *args, **kwargs)

    def check_dimensions(self, pos):
        """Functions that checks dimensions of the received position"""
        assert pos.ndim == 2 and pos.shape[1] == 3, "We expect pos.shape (natoms,3), but we have {}".format(
            pos.shape
        )

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

    def get_diffusion_coefficient(self, pos):
        """Function that computes the array of diffusion coefficients."""
        self.check_dimensions(pos)
        natoms = pos.shape[0]
        diffusion_coefficient = np.zeros((natoms, 3 * natoms))
        for iatom, q in enumerate(pos[:, 0]):
            #diffusion_coefficient[iatom, 3 * iatom] = self.root_eta0 * self.dSD_dq(q)
            
            #If using A matrix style implementation, sigma should not contain friction strength
            diffusion_coefficient[iatom, 3 * iatom] =  self.dSD_dq(q)  
        return diffusion_coefficient

    def get_friction_coupling(self, pos):
        """Returns the separable coupling g(q)=q*SD(q) for each bath channel."""
        self.check_dimensions(pos)
        return np.asarray([q * self.SD(q) for q in pos[:, 0]], dtype=float)

    def get_scaled_friction_coupling(self, pos):
        """Returns sqrt(eta0) * g(q), consistent with diffusion_coefficient."""
        return self.root_eta0 * self.get_friction_coupling(pos)

    def get_friction_tensor(self, pos):
        """Function that computes spatially dependent friction tensor"""

        self.check_dimensions(pos)
        natoms = pos.shape[0]
        friction_tensor = np.zeros((3 * natoms, 3 * natoms))
        for iatom, q in enumerate(pos[:, 0]):
            friction_tensor[3 * iatom, 3 * iatom] = self.eta0 * self.dSD_dq(q) ** 2
        return friction_tensor

    def get_diffusion_and_friction(self, pos):
        """Function that computes the vector of diffusion coefficients
        and its outer product with itself, i.e., the static friction tensor.
        """
        diffusion_coefficient = self.get_diffusion_coefficient(pos)
        friction_tensor = diffusion_coefficient.T @ diffusion_coefficient

        return diffusion_coefficient, friction_tensor

    def compute_structure(self, cell, pos):
        """DoubleWell potential l"""

        pot, force, vir, extras = super(
            DoubleWell_with_friction_driver, self
        ).compute_structure(cell, pos)

        diffusion_coefficient, friction_tensor = self.get_diffusion_and_friction(pos)
        friction_coupling = self.get_friction_coupling(pos)
        scaled_friction_coupling = self.get_scaled_friction_coupling(pos)

        extras = json.dumps(
            {
                "friction": friction_tensor.tolist(),
                "diffusion_coefficient": diffusion_coefficient.tolist(),
                "friction_coupling": friction_coupling.tolist(),
                "scaled_friction_coupling": scaled_friction_coupling.tolist(),
            }
        )
        return pot, force, vir, extras

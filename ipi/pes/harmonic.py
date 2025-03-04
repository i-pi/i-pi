"""Harmonic potential"""

from .dummy import Dummy_driver
import numpy as np

__DRIVER_NAME__ = "harmonic"
__DRIVER_CLASS__ = "Harmonic_driver"


class Harmonic_driver(Dummy_driver):
    """
    Harmonic driver, generating either an isotropic or 3D confining
    potential for each atom.

    Command-line:
        i-pi-py_driver -m harmonic -u -o 1.3 [...]
        i-pi-py_driver -m harmonic -u -o 1.3,2.1,2.3 [...]

    Init parameters:
        :param k1: float, spring constant of the oscillator (in a.u.)
        :param k2: float, spring constant of the oscillator along y (in a.u.), optional
        :param k3: float, spring constant of the oscillator along z (in a.u.), optional
    """

    def __init__(self, k1, k2=None, k3=None, *args, **kwargs):

        if k2 == None or k3 == None:
            self.k = k1
            self.type = "isotropic"
        else:
            self.k = np.asarray([k1, k2, k3])
            self.type = "non-isotropic"

        super().__init__(*args, **kwargs)

    def __call__(self, cell, pos):
        """Silly harmonic potential"""
        if self.type == "isotropic":
            pot = 0.5 * self.k * (pos**2).sum()
            force = -self.k * pos
            vir = cell * 0.0  # makes a zero virial with same shape as cell
            extras = "nada"
        else:
            pot = 0
            pos3 = pos.reshape(-1, 3)
            force3 = np.zeros(pos.shape)
            for i in range(3):
                pot += 0.5 * self.k[i] * (pos3[:, i] ** 2).sum()
                force3[:, i] = (
                    -self.k[i] * pos3[:, i]
                )  # makes a zero force with same shape as pos
            vir = cell * 0.0  # makes a zero virial with same shape as cell
            extras = "nada"
            force = force3.reshape(pos.shape)
        return pot, force, vir, extras

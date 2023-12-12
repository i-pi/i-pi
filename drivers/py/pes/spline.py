""" spline potential """

import numpy as np
import sys
from .dummy import Dummy_driver
from scipy import interpolate
import json

"""Spline driver. This is not a serious interpolation, use it if you know what you are doing. """
factor_coord = 5
mass = 1836
friction = False
SI = True
fric_value = 0.165

__DRIVER_NAME__ = "spline"
__DRIVER_CLASS__ = "Spline_driver"


class Spline_driver(Dummy_driver):
    def __init__(self, args=None, verbose=None):
        self.error_msg = """\nspline driver requires specification of filename that contains 5 columns (pos, f1,f2,f3,e) to perform 3x1D spline.\nExample: python driver.py -m spline -u -o <filename>\n"""
        super(Spline_driver, self).__init__(args)
        self.get_spline()
        if friction and not SI:
            self.get_spline_fric()
        self.k = 1836 * (3800.0 / 219323.0) ** 2

    def check_arguments(self):
        """Function that checks the arguments required to run the driver"""

        try:
            data = np.loadtxt(self.args[0]).reshape(-1, 5)
        except ValueError:
            sys.exit(self.error_msg)
        self.data = data

    def get_spline(self):
        """Function that creates the 1D spline and its derivative"""

        self.spline_f = []
        for i in range(3):
            self.spline_f.append(
                interpolate.interp1d(self.data[:, 0], self.data[:, i + 1], kind="cubic")
            )

        self.spline_e = interpolate.interp1d(
            factor_coord * self.data[:, 0], self.data[:, 4], kind="cubic"
        )

    def get_spline_fric(self):
        """Function that creates the 1D spline for the friction"""
        data = np.loadtxt(self.args[1]).reshape(-1, 10)
        self.spline_fric = []
        for i in range(9):
            self.spline_fric.append(
                interpolate.interp1d(
                    factor_coord * data[:, 0], mass * data[:, i + 1], kind="cubic"
                )
            )

    def get_energy(self, pos):
        x = self.full2oneD(pos)

        pot = self.spline_e(x)
        pot += 0.5 * self.k * (pos[0, 1] ** 2 + pos[0, 2] ** 2)
        return pot

    def get_forces(self, pos):
        x = self.full2oneD(pos)
        force = np.zeros(pos.shape)
        d = 0.001
        force[0, 0] = -(self.spline_e(x + d) - self.spline_e(x - d)) / (2 * d)
        # force[0, 0] = self.spline_f[0](x) #+ self.spline_f[1](x) + self.spline_f[2](x)
        force[0, 1] = -self.k * pos[0, 1]
        force[0, 2] = -self.k * pos[0, 2]

        return force

    def get_friction(self, pos):
        x = self.full2oneD(pos)

        friction_tensor = np.zeros(9)

        if SI:
            friction_tensor = np.eye(3) * fric_value
        else:
            for i in range(9):
                friction_tensor[i] = self.spline_fric[i](x)
            friction_tensor = friction_tensor.reshape((3, 3))

        w = np.linalg.eigvals(friction_tensor)
        assert (w >= 0).all()
        assert (friction_tensor - friction_tensor.T == np.zeros((3, 3))).all()
        return friction_tensor

    def full2oneD(self, pos):
        """Function that gets the 1D coordinates from the pos vector"""
        return factor_coord * pos[0, 0]

    def check_dimensions(self, pos):
        """Functions that checks dimensions of the received position"""
        assert pos.shape == (1, 3)

    def __call__(self, cell, pos):
        """Evaluate energy, forces and friction"""
        self.check_dimensions(pos)
        vir = cell * 0.0  # makes a zero virial with same shape as cell

        pot = self.get_energy(pos)
        force = self.get_forces(pos)

        if friction:
            friction_tensor = self.get_friction(pos)
            extras = json.dumps({"friction": friction_tensor.tolist()})
        else:
            extras = ""

        return pot, force, vir, extras

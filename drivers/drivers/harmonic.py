""" Harmonic potential """

import sys
from .dummy import Dummy_driver


class Harmonic_driver(Dummy_driver):
    def __init__(self, args=None):

        self.error_msg = """\nHarmonic driver requires specification of force constant.\nExample: python driver.py -m harmonic -u -o 1.3\n"""
        super(Harmonic_driver, self).__init__(args)

    def check_arguments(self):
        """ Function that checks the arguments required to run the driver """
        try:
            k = list(map(float, self.args.split()))
        except ValueError:
            sys.exit(self.error_msg)

        if len(k) == 1:
            self.k = k[0]
        else:
            sys.exit(self.error_msg)

    def __call__(self, cell, pos):
        """ Silly harmonic potential"""
        pot = 0.5 * self.k * (pos ** 2).sum()
        force = -self.k * pos  # makes a zero force with same shape as pos
        vir = cell * 0.0  # makes a zero virial with same shape as cell
        extras = "nada"
        return pot, force, vir, extras

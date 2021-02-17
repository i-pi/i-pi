""" Harmonic potential """

import sys
import numpy as np
from .dummy import Dummy_driver

try:
    from rascal.models.IP_ipi_interface import IPICalculator as IPICalc
except:
    IPICalc = None


class Rascal_driver(Dummy_driver):
    def __init__(self, args=None):

        self.error_msg = """Rascal driver requires specification of a .json model file fitted with librascal, 
                            and a template file that describes the chemical makeup of the structure. 
                            Example: python driver.py -m rascal -u -o model.json,template.xyz"""

        super().__init__(args)

        if IPICalc is None:
            raise ImportError("Couldn't load librascal bindings")

    def check_arguments(self):
        """ Function that checks the arguments required to run the driver """
        try:
            arglist = self.args.split(",")
        except ValueError:
            sys.exit(self.error_msg)

        print("ARGLIST ", arglist)
        if len(arglist) == 2:
            self.model = arglist[0]
            self.template = arglist[1]
        else:
            raise
            sys.exit(self.error_msg)

        self.rascal_calc = IPICalc(self.model, self.template)

    def __call__(self, cell, pos):
        """ Silly harmonic potential, with unit frequency in a.u."""

        cell_tuple = np.zeros((2, 3, 3), float)
        cell_tuple[0] = cell
        cell_tuple[1] = cell
        pot, force, vir = self.rascal_calc.calculate(pos, cell_tuple)
        extras = "nada"
        return pot, force, vir, extras

import sys

try:
    from .dummy import Dummy_driver
except:
    from dummy import Dummy_driver

import numpy as np
from ipi.utils import units

# Unit conversions
A2au    = units.unit_to_internal("length","angstrom", 1.0)
ev2au   = units.unit_to_internal("energy","electronvolt", 1.0)


__DRIVER_NAME__  = "MorseHarmonic"
__DRIVER_CLASS__ = "MorseHarmonic_driver"

class MorseHarmonic_driver(Dummy_driver):
    """
    Morse + harmonic driver accepts 0 or 4 arguments.
    Example: i-pi-py_driver -m MorseHarmonic \
             -o De(eV) a(1/Angstrom) z0(Angstrom) k(a.u)
    """
    def __init__(self, De=None, a=None, z0=None, k=None, *args, **kwargs):
        # Default parameters if none provided (from JPCC 2024, typical H-metal adsorption)
        if De is None or a is None or m is None or z0 is None or k is None:
            # Set default values
            
            De     = 0.00735         # well depth in Hartree 
            a      = 1.61            # range parameter in a.u.^-1
            z0     = 2.1             # equilibrium position in a.u. (~1.1 Å)
            k      = 0.343            # harmonic force constant

            # Print default Morse parameters
            print("Using default Morse parameters for H adsorption:")
            print(f"  De (depth) = {De:.6f} Ha ({De/ev2au:.6f} eV)")
            print(f"  a (range)  = {a:.6f} a.u.^-1")
            print(f"  z0 (eq.)   = {z0:.6f} a.u. ({z0/A2au:.6f} Å)")
        else:
            try:
                # Convert user inputs: De from eV, a from 1/Å, z0 from Å
                De      = De  * ev2au
                a       = a   / A2au
                z0      = z0 * A2au
            except:
                sys.exit(self.__doc__)

        # Harmonic spring constant (to match reference 3D harmonic in x & y)
        

        # Store Morse parameters
        self.De = De
        self.a  = a
        self.z0 = z0
        self.k  = k
        super().__init__(*args, **kwargs)

    def __call__(self, cell, pos):
        """Compute total potential and forces: Morse in z, harmonic in x & y"""
        pos3   = pos.reshape(-1, 3)
        force3 = np.zeros_like(pos3)
        pot    = 0.0

        # Morse potential along z (idx 2)
        z        = pos3[:, 2]
        exp_term = np.exp(-self.a * (z - self.z0))
        U_morse  = self.De * (1.0 - exp_term) ** 2
        pot     += U_morse.sum()
        # Force in z: -dU/dz
        force3[:, 2] = -2.0 * self.De * self.a * exp_term * (1.0 - exp_term)

        # Independent harmonic potentials in x (idx 0) & y (idx 1)
        for i in (0, 1):
            coord      = pos3[:, i]
            pot       += 0.5 * self.k * (coord ** 2).sum()
            force3[:, i] = -self.k * coord

        # Zero virial and dummy extras
        vir    = cell * 0.0
        extras = "empty"

        # Reshape forces back to original shape
        return pot, force3.reshape(pos.shape), vir, extras


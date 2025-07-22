try:
    from .dummy import Dummy_driver
except:
    from dummy import Dummy_driver

import numpy as np
from ipi.utils import units

# Unit conversions
A2au = units.unit_to_internal("length", "angstrom", 1.0)
ev2au = units.unit_to_internal("energy", "electronvolt", 1.0)


__DRIVER_NAME__ = "MorseHarmonic"
__DRIVER_CLASS__ = "MorseHarmonic_driver"


class MorseHarmonic_driver(Dummy_driver):
    """
    Quartic expansion of Morse potential around equilibrium + harmonic driver.

    Accepts 0 or 4 arguments.

    Example usage:
        i-pi-py_driver -m MorseHarmonic -o a=a_value,k=k_value,z0=z0_value,De=De_value

    Parameter descriptions (with expected units):

    - De (Ha)         : Well depth in Hartree. Default: 0.00735 Ha
    - a  (a.u.^-1)    : Range parameter in inverse Bohr. Default: 1.61 a.u.^-1
    - z0 (a.u.)       : Equilibrium position in Bohr (~1.1 Å). Default: 2.1 a.u.
    - k  (a.u.)       : Harmonic force constant. Default: 0.343 a.u.
    """

    def __init__(self, De=None, a=None, z0=None, k=None, *args, **kwargs):
        # Default parameters if none provided (from JPCC 2024, typical H-metal adsorption)
        if De is None or a is None or z0 is None or k is None:
            # Set default values

            De = 0.00735  # well depth in Hartree
            a = 1.61  # range parameter in a.u.^-1
            z0 = 2.1  # equilibrium position in a.u. (~1.1 Å)
            k = 0.343  # harmonic force constant

            # Print default Morse parameters
            print("Using default Morse parameters for H adsorption:")
            print(f"  De (depth) = {De:.6f} Ha ({De/ev2au:.6f} eV)")
            print(f"  a (range)  = {a:.6f} a.u.^-1")
            print(f"  z0 (eq.)   = {z0:.6f} a.u. ({z0/A2au:.6f} Å)")
        else:
            try:
                # Convert user inputs: De from eV, a from 1/Å, z0 from Å, k from eV/Å^2
                De = De * ev2au
                a = a / A2au
                z0 = z0 * A2au
                k = k * ev2au / A2au**2
            except:
                raise ValueError(
                    "Error: Invalid input parameters for MorseHarmonic_driver."
                )

        # Store Morse parameters
        self.De = De
        self.a = a
        self.z0 = z0
        self.k = k
        super().__init__(*args, **kwargs)

    def potential(self, pos: np.ndarray):
        pot = np.zeros(pos.shape[:-1])
        pos3 = pos.reshape(-1, 3)
        pot3 = np.reshape(pot, -1)
        xy = pos3[:, :2]
        z = (pos3[:, 2] - self.z0) * self.a
        pot3[:] = self.De * (z**2 - z**3 + (7.0 / 12.0) * z**4)
        pot3 += self.k / 2 * np.sum(xy**2, axis=-1)
        return pot

    def force(self, pos: np.ndarray):
        force = np.zeros_like(pos)
        pos3 = pos.reshape(-1, 3)
        force3 = np.reshape(force, pos3.shape)
        xy = pos3[:, :2]
        z = (pos3[:, 2] - self.z0) * self.a
        force3[:, :2] = -self.k * xy
        force3[:, 2] = -self.De * self.a * (2 * z - 3 * z**2 + (7.0 / 3.0) * z**3)
        return force

    def both(self, pos: np.ndarray):
        force = np.zeros_like(pos)
        pot = np.zeros(pos.shape[:-1])
        pos3 = pos.reshape(-1, 3)
        pot3 = np.reshape(pot, -1)
        force3 = np.reshape(force, pos3.shape)
        xy = pos3[:, :2]
        z = (pos3[:, 2] - self.z0) * self.a
        pot3[:] = self.De * (z**2 - z**3 + (7.0 / 12.0) * z**4)
        pot3 += self.k / 2 * np.sum(xy**2, axis=-1)
        force3[:, :2] = -self.k * xy
        force3[:, 2] = -self.De * self.a * (2 * z - 3 * z**2 + (7.0 / 3.0) * z**3)
        return pot, force

    def __call__(self, cell: np.ndarray, pos: np.ndarray):
        """Compute total potential and forces: Morse in z, harmonic in x & y"""
        pot, force = self.both(pos)
        # Zero virial and dummy extras
        vir = cell * 0.0
        extras = "empty"
        # Reshape forces back to original shape
        return np.sum(pot), force, vir, extras

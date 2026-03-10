"""Creates objects that deal with classical simulations.

Generates an atoms class either from a set of positions and momenta.
This class is only used if no beads tag is present in the xml file.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np

from ipi.engine.atoms import *
from ipi.utils.inputvalue import *
from ipi.utils.depend import *


__all__ = ["InputAtoms"]


class InputAtoms(Input):
    """Atoms input class.

    Handles generating the appropriate atoms class from the xml input file,
    and generating the xml checkpoint tags and data from an instance of the
    object.

    Attributes:
       natoms: An optional integer giving the number of atoms. Defaults to 0.
       q: An optional array giving the atom positions. Defaults to an empty
          array with no elements.
       p: An optional array giving the atom momenta. Defaults to an empty
          array with no elements.
       m: An optional array giving the atom masses. Defaults to an empty
          array with no elements.
       names: An optional array giving the atom names. Defaults to an empty
          array with no elements
    """

    fields = {
        "natoms": (
            InputValue,
            {"dtype": int, "default": 0, "help": "The number of atoms."},
        ),
        "q": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "The positions of the atoms, in the format [x1, y1, z1, x2, ... ].",
                "dimension": "length",
            },
        ),
        "p": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "The momenta of the atoms, in the format [px1, py1, pz1, px2, ... ].",
                "dimension": "momentum",
            },
        ),
        "m": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "The masses of the atoms, in the format [m1, m2, ... ].",
                "dimension": "mass",
            },
        ),
        "names": (
            InputArray,
            {
                "dtype": str,
                "default": input_default(
                    factory=np.zeros, args=(0,), kwargs={"dtype": np.dtype("|U6")}
                ),
                "help": "The names of the atoms, in the format [name1, name2, ... ].",
            },
        ),
    }

    default_help = "Deals with a single replica of the system or classical simulations."
    default_label = "ATOMS"

    def store(self, atoms):
        """Takes an Atoms instance and stores a minimal representation of it.

        Args:
           atoms: An Atoms object from which to initialise from.
           filename: An optional string giving a filename to take the atom
              positions from. Defaults to ''.
        """

        super(InputAtoms, self).store()
        self.natoms.store(atoms.natoms)
        self.q.store(dstrip(atoms.q))
        self.p.store(dstrip(atoms.p))
        self.m.store(dstrip(atoms.m))
        self.names.store(dstrip(atoms.names))

    def fetch(self):
        """Creates an atoms object.

        Returns:
           An atoms object of the appropriate type and with the appropriate
           properties given the attributes of the InputAtoms object.
        """

        super(InputAtoms, self).fetch()
        atoms = Atoms(self.natoms.fetch())
        atoms.q = self.q.fetch()
        atoms.p = self.p.fetch()
        atoms.m = self.m.fetch()
        atoms.names = self.names.fetch()
        return atoms

    def write(self, name="", indent=""):
        """Overloads Input write() function so that nothing is written if
        no atoms are present. This occurs if the beads object has been specified,
        so that the classical atoms object is not initialized.

        Returns:
           A string giving the appropriate xml tags for the checkpoint file.
        """

        if self.natoms.fetch() > 0:
            return super(InputAtoms, self).write(name=name, indent=indent)
        else:
            return ""

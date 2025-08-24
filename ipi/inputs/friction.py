import numpy as np

from ipi.utils.inputvalue import (
    Input,
    InputValue,
    InputArray,
    input_default,
)

class InputFriction(Input):
    attribs = {
        "spectral_density": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.ones, args=(0,)),
                "help": "A two column data is expected. First column: w (cm^-1) frequency. Second column: J(w) spectral density. See Eq. 6 in Phys. Rev. Lett. 134,226201(2025).",
            }
        ),
        "frequency": (
            InputValue,
            {
                "dtype": float,
                "default": 0.0,
                "help": "Energy at which the friction tensor is evaluated in the client code",
                "dimension": "energy",
            },
        ),
    }

    default_help = "Simulates the elctronic friction"
    default_label = "FRICTION"

    def store(self,friction):
        """Takes a friction instance and store a minimal representation of it.
        
        Args: 
            friction: A friction object.
        """

        super(InputFriction, self).store(friction)
        self.spectral_density.store(friction.spectral_density)
        self.frequency.store(friction.frequency)
        
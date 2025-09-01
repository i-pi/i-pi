import numpy as np
from ipi.engine.friction import Friction
from ipi.utils.inputvalue import (
    Input,
    InputValue,
    InputArray,
    input_default,
)


class InputFriction(Input):
    attribs = {}

    fields = {
        "spectral_density": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.ones, args=(0,)),
                "help": "A two column data is expected. First column: w (cm^-1) frequency. Second column: J(w) spectral density. See Eq. 6 in Phys. Rev. Lett. 134,226201(2025).",
            },
        ),
    }

    default_help = "Simulates the elctronic friction"
    default_label = "FRICTION"

    def store(self, friction: Friction) -> None:
        """Takes a friction instance and store a minimal representation of it.

        Args:
            friction: A friction object.
        """

        super(InputFriction, self).store(friction)
        self.spectral_density.store(friction.spectral_density)
    
    def fetch(self) -> Friction:
        return Friction(
            spectral_density=self.spectral_density.fetch(),
        )
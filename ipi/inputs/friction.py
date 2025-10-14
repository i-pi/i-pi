import numpy as np
from ipi.engine.friction import Friction, FrictionProtocol
from ipi.utils.inputvalue import (
    Input,
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
        "alpha": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.ones, args=(0,)),
                "help": "A two column data is expected. First column: normal mode frequency. Second column: alpha. See Eq. 8b in Phys. Rev. Lett. 134,226201(2025).",
            },
        ),
    }

    default_help = "Simulates the electronic friction"
    default_label = "FRICTION"

    def store(self, friction: FrictionProtocol) -> None:
        """Takes a friction instance and store a minimal representation of it.

        Args:
            friction: A friction object.
        """

        super(InputFriction, self).store(friction)
        if isinstance(friction, Friction):
            self.spectral_density.store(friction.spectral_density)
            self.alpha.store(friction.alpha)

    def fetch(self) -> FrictionProtocol:
        return Friction(
            spectral_density=self.spectral_density.fetch(),
            alpha=self.alpha.fetch(),
        )

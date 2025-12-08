import numpy as np
from ipi.engine.friction import Friction
from ipi.utils.inputvalue import (
    Input,
    InputArray,
    input_default,
    InputValue,
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
        "alpha_input": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.ones, args=(0,)),
                "help": "A two column data is expected. First column: normal mode frequency. Second column: alpha. See Eq. 8b in Phys. Rev. Lett. 134,226201(2025).",
            },
        ),
        "position_dependent": (
            InputValue,
            {
                "dtype": bool,
                "default": False,
                "help": "If True, position dependent friction is used. Not implemented yet.",
            },
        ),
        "non_markovian": (
            InputValue,
            {
                "dtype": bool,
                "default": False,
                "help": "If True, non-markovian friction is used. Not implemented yet.",
            },
        ),
    }

    default_help = "Simulates the electronic friction"
    default_label = "FRICTION"

    def store(self, friction: Friction) -> None:
        """Takes a friction instance and store a minimal representation of it.

        Args:
            friction: A friction object.
        """

        super(InputFriction, self).store(friction)
        if isinstance(friction, Friction):
            self.spectral_density.store(friction.spectral_density)
            self.alpha_input.store(friction.alpha_input)
            self.position_dependent.store(friction.position_dependent)
            self.non_markovian.store(friction.non_markovian)

    def fetch(self) -> Friction:
        return Friction(
            spectral_density=self.spectral_density.fetch(),
            alpha_input=self.alpha_input.fetch(),
            position_dependent=self.position_dependent.fetch(),
            non_markovian=self.non_markovian.fetch(),
        )

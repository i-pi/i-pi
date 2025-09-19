import numpy as np
from ipi.engine.friction import Friction, FrictionProtocol
from ipi.utils.inputvalue import (
    Input,
    InputArray,
    InputValue,
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
        "analytical_equation": (
            InputValue,
            {
                "dtype": str,
                "default": "",
                "help": "",
            }
        ),
        "omega_cutoff": (
            InputValue,
            {
                "dtype": float,
                "default": 0.0,
                "help": "",
            }
        ),
        "eta":(
            InputValue,
            {
                "dtype": float,
                "default": 0.0,
                "help":"",
            }
        )
    }

    default_help = "Simulates the elctronic friction"
    default_label = "FRICTION"

    def store(self, friction: FrictionProtocol) -> None:
        """Takes a friction instance and store a minimal representation of it.

        Args:
            friction: A friction object.
        """

        super(InputFriction, self).store(friction)
        if isinstance(friction, Friction):
            self.spectral_density.store(friction.spectral_density)
            self.analytical_equation.store("")
        if hasattr(friction, "omega_cutoff"):
            self.omega_cutoff.store(friction.omega_cutoff)

    def fetch(self) -> FrictionProtocol:
        analytical_equation = self.analytical_equation.fetch()
        if analytical_equation:
            from ipi.utils.frictiontools import Fricition_eq133, Fricition_eq134
            if analytical_equation == "1.33":
                return Fricition_eq133(self.omega_cutoff.fetch(), self.eta.fetch())
            if analytical_equation == "1.34":
                return Fricition_eq134(self.omega_cutoff.fetch(), self.eta.fetch())
            raise ValueError(f"Unknown analytical equation for friction {analytical_equation}")
        return Friction(
            spectral_density=self.spectral_density.fetch(),
        )

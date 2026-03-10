import numpy as np

from ipi.engine.friction import Friction
from ipi.utils.inputvalue import Input, InputArray, InputValue, input_default


class InputFriction(Input):
    attribs = {}

    fields = {
        "variable_friction": (
            InputValue,
            {
                "dtype": bool,
                "default": True,
                "help": "If true, read position-dependent sigma from force extras; if false, use sigma_static.",
            },
        ),
        "bath_mode": (
            InputValue,
            {"dtype": str, "default": "non-markovian", "help": "..."},
        ),
        "debug_mf_mode": (
            InputValue,
            {"dtype": str, "default": "none", "help": "..."},
        ),
        "Lambda": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "...",
            },
        ),
        "debug_alpha_input": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "...",
            },
        ),
        "sigma_static": (
            InputValue,
            {
                "dtype": float,
                "default": 1.0,
                "help": "Constant linear coupling amplitude used when variable_friction=False.",
            },
        ),
        "sigma_key": (
            InputValue,
            {
                "dtype": str,
                "default": "sigma",
                "help": "Force-extras key for variable friction payload. Expected shape is (nbeads, nbath, 3*natoms) or (nbeads, nbath, 3*len(friction_atoms)).",
            },
        ),

        # 0-based atom indices to include in friction; empty means all atoms.
        "friction_atoms": (
            InputArray,
            {
                "dtype": int,
                "default": input_default(factory=lambda n: np.zeros(n, dtype=int), args=(0,)),
                "help": "0-based atom indices included in friction. Empty means all atoms. If set, driver may return sigma only for this subset (3*len(friction_atoms) DOFs), and the rest are treated as zero-friction.",
            },
        ),
    }

    default_help = "Friction operator configuration (MF + markovian/non-markovian bath). For variable friction, sigma must be provided in force extras."
    default_label = "FRICTION"

    def store(self, friction: Friction) -> None:
        super(InputFriction, self).store(friction)

        if not isinstance(friction, Friction):
            return

        self.variable_friction.store(friction.variable_friction)
        self.bath_mode.store(friction.bath_mode)
        self.debug_mf_mode.store(friction.debug_mf_mode)

        self.Lambda.store(friction.Lambda)
        self.debug_alpha_input.store(friction.debug_alpha_input)
        self.sigma_static.store(friction.sigma_static)

        self.sigma_key.store(friction.sigma_key)

        fa = getattr(friction, "friction_atoms", None)
        if fa is None:
            fa = np.zeros(0, dtype=int)
        self.friction_atoms.store(np.array(fa, dtype=int).flatten())

    def fetch(self) -> Friction:
        return Friction(
            variable_friction=self.variable_friction.fetch(),
            bath_mode=self.bath_mode.fetch(),
            debug_mf_mode=self.debug_mf_mode.fetch(),
            
            Lambda=self.Lambda.fetch(),
            debug_alpha_input=self.debug_alpha_input.fetch(),
            sigma_static=self.sigma_static.fetch(),

            sigma_key=self.sigma_key.fetch(),

            friction_atoms=self.friction_atoms.fetch(),
        )

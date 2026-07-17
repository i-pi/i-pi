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
            {
                "dtype": str,
                "default": "on",
                "help": "Caldeira-Leggett mean-field/counterterm toggle. Use 'on' to apply the MF force and 'off' to disable it. Legacy 'none' is treated as 'off'.",
            },
        ),
        "Lambda": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "...",
            },
        ),
        "Ap": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=((0, 0),)),
                "help": "Momentum-plus-auxiliary drift matrix for non-Markovian GLE, Ap = [[0, theta^T], [-theta, A]].",
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
                "help": "Force-extras key for variable friction payload. Expected shape is (nbeads, nbath, 3*natoms).",
            },
        ),
        "coupling_key": (
            InputValue,
            {
                "dtype": str,
                "default": "friction_coupling",
                "help": "Force-extras key for bead-resolved coupling values F(q), expected shape (nbeads, nbath).",
            },
        ),
        "coupling_mode": (
            InputValue,
            {
                "dtype": str,
                "default": "driver",
                "help": "How to obtain variable-friction coupling values. 'driver' reads coupling_key from extras; 'centroid_endpoint_trapezoid' computes bead-centroid endpoint coupling from bead Sigma and centroid Sigma.",
            },
        ),
        "centroid_sigma_key": (
            InputValue,
            {
                "dtype": str,
                "default": "centroid_sigma",
                "help": "Force-extras key for the centroid Sigma payload used by coupling_mode='centroid_endpoint_trapezoid'.",
            },
        ),
        "coupling_friction_atom": (
            InputValue,
            {
                "dtype": int,
                "default": -1,
                "help": "Optional 0-based atom index used for centroid-relative coupling displacements. If negative, infer from sigma_meta.friction_atoms when exactly one atom is present.",
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
        self.Ap.store(friction.Ap)
        self.debug_alpha_input.store(friction.debug_alpha_input)
        self.sigma_static.store(friction.sigma_static)

        self.sigma_key.store(friction.sigma_key)
        self.coupling_key.store(friction.coupling_key)
        self.coupling_mode.store(friction.coupling_mode)
        self.centroid_sigma_key.store(friction.centroid_sigma_key)
        self.coupling_friction_atom.store(friction.coupling_friction_atom)

    def fetch(self) -> Friction:
        return Friction(
            variable_friction=self.variable_friction.fetch(),
            bath_mode=self.bath_mode.fetch(),
            debug_mf_mode=self.debug_mf_mode.fetch(),
            
            Lambda=self.Lambda.fetch(),
            Ap=self.Ap.fetch(),
            debug_alpha_input=self.debug_alpha_input.fetch(),
            sigma_static=self.sigma_static.fetch(),

            sigma_key=self.sigma_key.fetch(),
            coupling_key=self.coupling_key.fetch(),
            coupling_mode=self.coupling_mode.fetch(),
            centroid_sigma_key=self.centroid_sigma_key.fetch(),
            coupling_friction_atom=self.coupling_friction_atom.fetch(),
        )

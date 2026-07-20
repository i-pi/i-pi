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
                "help": "If true, read a position-dependent canonical coupling Jacobian and Gamma from force extras; if false, use sigma_static.",
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
        "coupling_jacobian_key": (
            InputValue,
            {
                "dtype": str,
                "default": "friction_coupling_jacobian",
                "help": "Force-extras key for dF_channel/dQ with shape (nbeads, nchannels, nactive_dof).",
            },
        ),
        "gamma_key": (
            InputValue,
            {
                "dtype": str,
                "default": "friction_gamma",
                "help": "Force-extras key for canonical Gamma with shape (nbeads, nactive_dof, nactive_dof).",
            },
        ),
        "friction_meta_key": (
            InputValue,
            {
                "dtype": str,
                "default": "friction_meta",
                "help": "Force-extras key for canonical schema, active DOFs, channel labels, and units.",
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
                "help": "How to obtain variable-friction coupling values. 'driver' reads coupling_key; 'centroid_endpoint_trapezoid' uses canonical bead and centroid coupling Jacobians.",
            },
        ),
        "centroid_coupling_jacobian_key": (
            InputValue,
            {
                "dtype": str,
                "default": "centroid_friction_coupling_jacobian",
                "help": "Force-extras key for the centroid coupling Jacobian used by endpoint trapezoid coupling.",
            },
        ),
        "centroid_friction_meta_key": (
            InputValue,
            {
                "dtype": str,
                "default": "centroid_friction_meta",
                "help": "Force-extras key for centroid canonical friction metadata.",
            },
        ),
        "coupling_friction_atom": (
            InputValue,
            {
                "dtype": int,
                "default": -1,
                "help": "Optional 0-based atom index used for centroid-relative coupling displacements. If negative, infer from friction_meta.active_atoms when exactly one atom is present.",
            },
        ),
    }

    default_help = "Friction operator configuration (MF + Markovian/non-Markovian bath) using a canonical driver coupling Jacobian and Gamma."
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

        self.coupling_jacobian_key.store(friction.coupling_jacobian_key)
        self.gamma_key.store(friction.gamma_key)
        self.friction_meta_key.store(friction.friction_meta_key)
        self.coupling_key.store(friction.coupling_key)
        self.coupling_mode.store(friction.coupling_mode)
        self.centroid_coupling_jacobian_key.store(
            friction.centroid_coupling_jacobian_key
        )
        self.centroid_friction_meta_key.store(friction.centroid_friction_meta_key)
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

            coupling_jacobian_key=self.coupling_jacobian_key.fetch(),
            gamma_key=self.gamma_key.fetch(),
            friction_meta_key=self.friction_meta_key.fetch(),
            coupling_key=self.coupling_key.fetch(),
            coupling_mode=self.coupling_mode.fetch(),
            centroid_coupling_jacobian_key=self.centroid_coupling_jacobian_key.fetch(),
            centroid_friction_meta_key=self.centroid_friction_meta_key.fetch(),
            coupling_friction_atom=self.coupling_friction_atom.fetch(),
        )

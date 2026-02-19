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
        "mf_mode": (
            InputValue,
            {"dtype": str, "default": "reconstruct", "help": "..."},
        ),
        "Lambda": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "...",
            },
        ),
        "alpha_input": (
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
        "ou_nterms": (InputValue, {"dtype": int, "default": 4, "help": "..."}),
        "ou_tmax": (InputValue, {"dtype": float, "default": 200.0, "help": "..."}),
        "ou_nt": (InputValue, {"dtype": int, "default": 2000, "help": "..."}),
        "ou_print": (InputValue, {"dtype": bool, "default": True, "help": "..."}),
        "ou_propagator": (InputValue, {"dtype": str, "default": "exact", "help": "..."}),

        "sigma_key": (
            InputValue,
            {
                "dtype": str,
                "default": "sigma",
                "help": "Force-extras key for variable friction payload. Expected shape is (nbeads, nbath, 3*natoms) or (nbeads, nbath, 3*len(friction_atoms)).",
            },
        ),
        "friction_key": (InputValue, {"dtype": str, "default": "friction", "help": "..."}),


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
        self.mf_mode.store(friction.mf_mode)

        self.Lambda.store(friction.Lambda)
        self.alpha_input.store(friction.alpha_input)
        self.sigma_static.store(friction.sigma_static)

        self.ou_nterms.store(friction.ou_nterms)
        self.ou_tmax.store(friction.ou_tmax)
        self.ou_nt.store(friction.ou_nt)
        self.ou_print.store(friction.ou_print)
        self.ou_propagator.store(friction.ou_propagator)

        self.sigma_key.store(friction.sigma_key)
        self.friction_key.store(friction.friction_key)


        fa = getattr(friction, "friction_atoms", None)
        if fa is None:
            fa = np.zeros(0, dtype=int)
        self.friction_atoms.store(np.array(fa, dtype=int).flatten())

    def fetch(self) -> Friction:
        return Friction(
            variable_friction=self.variable_friction.fetch(),
            bath_mode=self.bath_mode.fetch(),
            mf_mode=self.mf_mode.fetch(),
            
            Lambda=self.Lambda.fetch(),
            alpha_input=self.alpha_input.fetch(),
            sigma_static=self.sigma_static.fetch(),

            ou_nterms=self.ou_nterms.fetch(),
            ou_tmax=self.ou_tmax.fetch(),
            ou_nt=self.ou_nt.fetch(),
            ou_print=self.ou_print.fetch(),
            ou_propagator=self.ou_propagator.fetch(),

            sigma_key=self.sigma_key.fetch(),
            friction_key=self.friction_key.fetch(),


            friction_atoms=self.friction_atoms.fetch(),
        )

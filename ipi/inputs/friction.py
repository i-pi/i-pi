import numpy as np

from ipi.engine.friction import Friction
from ipi.utils.inputvalue import Input, InputArray, InputValue, input_default


class InputFriction(Input):
    attribs = {}

    fields = {
        "use_linear_coupling": (
            InputValue,
            {"dtype": bool, "default": False, "help": "..."},
        ),
        "bath_mode": (
            InputValue,
            {"dtype": str, "default": "non-markovian", "help": "..."},
        ),
        "mf_mode": (
            InputValue,
            {"dtype": str, "default": "reconstruct", "help": "..."},
        ),
        "spectral_density": (
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
        "friction_static": (
            InputValue,
            {"dtype": float, "default": 1.0, "help": "..."},
        ),
        "ou_fit_kind": (InputValue, {"dtype": str, "default": "damped_cosine", "help": "..."}),
        "ou_nterms": (InputValue, {"dtype": int, "default": 4, "help": "..."}),
        "ou_tmax": (InputValue, {"dtype": float, "default": 200.0, "help": "..."}),
        "ou_nt": (InputValue, {"dtype": int, "default": 2000, "help": "..."}),
        "ou_print": (InputValue, {"dtype": bool, "default": True, "help": "..."}),
        "ou_propagator": (InputValue, {"dtype": str, "default": "exact", "help": "..."}),

        "sigma_key": (InputValue, {"dtype": str, "default": "sigma", "help": "..."}),
        "friction_key": (InputValue, {"dtype": str, "default": "friction", "help": "..."}),

        # IMPORTANT: XML-safe defaults (no None)
        "sigma_rep": (
            InputValue,
            {"dtype": str, "default": "", "help": "'' => auto"},
        ),
        "sigma_index": (
            InputValue,
            {"dtype": int, "default": 0, "help": "0 => auto"},
        ),

        # IMPORTANT: empty int array default (not ones)
        "friction_atoms": (
            InputArray,
            {
                "dtype": int,
                "default": input_default(factory=lambda n: np.zeros(n, dtype=int), args=(0,)),
                "help": "...",
            },
        ),

        "broadcast_single_bead_friction": (
            InputValue,
            {"dtype": bool, "default": True, "help": "..."},
        ),
    }

    default_help = "Friction operator configuration (MF + markovian/non-markovian bath)."
    default_label = "FRICTION"

    def store(self, friction: Friction) -> None:
        super(InputFriction, self).store(friction)

        if not isinstance(friction, Friction):
            return

        self.use_linear_coupling.store(friction.use_linear_coupling)
        self.bath_mode.store(friction.bath_mode)
        self.mf_mode.store(friction.mf_mode)

        self.spectral_density.store(friction.spectral_density)
        self.alpha_input.store(friction.alpha_input)
        self.friction_static.store(friction.friction_static)

        self.ou_fit_kind.store(friction.ou_fit_kind)
        self.ou_nterms.store(friction.ou_nterms)
        self.ou_tmax.store(friction.ou_tmax)
        self.ou_nt.store(friction.ou_nt)
        self.ou_print.store(friction.ou_print)
        self.ou_propagator.store(friction.ou_propagator)

        self.sigma_key.store(friction.sigma_key)
        self.friction_key.store(friction.friction_key)

        # additions: NEVER store None into str/int/array typed inputs
        rep = getattr(friction, "sigma_rep", "")
        if rep is None:
            rep = ""
        self.sigma_rep.store(rep)

        idx = getattr(friction, "sigma_index", 0)
        if idx is None:
            idx = 0
        self.sigma_index.store(int(idx))

        fa = getattr(friction, "friction_atoms", None)
        if fa is None:
            fa = np.zeros(0, dtype=int)
        self.friction_atoms.store(np.array(fa, dtype=int).flatten())

        self.broadcast_single_bead_friction.store(
            getattr(friction, "broadcast_single_bead_friction", True)
        )

    def fetch(self) -> Friction:
        # CRITICAL: do NOT pass None. Keep sentinels.
        return Friction(
            use_linear_coupling=self.use_linear_coupling.fetch(),
            bath_mode=self.bath_mode.fetch(),
            mf_mode=self.mf_mode.fetch(),
            spectral_density=self.spectral_density.fetch(),
            alpha_input=self.alpha_input.fetch(),
            friction_static=self.friction_static.fetch(),
            ou_fit_kind=self.ou_fit_kind.fetch(),
            ou_nterms=self.ou_nterms.fetch(),
            ou_tmax=self.ou_tmax.fetch(),
            ou_nt=self.ou_nt.fetch(),
            ou_print=self.ou_print.fetch(),
            ou_propagator=self.ou_propagator.fetch(),
            sigma_key=self.sigma_key.fetch(),
            friction_key=self.friction_key.fetch(),

            # additions
            sigma_rep=self.sigma_rep.fetch(),          # "" means auto
            sigma_index=self.sigma_index.fetch(),      # 0 means auto
            friction_atoms=self.friction_atoms.fetch(),
        )

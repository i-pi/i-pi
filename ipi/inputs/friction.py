import numpy as np

from ipi.engine.friction import Friction
from ipi.utils.inputvalue import Input, InputArray, InputValue, input_default


class InputFriction(Input):
    """
    Input handler for the Friction operator (used by dynamics mode '*-f').

    Markovian bath:
      - does NOT depend on spectral_density
      - general coupling: uses sigma=dF/dr from forces.extras and Gamma = sigma^T sigma
      - linear coupling: uses friction_static as Markovian gamma

    MF policy:
      - if bath_mode='markovian': MF is OFF unless alpha_input is provided
      - if bath_mode='non-markovian': alpha can come from alpha_input or spectral_density
    """

    attribs = {}

    fields = {
        "use_linear_coupling": (
            InputValue,
            {
                "dtype": bool,
                "default": False,
                "help": (
                    "If True, use linear coupling special case (sigma not needed for bath). "
                    "If False, use sigma=dF/dr from forces.extras for general coupling."
                ),
            },
        ),

        "bath_mode": (
            InputValue,
            {
                "dtype": str,
                "default": "non-markovian",
                "help": "One of: 'none', 'markovian', 'non-markovian'.",
            },
        ),

        "mf_mode": (
            InputValue,
            {
                "dtype": str,
                "default": "reconstruct",
                "help": (
                    "One of: 'none', 'linear', 'reconstruct'. "
                    "Note: in bath_mode='markovian', MF is disabled unless alpha_input is provided."
                ),
            },
        ),

        "spectral_density": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.ones, args=(0,)),
                "help": (
                    "Two-column data [omega, J(omega)] defining Lambda(omega)=J(omega)/omega. "
                    "Required for bath_mode='non-markovian' (OU fitting). "
                    "Optional for MF alpha if alpha_input not provided. "
                    "Not required for bath_mode='markovian'."
                ),
            },
        ),

        "alpha_input": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.ones, args=(0,)),
                "help": (
                    "Optional table [omega_k, alpha(omega_k)] with one row per RP normal mode. "
                    "If provided and omega_k matches nm.omegak, used for MF. "
                    "Required if you want MF in bath_mode='markovian'."
                ),
            },
        ),

        "friction_static": (
            InputValue,
            {
                "dtype": float,
                "default": 1.0,
                "help": (
                    "General coupling: scales coupling strength via sigma_eff = friction_static * sigma "
                    "(so friction and noise scale consistently). "
                    "Linear coupling + Markovian: interpreted as Markovian gamma."
                ),
            },
        ),

        # OU fitting (non-markovian)
        "ou_fit_kind": (
            InputValue,
            {
                "dtype": str,
                "default": "damped_cosine",
                "help": "Kernel fit family: 'exp' or 'damped_cosine'.",
            },
        ),
        "ou_nterms": (
            InputValue,
            {"dtype": int, "default": 4, "help": "Number of fit terms per normal mode."},
        ),
        "ou_tmax": (
            InputValue,
            {"dtype": float, "default": 200.0, "help": "Max time for K^(n)(t) fit grid."},
        ),
        "ou_nt": (
            InputValue,
            {"dtype": int, "default": 2000, "help": "Number of time samples for fitting."},
        ),
        "ou_print": (
            InputValue,
            {"dtype": bool, "default": True, "help": "Print fit parameters at bind()."},
        ),
        "ou_propagator": (
            InputValue,
            {
                "dtype": str,
                "default": "exact",
                "help": "OU auxiliary propagation: 'exact' (stable) or 'euler' (legacy/debug).",
            },
        ),

        # extras keys
        "sigma_key": (
            InputValue,
            {
                "dtype": str,
                "default": "sigma",
                "help": "forces.extras key for sigma = dF/dr (general coupling).",
            },
        ),
        "friction_key": (
            InputValue,
            {
                "dtype": str,
                "default": "friction",
                "help": "forces.extras key for friction tensor (currently unused in sigma-based modes).",
            },
        ),
    }

    default_help = "Friction operator configuration (MF + markovian/non-markovian bath)."
    default_label = "FRICTION"

    def store(self, friction: Friction) -> None:
        super(InputFriction, self).store(friction)
        if isinstance(friction, Friction):
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

    def fetch(self) -> Friction:
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
        )

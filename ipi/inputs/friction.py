from ipi.utils.array_backend import to_numpy
import json

import numpy as np

from ipi.engine.friction import Friction
from ipi.utils.inputvalue import Input, InputArray, InputValue, input_default


class InputFriction(Input):
    attribs = {}

    fields = {
        "forcefield": (
            InputValue,
            {
                "dtype": str,
                "default": "",
                "help": "Optional forcefield name for one whole-ring request returning "
                "the correction energy/forces and bath factors. Uses ordinary i-PI "
                "socket messages with bead-major flattened coordinates. This forcefield "
                "must not also appear in the potential-energy forces list.",
            },
        ),
        "s": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "dimension": "ms-momentum",
                "help": "Checkpointed non-Markovian auxiliary momenta, shape "
                "(nbeads,naux,nchannels), or (nbeads,naux,3*natoms) for scalar "
                "friction. Omit for a fresh thermal bath. Restoring these values "
                "preserves memory without drawing new random numbers.",
            },
        ),
        "ediss": (
            InputValue,
            {
                "dtype": float,
                "default": 0.0,
                "dimension": "energy",
                "help": "Checkpointed cumulative energy transferred to the friction bath.",
            },
        ),
        "bath_contract": (
            InputValue,
            {
                "dtype": str,
                "default": "",
                "help": "Checkpointed JSON signature of active Cartesian DOFs, "
                "bath channel labels and units. Validates channel identity on restart.",
            },
        ),
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
            {
                "dtype": str,
                "default": "non-markovian",
                "help": "Electronic-friction bath: none, markovian (one bead only), "
                "or non-markovian (requires fitted Ap matrices per normal mode).",
            },
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
                "help": "Physical friction spectrum in atomic units, shape "
                "(nfrequency,2) with rows [omega, Lambda(omega)]. Used to "
                "integrate MF coefficients independently of Ap. Frequencies "
                "must start at zero and be strictly increasing; cover the "
                "spectrum sufficiently densely and to a converged cutoff. "
                "Required with MF on unless debug_alpha_input is supplied.",
            },
        ),
        "Ap": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=((0, 0),)),
                "help": "Mode-resolved momentum-plus-auxiliary drift matrices "
                "for the Trenins--Rossi non-Markovian RPMD-GLE. For P beads "
                "use shape (P,1+naux,1+naux), with "
                "Ap[n]=[[0,theta_n^T],[-theta_n,A_n]]. Matrices must be in "
                "the native i-PI normal-mode order (centroid n=0 first). A "
                "2-D matrix is valid only for P=1. These matrices are fitted "
                "externally; i-PI does not fit or broadcast them. They specify "
                "dynamics only, not the MF coefficients; supply Lambda for MF.",
            },
        ),
        "debug_alpha_input": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "Explicit debug override for MF coefficients, shape "
                "(nbeads,2) with rows [normal-mode frequency, alpha]. Takes "
                "precedence over Lambda and must use i-PI normal-mode order.",
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
    }

    default_help = "Friction operator configuration (MF + Markovian/non-Markovian bath) using a canonical driver coupling Jacobian and Gamma."
    default_label = "FRICTION"

    def store(self, friction: Friction) -> None:
        super(InputFriction, self).store(friction)

        if not isinstance(friction, Friction):
            return

        self.forcefield.store(friction.forcefield)
        self.variable_friction.store(friction.variable_friction)
        self.bath_mode.store(friction.bath_mode)
        self.debug_mf_mode.store(friction.debug_mf_mode)

        self.Lambda.store(friction.Lambda)
        self.Ap.store(friction.Ap)
        self.debug_alpha_input.store(friction.debug_alpha_input)
        self.sigma_static.store(friction.sigma_static)
        bath_s = getattr(friction.bath, "s", None)
        self.s.store(friction._initial_s if bath_s is None else to_numpy(bath_s))
        self.ediss.store(friction.ediss)
        signature = friction._friction_contract_signature
        self.bath_contract.store("" if signature is None else json.dumps(signature))

        self.coupling_jacobian_key.store(friction.coupling_jacobian_key)
        self.gamma_key.store(friction.gamma_key)
        self.friction_meta_key.store(friction.friction_meta_key)
        self.coupling_key.store(friction.coupling_key)

    def fetch(self) -> Friction:
        return Friction(
            forcefield=self.forcefield.fetch(),
            variable_friction=self.variable_friction.fetch(),
            bath_mode=self.bath_mode.fetch(),
            debug_mf_mode=self.debug_mf_mode.fetch(),
            Lambda=self.Lambda.fetch(),
            Ap=self.Ap.fetch(),
            debug_alpha_input=self.debug_alpha_input.fetch(),
            sigma_static=self.sigma_static.fetch(),
            s=self.s.fetch(),
            ediss=self.ediss.fetch(),
            bath_contract=self.bath_contract.fetch(),
            coupling_jacobian_key=self.coupling_jacobian_key.fetch(),
            gamma_key=self.gamma_key.fetch(),
            friction_meta_key=self.friction_meta_key.fetch(),
            coupling_key=self.coupling_key.fetch(),
        )

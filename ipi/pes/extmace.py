"""MACE client driver with fixed-electric-field coupling.

The electric field is supplied by i-PI in the request extra data.  The
calculator adds ``-mu.E`` to the MACE energy before differentiating, so the
returned forces and stress include the field coupling.  This module
intentionally implements only fixed-E client coupling.
"""

import json

import numpy as np

from ipi.pes._mace import BatchedMACE, MACE_driver
from ipi.utils.units import unit_to_user

__DRIVER_NAME__ = "extmace"
__DRIVER_CLASS__ = "Extended_MACE_driver"


class ExtendedMACECalculator(BatchedMACE):
    """Batched MACE calculator that applies an electric field from i-PI."""

    def __init__(self, *args, **kwargs):
        self.extras = {}
        super().__init__(*args, **kwargs)

    def _electric_field(self, reference):
        if "electric_field" not in self.extras:
            raise ValueError(
                "The request extra data must contain 'electric_field' for extmace."
            )
        field = unit_to_user(
            "electric-field", "v/ang", np.asarray(self.extras["electric_field"])
        )
        field = np.asarray(field, dtype=float)
        if field.shape != (3,):
            raise ValueError(
                f"'electric_field' must have shape (3,), got {field.shape}."
            )
        # MACE dipoles are in e Angstrom and energies in eV; V/Angstrom gives
        # the corresponding eV coupling directly.
        import torch

        return torch.as_tensor(field, device=reference.device, dtype=reference.dtype)

    def augment_output(self, data, batch, training):
        if "electric_field" in self.extras:
            field = self._electric_field(data["energy"])
            if field.any():
                dipole = self._response_dipole(data)
                data["energy"] = data["energy"] - dipole @ field
        return super().augment_output(data, batch, training)


class Extended_MACE_driver(MACE_driver):
    """MACE Python driver exposing fixed-E client coupling."""

    requires_extra = True

    def check_parameters(self):
        self.batched_calculator = ExtendedMACECalculator(
            model_paths=self.model, device=self.device, **self.mace_kwargs
        )

    def store_extra(self, extra):
        """Receive the per-request JSON payload used by ``FFDielectric``."""
        if not extra:
            self.extra = {}
        elif isinstance(extra, str):
            self.extra = json.loads(extra)
        else:
            self.extra = dict(extra)

    def compute(self, cell, pos):
        self.batched_calculator.extras = getattr(self, "extra", {}) or {}
        return super().compute(cell, pos)

    def post_process(self, properties, structure):
        result = list(super().post_process(properties, structure))
        if "electric_field" in (getattr(self, "extra", {}) or {}):
            extras = {} if not result[3] else json.loads(result[3])
            extras["applied_fields"] = ["electric_field"]
            result[3] = json.dumps(extras)
        return tuple(result)

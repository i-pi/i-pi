#!/usr/bin/env python
"""MACE calculator extensions for electric fields and dielectric response.

The standard batching, model evaluation, force/stress calculation, output
handling, and command-line interface live in :mod:`ipi.pes._mace`.  This module
only adds electric-field coupling and the associated response tensors.

Runnable static-field and resonant-field examples for both client-side and
server-side coupling are provided under ``examples/features/ffdieletric``.
"""

import json
from typing import Dict, Optional, Tuple

import numpy as np
import torch

from ipi.pes._mace import (
    BatchedMACE,
    MACE_driver,
    ase_like_properties as mace_ase_like_properties,
    compute_dielectric_gradients,
    proper_dipole,
    run_cli,
)
from ipi.utils.units import unit_to_user

__DRIVER_NAME__ = "extmace"
__DRIVER_CLASS__ = "Extended_MACE_driver"

__all__ = [
    "Extended_MACE_driver",
    "ExtendedMACECalculator",
    "add_bec_inplace",
    "compute_dielectric_gradients",
    "proper_dipole",
]

_EXTENDED_ASE_LIKE_PROPERTIES = {
    "node_energy": ("natoms",),
    "atomic-oxn-dipole": ("natoms", 3),
    "BEC": ("natoms", 3, 3),
    "BECx": ("natoms", 3),
    "BECy": ("natoms", 3),
    "BECz": ("natoms", 3),
}

# Keep these public module-level names for backwards compatibility.
ase_like_properties = mace_ase_like_properties.copy()
ase_like_properties.update(_EXTENDED_ASE_LIKE_PROPERTIES)
to_ignore_properties = ["interaction_energy", "node_feats"]


def add_bec_inplace(data: Dict[str, torch.Tensor], bec: torch.Tensor) -> None:
    """Store a 3xNx3 Born-charge tensor as ASE-compatible per-atom arrays."""

    data["BECx"] = bec[0, :, :]
    data["BECy"] = bec[1, :, :]
    data["BECz"] = bec[2, :, :]


class Extended_MACE_driver(MACE_driver):
    """In-process MACE driver with electric-field response support."""

    def check_parameters(self):
        """Initialize the extended calculator using extras sent by i-PI."""

        self.batched_calculator = ExtendedMACECalculator(
            model_paths=self.model,
            device=self.device,
            **self.mace_kwargs,
        )

    def compute(self, cell, pos):
        """Evaluate using the extra information attached to this request."""

        self.batched_calculator.extras = self.extra or {}
        return super().compute(cell, pos)

    def post_process(self, properties, structure):
        """Return normal MACE extras and acknowledge any applied field."""
        energy, forces, virial, extras = super().post_process(properties, structure)
        if "Efield" not in (self.extra or {}):
            return energy, forces, virial, extras

        extras_dict = {} if not extras else json.loads(extras)
        extras_dict["applied_fields"] = ["electric_field"]
        return energy, forces, virial, json.dumps(extras_dict)


class ExtendedMACECalculator(BatchedMACE):
    """Batched MACE calculator extended with electric-field coupling."""

    ignored_properties = frozenset(to_ignore_properties)

    def __init__(
        self,
        instructions: Optional[dict] = None,
        ase_like_properties: Optional[Dict[str, Tuple]] = None,
        *args,
        **kwargs,
    ):
        # BatchedMACE normalizes and consumes some instruction values, so copy
        # nested mutable state before adding the extended settings.
        instructions = {} if instructions is None else instructions.copy()
        instructions["forward_kwargs"] = instructions.get("forward_kwargs", {}).copy()

        if "ensemble" in instructions:
            raise ValueError(
                "The extmace 'ensemble' instruction has been removed. "
                "Use <electric_field> in <ffdielectric> and select where the "
                "coupling is applied with its 'where' attribute."
            )

        self.extras = {}

        properties = _EXTENDED_ASE_LIKE_PROPERTIES.copy()
        if ase_like_properties is not None:
            properties.update(ase_like_properties)

        super().__init__(instructions, properties, *args, **kwargs)

    def augment_output(
        self,
        data: Dict[str, torch.Tensor],
        batch: Dict[str, torch.Tensor],
        training: bool,
        compute_bec: bool,
    ) -> Dict[str, torch.Tensor]:
        """Add electric-field contributions and dielectric response tensors."""

        if "Dfield" in self.extras:
            raise NotImplementedError(
                "Electric-displacement coupling is not implemented in extmace yet."
            )

        if "Efield" in self.extras:
            mu = self._response_dipole(data)
            electric_field = self._electric_field(mu)
            # Differentiating the field-coupled energy supplies the field
            # contributions to forces and stress automatically.
            data["energy"] -= mu @ electric_field
            data = self.get_forces_stress(data, batch, training)
        else:
            data = self.get_forces_stress(data, batch, training)

        if compute_bec:
            bec, _ = self.add_dielectric_response(data, batch)
            add_bec_inplace(data, bec)

        # Preserve the extended driver's historical per-atom interaction
        # energies while the base MACE driver avoids transferring this output.
        if data.get("node_energy") is not None and "node_e0" in batch:
            data["node_energy"] = data["node_energy"] - batch["node_e0"]

        return data

    def _electric_field(self, reference: torch.Tensor) -> torch.Tensor:
        extras = self.extras
        if not extras or "Efield" not in extras:
            raise ValueError(
                "The extra information dictionary must contain 'Efield' for "
                "client-side electric-field coupling."
            )

        electric_field = np.asarray(extras["Efield"])
        electric_field = unit_to_user("electric-field", "v/ang", electric_field)
        electric_field = torch.as_tensor(
            electric_field,
            device=reference.device,
            dtype=reference.dtype,
        )
        if electric_field.shape != (3,):
            raise ValueError(
                f"'Efield' must have shape (3,), got {tuple(electric_field.shape)}."
            )
        return electric_field


if __name__ == "__main__":
    run_cli(
        calculator_class=ExtendedMACECalculator,
        calculator_name="ExtendedMACECalculator",
    )

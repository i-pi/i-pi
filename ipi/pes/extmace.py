#!/usr/bin/env python
"""MACE calculator extensions for electric fields and dielectric response.

The standard batching, model evaluation, force/stress calculation, output
handling, and command-line interface live in :mod:`ipi.pes._mace`.  This module
only adds the electric-field ensembles and the associated response tensors.
"""

from typing import Dict, List, Optional, Tuple

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
from ipi.pes.tools import JSONLogger, Parent
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


class ExtendedMACECalculator(BatchedMACE):
    """Batched MACE calculator extended with electric-field ensembles."""

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

        self.ensemble = str(instructions.get("ensemble", "none")).upper()
        instructions["ensemble"] = self.ensemble
        if self.ensemble not in {"NONE", "E"}:
            raise ValueError(f"Ensemble {self.ensemble} not implemented (yet).")

        instructions.pop("log", None)
        self.results_logger = JSONLogger(instructions.pop("log_results", None))
        self.extras = {}

        properties = _EXTENDED_ASE_LIKE_PROPERTIES.copy()
        if ase_like_properties is not None:
            properties.update(ase_like_properties)

        super().__init__(instructions, properties, *args, **kwargs)

    def compute_batched(self, atoms) -> List[Parent]:
        """Evaluate structures and optionally retain the legacy result log."""

        results = super().compute_batched(atoms)
        for index, result in enumerate(results):
            self.results_logger.save(result, f"results.{index}.json")
        return results

    def augment_output(
        self,
        data: Dict[str, torch.Tensor],
        batch: Dict[str, torch.Tensor],
        training: bool,
        compute_bec: bool,
    ) -> Dict[str, torch.Tensor]:
        """Add electric-field contributions and dielectric response tensors."""

        if self.ensemble == "NONE":
            data = self.get_forces_stress(data, batch, training)
        else:
            mu = self._proper_model_dipole(data)
            electric_field = self._electric_field(mu)
            # Differentiating the field-coupled energy supplies the field
            # contributions to forces and stress automatically.
            data["energy"] -= mu @ electric_field
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
                "The extra information dictionary must contain 'Efield' when "
                f"using ensemble '{self.ensemble}'."
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

#!/usr/bin/env python
"""MACE calculator extensions for electric fields and dielectric response.

The standard batching, model evaluation, force/stress calculation, output
handling, and command-line interface live in :mod:`ipi.pes._mace`.  This module
only adds electric-field coupling and the associated response tensors.

Runnable static-field and resonant-field examples for both client-side and
server-side coupling are provided under ``examples/features/ffdieletric``.
"""

import json
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch
from ase.outputs import _defineprop, all_outputs

from ipi.pes._mace import (
    BatchedMACE,
    MACE_driver,
    ase_like_properties as mace_ase_like_properties,
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
    "piezoelectric": (3, 3, 3),
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


def _normalize_dipole_shape(mu: torch.Tensor, strain: torch.Tensor) -> torch.Tensor:
    """Normalize removable singleton dimensions in a batched model dipole."""

    if not isinstance(mu, torch.Tensor):
        raise TypeError(f"The MACE dipole must be a torch.Tensor, got {type(mu)}.")
    if not isinstance(strain, torch.Tensor):
        raise TypeError(
            f"The MACE displacement must be a torch.Tensor, got {type(strain)}."
        )
    if strain.ndim < 2 or tuple(strain.shape[-2:]) != (3, 3):
        raise ValueError(
            "The MACE displacement must end in shape (3, 3), got "
            f"{tuple(strain.shape)}."
        )

    expected_shape = (*strain.shape[:-2], 3)
    if tuple(mu.shape) == expected_shape:
        return mu

    # Some MACE/PyTorch combinations wrap per-structure outputs in an extra
    # singleton dimension, e.g. (B, 1, 3) rather than (B, 3). Reshaping is safe
    # only when removing singleton dimensions recovers the expected layout.
    observed_non_singleton = tuple(size for size in mu.shape if size != 1)
    expected_non_singleton = tuple(size for size in expected_shape if size != 1)
    if observed_non_singleton != expected_non_singleton:
        raise ValueError(
            "The MACE dipole shape is incompatible with its displacement: "
            f"expected {expected_shape}, got {tuple(mu.shape)}. Only extra "
            "singleton dimensions can be normalized."
        )

    return mu.reshape(expected_shape)


def proper_dipole(mu: torch.Tensor, strain: torch.Tensor) -> torch.Tensor:
    """Return the dipole corrected for an infinitesimal symmetric strain."""

    mu = _normalize_dipole_shape(mu, strain)
    symmetric_strain = 0.5 * (strain + strain.transpose(-1, -2))
    return mu - torch.einsum("...il,...l->...i", symmetric_strain, mu)


# Function adapted from https://github.com/davkovacs/mace/tree/mu_alpha
def compute_dielectric_gradients(
    dielectric: torch.Tensor,
    inputs: List[torch.Tensor],
    clean: Optional[bool] = False,
) -> List[torch.Tensor]:
    """Differentiate a dielectric tensor with respect to input tensors."""

    d_dielectric_dr = [
        [None for _ in range(dielectric.shape[-1])] for _ in range(len(inputs))
    ]
    grad_outputs: List[torch.Tensor] = [
        torch.ones((dielectric.shape[0], 1)).to(dielectric.device)
    ]
    for i in range(dielectric.shape[-1]):
        gradients = torch.autograd.grad(
            outputs=[dielectric[:, i].unsqueeze(-1)],
            inputs=inputs,
            grad_outputs=grad_outputs,
            retain_graph=(i < dielectric.shape[-1] - 1) or not clean,
            create_graph=False,
            allow_unused=False,
        )
        assert len(gradients) == len(inputs), "coding error"
        for j, (gradient, input_tensor) in enumerate(zip(gradients, inputs)):
            assert gradient.shape == input_tensor.shape, "coding error"
            d_dielectric_dr[j][i] = gradient.detach()
        del gradients
    del grad_outputs
    return [torch.stack(out, dim=0) for out in d_dielectric_dr]


class Extended_MACE_driver(MACE_driver):
    """In-process MACE driver with electric-field response support."""

    def __init__(self, *args, use_proper_dipole=None, **kwargs):
        """Initialize the MACE driver and its electrical-response options."""

        self._use_proper_dipole = use_proper_dipole
        super().__init__(*args, **kwargs)

    def check_parameters(self):
        """Initialize the extended calculator using extras sent by i-PI."""

        if self._use_proper_dipole is not None:
            self.mace_kwargs["use_proper_dipole"] = self._use_proper_dipole

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
    supported_instruction_keys = BatchedMACE.supported_instruction_keys | {
        "compute_BEC"
    }

    def __init__(
        self,
        instructions: Optional[dict] = None,
        ase_like_properties: Optional[Dict[str, Tuple]] = None,
        use_proper_dipole: bool = True,
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

        compute_bec = instructions.get("compute_BEC", False)
        if not isinstance(compute_bec, bool):
            raise TypeError("'instructions.compute_BEC' must be a boolean")
        if not isinstance(use_proper_dipole, bool):
            raise TypeError("'use_proper_dipole' must be a boolean")

        self.extras = {}

        properties = _EXTENDED_ASE_LIKE_PROPERTIES.copy()
        if ase_like_properties is not None:
            properties.update(ase_like_properties)

        super().__init__(instructions, properties, *args, **kwargs)

        # MACECalculator owns a lower-case ``compute_bec`` attribute, so keep
        # the i-PI response switch distinct from that upstream implementation.
        self.compute_bec_response = compute_bec
        self.use_proper_dipole = use_proper_dipole
        if self.compute_bec_response and "BEC" not in all_outputs:
            _defineprop("BEC", dtype=float, shape=("natoms", 3, 3))

    def requires_model_gradients(self) -> bool:
        """Retain the graph when electrical response derivatives are requested."""

        return self.compute_bec_response

    def _response_dipole(self, data: Dict[str, torch.Tensor]) -> torch.Tensor:
        """Return the selected model dipole for electrical response."""

        if data.get("dipole") is None:
            raise ValueError(
                "The selected MACE model does not provide the dipole required "
                "for electrical-response calculations."
            )
        if data.get("displacement") is None:
            raise ValueError(
                "The MACE output does not contain the displacement tensor "
                "required for electrical-response calculations."
            )

        mu = _normalize_dipole_shape(data["dipole"], data["displacement"])
        # Store the normalized value too, so a model compatibility dimension
        # cannot reach ModelResults when the raw dipole is also requested.
        data["dipole"] = mu
        if not self.use_proper_dipole:
            return mu
        return proper_dipole(mu, data["displacement"])

    def add_dielectric_response(
        self,
        data: Dict[str, torch.Tensor],
        batch: Dict[str, torch.Tensor],
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """Compute and store Born charges and the piezoelectric tensor."""

        bec, dmu_deta = self.compute_dmu_dR_deta(data, batch)

        if data.get("BEC") is None:
            # (mu_xyz,node,R_xyz) --> (node,mu_xyz,R_xyz)
            data["BEC"] = bec.moveaxis(0, 1)

        if data.get("piezoelectric") is None:
            # (mu_xyz,graph,eta_i,eta_j) --> (graph,mu_xyz,eta_i,eta_j)
            cell = batch.get("cell")
            if not isinstance(cell, torch.Tensor):
                raise ValueError("The MACE batch does not contain a tensor cell.")
            volume = torch.linalg.det(cell.view(-1, 3, 3)).abs()
            if volume.shape[0] != dmu_deta.shape[1]:
                raise ValueError(
                    "The number of cell volumes does not match the number of "
                    "dipole-strain derivatives."
                )
            data["piezoelectric"] = (
                dmu_deta.moveaxis(0, 1) / volume[:, None, None, None]
            )

        return bec, dmu_deta

    def compute_dmu_dR_deta(
        self,
        data: Dict[str, torch.Tensor],
        batch: Dict[str, torch.Tensor],
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """Differentiate the dipole with respect to positions and strain."""

        mu = self._response_dipole(data)
        positions = batch.get("positions")
        if not isinstance(positions, torch.Tensor):
            raise ValueError("The MACE batch does not contain tensor positions.")
        if not positions.requires_grad:
            raise ValueError(
                "MACE positions must retain gradients when 'compute_BEC: true' "
                "so that Born charges can be computed."
            )

        displacement = data["displacement"]
        if not isinstance(displacement, torch.Tensor):
            raise ValueError("The MACE output does not contain a tensor displacement.")
        if not displacement.requires_grad:
            raise ValueError(
                "The MACE displacement tensor must retain gradients when "
                "'compute_BEC: true' so that the piezoelectric tensor can be "
                "computed."
            )

        bec, dmu_deta = compute_dielectric_gradients(
            mu,
            [positions, displacement],
        )

        expected_bec_shape = (3, *positions.shape)
        if tuple(bec.shape) != expected_bec_shape:
            raise ValueError(
                "The computed Born charges have the wrong shape: expected "
                f"{expected_bec_shape}, got {tuple(bec.shape)}."
            )

        expected_strain_shape = (3, *displacement.shape)
        if tuple(dmu_deta.shape) != expected_strain_shape:
            raise ValueError(
                "The dipole-strain derivative has the wrong shape: expected "
                f"{expected_strain_shape}, got {tuple(dmu_deta.shape)}."
            )
        if not torch.allclose(dmu_deta, dmu_deta.transpose(-1, -2)):
            nonsymmetric_norm = torch.linalg.norm(
                dmu_deta - dmu_deta.transpose(-1, -2)
            ).item()
            raise ValueError(
                "The dipole-strain derivative is not symmetric: "
                f"the norm of its nonsymmetric part is {nonsymmetric_norm:.6e}."
            )

        return bec, dmu_deta

    def compute_dmu_dR(
        self,
        data: Dict[str, torch.Tensor],
        batch: Dict[str, torch.Tensor],
    ) -> torch.Tensor:
        """Compute Born effective charges, retaining the historical API."""

        return self.compute_dmu_dR_deta(data, batch)[0]

    def augment_output(
        self,
        data: Dict[str, torch.Tensor],
        batch: Dict[str, torch.Tensor],
        training: bool,
    ) -> Dict[str, torch.Tensor]:
        """Add electric-field contributions and dielectric response tensors."""

        if "Dfield" in self.extras:
            raise NotImplementedError(
                "Electric-displacement coupling is not implemented in extmace yet."
            )

        if "Efield" in self.extras:
            electric_field = self._electric_field(data["energy"])
            mu = self._response_dipole(data)
            # Differentiating the field-coupled energy supplies the field
            # contributions to forces and stress automatically.
            data["energy"] -= mu @ electric_field
            data = self.get_forces_stress(data, batch, training)
        else:
            data = self.get_forces_stress(data, batch, training)

        if self.compute_bec_response:
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

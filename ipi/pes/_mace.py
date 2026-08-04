#!/usr/bin/env python
"""An interface for the [MACE](https://github.com/ACEsuit/mace) calculator that
supports batched evaluation.

See ``examples/clients/mace-batched`` for a complete example.
"""

try:
    import mace  # noqa: F401

except Exception:
    msg = (
        "You are trying to use the MACE client of i-PI, "
        "but it seems that you have not installed MACE.\n"
        "Please install MACE via:\n"
        "    pip install mace-torch\n"
        "For more information, visit:\n"
        "    https://github.com/ACEsuit/mace"
    )

    print(msg)
    raise

import argparse
import json
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch
from ase import Atoms
from ase.io import read, write
from mace import data
from mace.calculators import MACECalculator
from mace.modules.utils import get_outputs
from mace.tools import torch_tools
from mace.tools.torch_geometric.batch import Batch
from mace.tools.torch_geometric.dataloader import DataLoader

from ipi.pes._ase import ASEDriver
from ipi.pes.tools import ModelResults, Parent

# --------------------------------------- #
__DRIVER_NAME__ = "mace"
__DRIVER_CLASS__ = "MACE_driver"

MAX_VOLUME = 1e12

_DEFAULT_ASE_LIKE_PROPERTIES = {
    "energy": (),
    "interaction_energy": (),
    "forces": ("natoms", 3),
    "displacement": (3, 3),
    "stress": (3, 3),
    "virials": (3, 3),
    "dipole": (3,),
    "charges": ("natoms",),
    "atomic_dipoles": ("natoms", 3),
    "polarizability": (3, 3),
    "polarizability_sh": (6,),
}
# Kept as a public module-level name for backwards compatibility.
ase_like_properties = _DEFAULT_ASE_LIKE_PROPERTIES

to_ignore_properties = ["interaction_energy", "node_feats", "node_energy"]


class MACE_driver(ASEDriver):
    """ASE driver for running MACE models with batched torch-based execution."""

    template: Atoms

    def __init__(
        self,
        template,
        model,
        device="cpu",
        mace_kwargs=None,
        *args,
        **kwargs,
    ):
        """
        Initialize the MACE driver.

        Parameters
        ----------
        template : str or Atoms
            Structure template.
        model : str or list[str]
            Path(s) to MACE model(s).
        device : str
            Torch device ("cpu" or "cuda").
        mace_kwargs : str or None
            Path to JSON file with MACE kwargs.
        """

        self.model = model
        self.device = device
        self.mace_kwargs = {}
        self.all_templates = None

        if mace_kwargs is not None:
            with open(mace_kwargs, "r") as f:
                self.mace_kwargs = json.load(f)
        template = read(template)
        super().__init__(template, *args, **kwargs)

    def check_parameters(self):
        """Initialize the MACE calculator from the provided model and settings."""

        self.batched_calculator = BatchedMACE(
            model_paths=self.model,
            device=self.device,
            **self.mace_kwargs,
        )

    def template2atoms(
        self,
        cell: List[np.ndarray],
        pos: List[np.ndarray],
    ) -> List[Atoms]:
        """
        Build ASE Atoms objects from batched cells and positions.
        """

        Nstructures = len(cell)

        if self.all_templates is None:
            self.all_templates = [self.template.copy() for _ in range(Nstructures)]
        elif len(self.all_templates) < Nstructures:
            Nalready = len(self.all_templates)
            self.all_templates.extend(
                self.template.copy() for _ in range(Nstructures - Nalready)
            )

        for _, (atoms, c, p) in enumerate(zip(self.all_templates, cell, pos)):
            atoms.set_positions(p)
            atoms.set_pbc(True)
            atoms.set_cell(c)
            volume = atoms.get_volume()
            if volume > MAX_VOLUME:
                raise ValueError(
                    f"Unphysical structure volume detected: {volume:.3f} Å³"
                )

        return self.all_templates[:Nstructures]

    def compute(self, cell, pos):
        """
        Run MACE on one or more structures and return post-processed results.
        """

        if isinstance(cell, list):
            # convert from atomic_unit to angstrom
            for n, (c, p) in enumerate(zip(cell, pos)):
                cell[n], pos[n] = self.convert_units(c, p)

            # modify cell and positions, keep the other arrays and info as in the template
            atoms = self.template2atoms(cell, pos)
            results = self.batched_calculator.compute_batched(atoms)  # Dict[str,List]

            # convert from angstrom to atomic_unit
            out = [self.post_process(r, a) for r, a in zip(results, atoms)]

            return out[0] if len(out) == 1 else out
        else:
            return self.compute([cell], [pos])


class BatchedMACE(MACECalculator):
    """
    Batched evaluation for MACE:
     - supports batched evaluation of many atomic structures
     - can be used with a i-PI driver or as a standalone
    """

    ignored_properties = frozenset(to_ignore_properties)
    supported_instruction_keys = frozenset(
        {"_comments", "batch_size", "forward_kwargs", "ignore"}
    )

    def __init__(
        self,
        instructions: dict = None,
        ase_like_properties: Optional[Dict[str, Tuple]] = None,
        *argc,
        **kwargs,
    ):
        self.instructions = instructions if instructions is not None else {}
        unsupported = set(self.instructions) - self.supported_instruction_keys
        if unsupported:
            names = ", ".join(sorted(unsupported))
            supported = ", ".join(sorted(self.supported_instruction_keys))
            raise ValueError(
                f"Unsupported {__DRIVER_NAME__} instruction(s): {names}. "
                f"Supported instructions are: {supported}."
            )
        if "forward_kwargs" not in self.instructions:
            self.instructions["forward_kwargs"] = {}

        self.batch_size = self.instructions.pop("batch_size", None)
        self.ase_like_properties = _DEFAULT_ASE_LIKE_PROPERTIES.copy()
        if ase_like_properties is not None:
            self.ase_like_properties.update(
                self._normalize_ase_like_properties(ase_like_properties)
            )
        if "arrays_keys" not in kwargs:
            kwargs["arrays_keys"] = {}
        if "oxn" not in kwargs["arrays_keys"]:
            kwargs["arrays_keys"].update({"oxn": "oxn"})
        super().__init__(*argc, **kwargs)
        if self.model_type == "PolarMACE":
            # PolarMACE exposes additional per-structure and per-atom tensors.
            # Register their shapes so callers can retain them instead of being
            # forced to list every PolarMACE-specific output under ``ignore``.
            density_dim = (
                getattr(self.models[0], "atomic_multipoles_max_l", 0) + 1
            ) ** 2
            self.ase_like_properties.update(
                {
                    "electrostatic_energy": (),
                    "electron_energy": (),
                    "spins": ("natoms",),
                    "density_coefficients": ("natoms", density_dim),
                    "spin_charge_density": ("natoms", 2, density_dim),
                }
            )
        # MACE versions before 0.3.16 select the model dtype but do not expose
        # it on the calculator. Infer it from the (possibly converted) model so
        # batched preprocessing can use the same dtype on all supported versions.
        if not hasattr(self, "default_dtype"):
            self.default_dtype = next(self.models[0].parameters()).dtype
        self._output_summary_printed = False
        assert not self.use_compile, "self.use_compile=True is not supported yet."

    @staticmethod
    def _normalize_ase_like_properties(
        properties: Dict[str, Tuple],
    ) -> Dict[str, Tuple]:
        """Validate property shapes and convert JSON lists to tuples."""
        if not isinstance(properties, dict):
            raise TypeError("'ase_like_properties' must be a dictionary")

        normalized = {}
        for key, shape in properties.items():
            if not isinstance(key, str):
                raise TypeError(
                    "Property names in 'ase_like_properties' must be strings"
                )
            if not isinstance(shape, (list, tuple)):
                raise TypeError(
                    f"The shape of property '{key}' must be a list or tuple"
                )

            shape = tuple(shape)
            if any(
                dimension != "natoms"
                and (not isinstance(dimension, int) or dimension < 0)
                for dimension in shape
            ):
                raise ValueError(
                    f"Invalid shape for property '{key}': {shape}. "
                    "Dimensions must be non-negative integers or 'natoms'."
                )
            if "natoms" in shape and shape[0] != "natoms":
                raise ValueError(
                    f"Invalid shape for property '{key}': 'natoms' must be the "
                    "first dimension"
                )
            normalized[key] = shape

        return normalized

    def requires_model_gradients(self) -> bool:
        """Whether output augmentation needs the model autograd graph."""

        return False

    def preprocess(self, atoms: List[Atoms]) -> Tuple[DataLoader, Dict[str, bool]]:
        """
        Preprocess the calculation: prepare the batch, result tensors, etc.
        """

        keyspec = data.KeySpecification(
            info_keys=self.info_keys, arrays_keys=self.arrays_keys
        )
        with torch_tools.default_dtype(self.default_dtype):
            configs = data.config_from_atoms_list(
                atoms, key_specification=keyspec, head_name=self.head
            )
            dataset = [
                data.AtomicData.from_config(
                    config,
                    z_table=self.z_table,
                    cutoff=self.r_max,
                    heads=self.available_heads,
                ).to(self.device)
                for config in configs
            ]

        training = self.use_compile or self.requires_model_gradients()

        if self.model_type in ["MACE", "EnergyDipoleMACE", "PolarMACE"]:
            compute_stress = not self.use_compile
        else:
            compute_stress = False

        assert compute_stress, "'compute_stress' is False"

        # -------------------#
        # Some extra parameters to the model
        if "compute_edge_forces" not in self.instructions["forward_kwargs"]:
            self.instructions["forward_kwargs"][
                "compute_edge_forces"
            ] = self.compute_atomic_stresses
        if "compute_stress" not in self.instructions["forward_kwargs"]:
            self.instructions["forward_kwargs"]["compute_stress"] = compute_stress

        forward_kwargs = self.instructions["forward_kwargs"].copy()

        forward_kwargs["compute_force"] = False
        forward_kwargs["compute_stress"] = False
        forward_kwargs["compute_displacement"] = compute_stress
        forward_kwargs["compute_virials"] = False
        forward_kwargs["compute_edge_forces"] = False

        batch_size = (
            self.batch_size if self.batch_size is not None else max(len(dataset), 1)
        )
        data_loader = DataLoader(
            dataset,
            batch_size=batch_size,
            shuffle=False,
            drop_last=False,
        )

        return (
            data_loader,
            {
                "training": training,
                "forward_kwargs": forward_kwargs,
            },
        )

    @staticmethod
    def batch2natoms(batch: Batch) -> List[int]:
        """
        Return the number of atoms in batch.
        """
        return [
            a.shape[0]
            for a in np.split(
                batch["positions"].cpu().numpy(),
                batch["ptr"][1:].cpu().numpy(),
                axis=0,
            )[:-1]
        ]

    def compute_batched(self, atoms: List[Atoms]) -> List[Parent]:
        """
        Evaluate the model(s) on a list of structures.
        """

        data_loader, options = self.preprocess(atoms)
        training = options["training"]
        forward_kwargs = options["forward_kwargs"]

        # model evaluation
        model_results = [
            ModelResults(self.ase_like_properties) for _ in range(self.num_models)
        ]
        # loop over models in the committee
        for batch_base in data_loader:
            # batch = batch_base.to(self.device)
            batch = self._clone_batch(batch_base).to_dict()
            Natoms = self.batch2natoms(batch)

            if self.model_type in ["MACE", "EnergyDipoleMACE", "PolarMACE"]:
                node_heads = batch["head"][batch["batch"]]
                num_atoms_arange = torch.arange(
                    batch["positions"].shape[0], device=batch["positions"].device
                )

                # this try-except is to be compatible with different MACE versions
                try:
                    # newer versions of MACE
                    node_e0 = self.models[0].atomic_energies_fn(batch["node_attrs"])[
                        num_atoms_arange, node_heads
                    ]
                except Exception:
                    # older versions of MACE
                    node_e0 = self.models[0].atomic_energies_fn(
                        batch["node_attrs"], node_heads
                    )
                batch["node_e0"] = node_e0

            for i, model in enumerate(self.models):
                out: dict[str, torch.Tensor] = model(
                    batch,
                    training=training,
                    **forward_kwargs,
                )

                out = self.augment_output(out, batch, training)

                # collect the results
                ignored = set(self.ignored_properties)
                if "ignore" in self.instructions:
                    ignored |= set(self.instructions["ignore"])

                results_tensors = {}

                for key, value in out.items():
                    if value is None or key in ignored:
                        continue

                    results_tensors[key] = value.detach().cpu().numpy()

                if not self._output_summary_printed:
                    self._print_output_summary(out, results_tensors)
                    self._output_summary_printed = True

                model_results[i].store(Natoms, results_tensors)

        # re-order results
        return ModelResults.mean(model_results)

    def _print_output_summary(
        self,
        model_output: Dict[str, torch.Tensor],
        cpu_output: Dict[str, np.ndarray],
    ) -> None:
        """Print the available and CPU-transferred outputs once per calculator."""

        def classify(keys):
            arrays = []
            info = []
            unregistered = []
            for key in sorted(keys):
                shape = self.ase_like_properties.get(key)
                if shape is None:
                    unregistered.append(key)
                elif "natoms" in shape:
                    arrays.append(key)
                else:
                    info.append(key)
            return arrays, info, unregistered

        available_keys = [
            key for key, value in model_output.items() if value is not None
        ]
        model_arrays, model_info, model_unregistered = classify(available_keys)
        cpu_arrays, cpu_info, cpu_unregistered = classify(cpu_output)

        def display(keys):
            return ", ".join(keys) if keys else "(none)"

        def display_with_shapes(keys, values):
            return (
                ", ".join(
                    f"{key} (shape {tuple(np.asarray(values[key]).shape)})"
                    for key in keys
                )
                if keys
                else "(none)"
            )

        print("-----------------------------------")
        print("  MACE output summary (printed once):")
        print(f"  Using device: {self.device}")
        print(
            "  Performance tip: avoid unnecessary GPU-to-CPU transfers by "
            "ignoring optional outputs you do not need."
        )
        print("  To do so, in a MACE settings JSON file, use:")
        print('    {"instructions": {"ignore": ["property_name", "..."]}}')
        print(
            "  Replace 'property_name' with an optional name from the " "lists below."
        )
        print("  Never ignore energy, forces, or stress; i-PI " "requires them.")
        print(
            "  You can pass that JSON file as 'mace_kwargs' in the i-PI force-field "
            "parameters, or with --mace_kwargs in the standalone CLI."
        )
        print()
        print(f"  Produced per-atom properties (ASE arrays): {display(model_arrays)}")
        print(
            "  Produced per-structure properties (ASE info): " f"{display(model_info)}"
        )
        print(
            "  Produced unregistered model outputs: " f"{display(model_unregistered)}"
        )
        print(
            f"  Copied to CPU per-atom properties (ASE arrays): {display(cpu_arrays)}"
        )
        print(
            "  Copied to CPU per-structure properties (ASE info): "
            f"{display(cpu_info)}"
        )
        print(
            "  Copied to CPU unregistered model outputs: "
            f"{display_with_shapes(cpu_unregistered, cpu_output)}"
        )
        if cpu_unregistered:
            print(
                "  Register an output under 'ase_like_properties' in the MACE "
                "settings JSON, or skip it with 'instructions.ignore'. If i-PI "
                "cannot infer the registration automatically, the error below "
                "will show the observed shape and suggested JSON."
            )
        print("-----------------------------------")

    def augment_output(
        self,
        data: Dict[str, torch.Tensor],
        batch: Batch,
        training: bool,
    ) -> Dict[str, torch.Tensor]:
        """Add quantities derived from the generic MACE model output."""

        return self.get_forces_stress(data, batch, training)

    def get_forces_stress(
        self, data: Dict[str, torch.Tensor], batch: Batch, training: bool
    ) -> Dict[str, torch.Tensor]:
        """
        Compute forces and stress from the energy.
        """

        for keyword in ["forces", "stress"]:  # "virials"
            if data.get(keyword) is not None:
                raise ValueError(f"'{keyword}' in 'data' should be None.")

        forces, _, stress, hessian, edge_forces = get_outputs(
            energy=data["energy"],
            positions=batch["positions"],
            cell=batch["cell"],
            displacement=data["displacement"],
            **self.instructions["forward_kwargs"],
            training=training,
        )

        to_assign = {
            "forces": forces,
            "stress": stress,
            "hessian": hessian,
            "edge_forces": edge_forces,
        }
        del data["virials"]  # virials should be computed from the stress tensor

        for keyword, value in to_assign.items():
            if keyword in data and data[keyword] is not None:
                raise ValueError(f"'{keyword}' in 'data' should be None.")
            if value is not None:
                data[keyword] = value

        return data


# -----------------------------------------------------------
# Reusable script entry point
# -----------------------------------------------------------
def run_cli(
    calculator_class=BatchedMACE,
    calculator_name="MACECalculator",
):
    """Run the standalone MACE structure-evaluation command-line interface."""

    argv = {
        "metavar": "\b",
    }

    parser = argparse.ArgumentParser(
        description=(
            f"Evaluate a MACE model on structures using {calculator_name}.\n"
            "Run with '-m mace.model -i dataset.extxyz -o output.extxyz'."
        )
    )

    parser.add_argument(
        "-m",
        "--model",
        type=str,
        required=True,
        help="Path to the trained MACE model file.",
        **argv,
    )
    parser.add_argument(
        "-d",
        "--device",
        type=str,
        required=False,
        default="cpu",
        help="torch device (default: %(default)s).",
        **argv,
    )
    parser.add_argument(
        "-mk",
        "--mace_kwargs",
        type=str,
        required=False,
        default=None,
        help="JSON file with extra input arguments for the calculator.",
        **argv,
    )
    parser.add_argument(
        "-i",
        "--input_structures",
        type=str,
        required=True,
        help="input file.",
        **argv,
    )
    parser.add_argument(
        "-o",
        "--output_structures",
        type=str,
        required=True,
        help="output file.",
        **argv,
    )
    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        required=False,
        default="MACE_",
        help="prefix for saved properties (default: %(default)s).",
        **argv,
    )
    parser.add_argument(
        "--ase_like_properties",
        "--ase-like-properties",
        dest="ase_like_properties",
        type=str,
        required=False,
        default=None,
        help=(
            "JSON file mapping additional MACE output names to ASE shapes. "
            "Use 'natoms' as the first dimension for arrays; other shapes are "
            'stored in info. Example: {"charges": ["natoms"], '
            '"polarizability": [3, 3], "free_energy": []}.'
        ),
        **argv,
    )

    args = parser.parse_args()

    print(f"Loading input structures from '{args.input_structures}'...")
    structures: List[Atoms] = read(args.input_structures, index=":")
    print(f"Loaded {len(structures)} structure(s).")

    # Load extra kwargs if provided
    mace_kwargs = {}
    if args.mace_kwargs is not None:
        print(f"Loading extra MACE kwargs from '{args.mace_kwargs}'...")
        with open(args.mace_kwargs, "r") as f:
            mace_kwargs = json.load(f)
        print("Loaded extra kwargs:", mace_kwargs)

    if args.ase_like_properties is not None:
        print(
            "Loading additional ASE-like properties from "
            f"'{args.ase_like_properties}'..."
        )
        with open(args.ase_like_properties, "r") as f:
            mace_kwargs["ase_like_properties"] = json.load(f)
        print(
            "Loaded additional ASE-like properties:",
            mace_kwargs["ase_like_properties"],
        )

    print(
        f"Initializing {calculator_name} with model '{args.model}' "
        f"on device '{args.device}'..."
    )
    calc = calculator_class(
        model_paths=args.model,
        device=args.device,
        **mace_kwargs,
    )
    print("Calculator initialized.")

    print("Evaluating structures with MACE model...")
    results: List[Parent] = calc.compute_batched(structures)
    assert len(structures) == len(results), "coding error"
    print("Evaluation complete.")

    print("Saving results into ASE Atoms objects...")
    for n, (atoms, results) in enumerate(zip(structures, results)):
        for key, value in results.items():
            shape = calc.ase_like_properties[key]
            if "natoms" in shape:
                atoms.arrays[f"{args.prefix}{key}"] = value
            else:
                atoms.info[f"{args.prefix}{key}"] = value
    print(f"Writing output structures to '{args.output_structures}'...")
    write(args.output_structures, images=structures, format="extxyz")
    print("All done!")


if __name__ == "__main__":
    run_cli()

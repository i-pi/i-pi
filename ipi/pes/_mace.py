#!/usr/bin/env python
"""An (extended) interface for the [MACE](https://github.com/ACEsuit/mace) calculator"""

import argparse
import json
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch
from ase import Atoms
from ase.io import read, write
from ase.outputs import _defineprop, all_outputs
from mace import data
from mace.calculators import MACECalculator
from mace.modules.utils import get_outputs
from mace.tools.torch_geometric.batch import Batch
from mace.tools.torch_geometric.dataloader import DataLoader

from ipi.pes._ase import ASEDriver
from ipi.pes.tools import JSONLogger, ModelResults, Parent
from ipi.utils.timing import Timer, timeit

# --------------------------------------- #
__DRIVER_NAME__ = "mace"
__DRIVER_CLASS__ = "MACE_driver"

MAX_VOLUME = 1e12

ase_like_properties = {
    "energy": (),
    "interaction_energy": (),
    "forces": ("natoms", 3),
    "displacement": (3, 3),
    "stress": (3, 3),
    "virials": (3, 3),
    "dipole": (3,),
    "atomic_dipoles": ("natoms", 3),
    "BEC": ("natoms", 9),  # ("natoms", 3, 3) is not supported by ASE
}

to_ignore_properties = ["interaction_energy", "node_feats", "node_energy"]


class MACE_driver(ASEDriver):
    """ASE driver for running MACE models with batched torch-based execution."""

    template: Atoms

    def __init__(
        self, template, model, device="cpu", mace_kwargs=None, *args, **kwargs
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
            get_extras=lambda: self.extra,
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


# --------------------------------------- #
class BatchedMACE(MACECalculator):
    """
    Batched evaluation for MACE:
     - supports batched evaluation of many atomic structures
     - can be used with a i-PI driver or as a standalone
    """

    def __init__(
        self,
        instructions: dict = None,
        get_extras: callable = None,
        *argc,
        **kwargs,
    ):
        if get_extras is not None:
            self.get_extras = get_extras

        self.instructions = instructions if instructions is not None else {}
        if "forward_kwargs" not in self.instructions:
            self.instructions["forward_kwargs"] = {}

        log = self.instructions.pop("log", None)
        self.logger = Timer(log is not None, log)
        log = self.instructions.pop("log_results", None)
        self.results_logger = JSONLogger(log)
        self.batch_size = self.instructions.pop("batch_size", 1)
        if "arrays_keys" not in kwargs:
            kwargs["arrays_keys"] = {}
        if "oxn" not in kwargs["arrays_keys"]:
            kwargs["arrays_keys"].update({"oxn": "oxn"})
        super().__init__(*argc, **kwargs)
        assert not self.use_compile, "self.use_compile=True is not supported yet."

    @timeit(name="preprocess")
    def preprocess(self, atoms: List[Atoms]) -> Tuple[DataLoader, Dict[str, bool]]:
        """
        Preprocess the calculation: prepare the batch, result tensors, etc.
        """

        keyspec = data.KeySpecification(
            info_keys=self.info_keys, arrays_keys=self.arrays_keys
        )
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

        compute_bec = False
        if "compute_BEC" in self.instructions:
            compute_bec = self.instructions["compute_BEC"]

        if compute_bec and "BEC" not in all_outputs:
            _defineprop("BEC", dtype=float, shape=("natoms", 3, 3))

        # Attention:
        # if we want to compute the Born Charges we need to call 'torch.autograd.grad' on the dipoles w.r.t. the positions.
        # However, since the forces are always computed, MACE always calls 'torch.autograd.grad' on the energy w.r.t. the positions.
        # This happens in 'compute_forces' in 'mace/modules/utils.py'.
        # If 'training' == False, in that function the computational graph will be destroy and the Born Charges can not be computed afterwards.
        # For this reason, we set 'training' == True so that the computational graph is preserved and we can call 'torch.autograd.grad' in 'compute_dielectric_gradients'.
        # If you don't believe me, please have a look at the keyword 'retain_graph' in 'mace/modules/utils.py' in the function 'compute_forces'.
        training = self.use_compile or compute_bec

        if self.model_type in ["MACE", "EnergyDipoleMACE"]:
            for n, batch in enumerate(dataset):
                # batch = next(iter(data_loader)).to(self.device)
                # batch = self._clone_batch(batch)
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
                dataset[n]["node_e0"] = node_e0
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

        data_loader = DataLoader(
            dataset,
            batch_size=self.batch_size,
            shuffle=False,
            drop_last=False,
        )

        return (
            data_loader,
            {
                "training": training,
                "compute_bec": compute_bec,
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

    @timeit(name="compute_batched", report=True)
    def compute_batched(self, atoms: List[Atoms]) -> List[Parent]:
        """
        Evaluate the model(s) on a list of structures.
        """

        data_loader, options = self.preprocess(atoms)
        training = options["training"]
        compute_bec = options["compute_bec"]
        forward_kwargs = options["forward_kwargs"]

        # model evaluation
        model_results = [
            ModelResults(ase_like_properties) for _ in range(self.num_models)
        ]
        # loop over models in the committee
        for batch_base in data_loader:
            # batch = batch_base.to(self.device)
            batch = self._clone_batch(batch_base).to_dict()
            Natoms = self.batch2natoms(batch)

            for i, model in enumerate(self.models):
                with self.logger.section("forward"):
                    out: dict[str, torch.Tensor] = model(
                        batch,
                        training=training,
                        **forward_kwargs,
                    )

                    out = self.augment_output(out, batch, training, compute_bec)

                # collect the results
                with self.logger.section("postprocess pt.1"):
                    ignored = set(to_ignore_properties)
                    if "ignore" in self.instructions:
                        ignored |= set(self.instructions["ignore"])

                    results_tensors = {}

                    for key, value in out.items():
                        if value is None or key in ignored:
                            continue

                        results_tensors[key] = value.detach().cpu().numpy()

                    model_results[i].store(Natoms, results_tensors)

        # re-order results
        with self.logger.section("postprocess pt.2"):
            out = ModelResults.mean(model_results)
            [
                self.results_logger.save(a, f"results.{n}.json")
                for n, a in enumerate(out)
            ]
        return out

    @timeit("augment_output")
    def augment_output(
        self,
        data: Dict[str, torch.Tensor],
        batch: Batch,
        training: bool,
        compute_bec: bool,
    ) -> Dict[str, torch.Tensor]:
        """Augment the output of the MACE model with derived properties such as forces, stress, and Born Effective Charges."""
        data = self.get_forces_stress(data, batch, training)

        if compute_bec and "BEC" not in data:
            bec = self.compute_dmu_dR(data, batch)
            # store to output results
            # (mu_xyz,node,R_xyz) --> (node,mu_xyz,R_xyz)
            data["BEC"] = bec.moveaxis(0, 1)

        return data

    @timeit("get_forces_stress")
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

    @timeit("compute_dmu_dR")
    def compute_dmu_dR(
        self, data: Dict[str, torch.Tensor], batch: Batch
    ) -> torch.Tensor:
        """
        Compute the derivative of the dipole (mu) w.r.t. the positions (R),
        i.e. the Born Effective Charges, or Atomic Polar Tensors.
        """

        if "dipole" not in data:
            raise ValueError(
                f"The keyword 'dipole' is not in the output data of the MACE model.\nThe data provided by the model is: {list(data.keys())}"
            )
        try:
            batch["positions"]
        except Exception:
            raise ValueError(
                f"The attribute 'positions' is not in the batch data provided to the MACE model.\nThe batch contains: {list(batch.keys())}"
            )
        dipole_components = 3
        mu = data["dipole"][:, :dipole_components]  # just for debugging
        pos = batch["positions"]
        if not isinstance(mu, torch.Tensor):
            raise ValueError(f"The dipole is not a torch.Tensor rather a {type(mu)}")
        if not isinstance(pos, torch.Tensor):
            raise ValueError(
                f"The positions are not a torch.Tensor rather a {type(pos)}"
            )

        bec = compute_dielectric_gradients(mu, [pos])[0]

        if not isinstance(bec, torch.Tensor):
            raise ValueError(
                f"The computed Born Charges are not a torch.Tensor rather a {type(bec)}"
            )
        if tuple(bec.shape) != (dipole_components, *pos.shape):
            raise ValueError(
                f"The computed Born Charges have the wrong shape. The shape {(dipole_components,*pos.shape)} was expected but got {tuple(bec.shape)}."
            )

        # Attention:
        # The tensor 'bec' has 3 dimensions.
        # Its shape is (3,*pos.shape).
        # This means that bec[0,3,2] will contain d mu_x / d R^3_z,
        # where mu_x is the x-component of the dipole and R^3_z is the z-component of the 4th (zero-indexed) atom i n the structure/batch.

        return bec


# --------------------------------------- #
# Function taken from https://github.com/davkovacs/mace/tree/mu_alpha
def compute_dielectric_gradients(
    dielectric: torch.Tensor, inputs: List[torch.Tensor], clean: Optional[bool] = False
) -> List[torch.Tensor]:
    """
    Compute gradients of the dielectric tensor with respect to a list of input tensors.

    Args:
        dielectric: Tensor whose gradients are computed.
        inputs: Tensors to differentiate with respect to (arbitrary shapes allowed).
        clean: If True, frees parts of the autograd graph when possible.

    Returns:
        List[torch.Tensor]: For each input tensor, the gradient d(dielectric)/d(input).
    """
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
            retain_graph=(i < dielectric.shape[-1] - 1)
            or not clean,  # small optimization
            create_graph=False,  # small optimization
            allow_unused=False,  # small optimization
        )
        assert len(gradients) == len(inputs), "coding error"
        for j, (gradient, input) in enumerate(zip(gradients, inputs)):
            assert gradient.shape == input.shape, "coding error"
            d_dielectric_dr[j][i] = gradient.detach()
        del gradients  # cleanup
    del grad_outputs  # cleanup
    return [torch.stack(out, dim=0) for out in d_dielectric_dr]


# -----------------------------------------------------------
# Script entry point
# -----------------------------------------------------------
if __name__ == "__main__":
    argv = {
        "metavar": "\b",
    }

    parser = argparse.ArgumentParser(
        description="Evaluate a MACE model on structures using MACECalculator.\n\
        Run with 'python extmace.py -m mace.model -i dataset.extxyz -o output.extxyz'"
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

    print(
        f"Initializing MACECalculator with model '{args.model}' on device '{args.device}'..."
    )
    calc = MACECalculator(model_paths=args.model, device=args.device, **mace_kwargs)
    print("Calculator initialized.")

    print("Evaluating structures with MACE model...")
    results: List[Parent] = calc.compute_batched(structures)
    assert len(structures) == len(results), "coding error"
    print("Evaluation complete.")

    print("Saving results into ASE Atoms objects...")
    for n, (atoms, results) in enumerate(zip(structures, results)):
        for key, value in results.items():
            shape = ase_like_properties[key]
            if "natoms" in shape:
                atoms.arrays[f"{args.prefix}{key}"] = value
            else:
                atoms.info[f"{args.prefix}{key}"] = value
    print(f"Writing output structures to '{args.output_structures}'...")
    write(args.output_structures, images=structures, format="extxyz")
    print("All done!")

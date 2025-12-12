#!/usr/bin/env python
"""An (extended) interface for the [MACE](https://github.com/ACEsuit/mace) calculator"""

import json
import torch
import numpy as np
from typing import List, Dict, Tuple, Optional

from mace import data
from mace.tools.torch_geometric.dataloader import DataLoader
from mace.tools.torch_geometric.batch import Batch
from mace.calculators import MACECalculator
from mace.modules.utils import get_outputs

from ase import Atoms
from ase.io import read
from ase.outputs import _defineprop, all_outputs

from ipi.pes.ase import ASEDriver
from ipi.pes.tools import Timer, timeit, JSONLogger, ModelResults
from ipi.utils.messages import warning, verbosity
from ipi.utils.units import unit_to_user


# --------------------------------------- #
__DRIVER_NAME__ = "extmace"
__DRIVER_CLASS__ = "Extended_MACE_driver"

DEBUG = False

ase_like_properties = {
    "energy": (),
    "node_energy": ("natoms",),
    "forces": ("natoms", 3),
    "displacement": (3, 3),
    "stress": (3, 3),
    "virials": (3, 3),
    "dipole": (3,),
    "atomic_dipoles": ("natoms", 3),
    "atomic-oxn-dipole": ("natoms", 3),
    "BEC": ("natoms", 3, 3),
    "piezoelectric": (3, 3, 3),
}


class Extended_MACE_driver(ASEDriver):
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

        self.batched_calculator = ExtendedMACECalculator(
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
            self.all_templates.append(
                [self.template.copy() for _ in range(Nstructures - Nalready)]
            )

        for n, (atoms, c, p) in enumerate(zip(self.all_templates, cell, pos)):
            atoms.set_positions(p)
            atoms.set_pbc(True)
            atoms.set_cell(c)

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
class ExtendedMACECalculator(MACECalculator):
    """
    Extended ase Calculator for MACE:
     - supports batched evaluation of many atomic structures
     - supports the inclusion of external electric fields
     - can be used with a i-PI driver or as a standalone
    """

    def __init__(
        self,
        instructions: dict = {},
        get_extras: callable = None,
        *argc,
        **kwargs,
    ):
        if get_extras is not None:
            self.get_extras = get_extras

        self.instructions = instructions
        if "forward_kwargs" not in self.instructions:
            self.instructions["forward_kwargs"] = {}
        if "ensemble" not in self.instructions:
            self.instructions["ensemble"] = "none"

        log = self.instructions.pop("log", None)
        self.logger = Timer(log is not None, log)
        log = self.instructions.pop("log_results", None)
        self.results_logger = JSONLogger(log)
        self.batch_size = self.instructions.pop("batch_size", 1)
        self.instructions = instructions
        super().__init__(*argc, **kwargs)
        assert not self.use_compile, "self.use_compile=True is not supported yet."

    def get_extras(self) -> dict:
        return {}

    @timeit(name="preprocess")
    def preprocess(self, atoms: List[Atoms]):
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
            )
            for config in configs
        ]

        compute_bec = False
        if "compute_BEC" in self.instructions:
            compute_bec = self.instructions["compute_BEC"]

        ensemble = str(self.instructions["ensemble"]).upper()
        if ensemble == "E-DEBUG":
            if not compute_bec:
                warning(
                    "'compute_bec' will be switched automatically to True since you specified 'ensemble' : 'E-debug'",
                    verbosity.high,
                )
            compute_bec = True

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
                num_atoms_arange = torch.arange(batch["positions"].shape[0])

                # this try-except is to be compatible with different MACE versions
                try:
                    # newer versions of MACE
                    node_e0 = self.models[0].atomic_energies_fn(batch["node_attrs"])[
                        num_atoms_arange, node_heads
                    ]
                except:
                    # older versions of MACE
                    node_e0 = self.models[0].atomic_energies_fn(
                        batch["node_attrs"], node_heads
                    )
                dataset[n]["node_e0"] = node_e0
            compute_stress = not self.use_compile
        else:
            compute_stress = False

        # -------------------#
        # Some extra parameters to the model
        if "compute_edge_forces" not in self.instructions["forward_kwargs"]:
            self.instructions["forward_kwargs"][
                "compute_edge_forces"
            ] = self.compute_atomic_stresses
        if "compute_stress" not in self.instructions["forward_kwargs"]:
            self.instructions["forward_kwargs"]["compute_stress"] = compute_stress

        forward_kwargs = self.instructions["forward_kwargs"].copy()

        if ensemble == "E":
            # disable all autograd flags to turn them on again in 'self.apply_ensemble'
            forward_kwargs["compute_force"] = False
            forward_kwargs["compute_stress"] = False
            forward_kwargs["compute_virials"] = False
            forward_kwargs["compute_edge_forces"] = False
            forward_kwargs["compute_displacement"] = compute_stress
            # disable Born Effective Charges computation if implemented in the model
            # if self.model_type in ["EnergyDipoleMACE"]:
            #     forward_kwargs["compute_bec"] = False
        else:
            # disable always these quantities
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
            for a in np.split(batch["positions"], batch["ptr"][1:], axis=0)[:-1]
        ]

    @timeit(name="compute_batched", report=True)
    def compute_batched(self, atoms: List[Atoms]):
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
            batch = batch_base.to(self.device)
            batch = self._clone_batch(batch_base).to_dict()
            Natoms = self.batch2natoms(batch)

            for i, model in enumerate(self.models):
                with self.logger.section("forward"):
                    out = model(
                        batch,
                        training=training,
                        **forward_kwargs,
                    )

                # apply the external electric/dielectric field
                out = self.apply_ensemble(out, batch, training, compute_bec)
                out["node_energy"] -= batch["node_e0"]

                # collect the results
                results_tensors = {}
                for key, value in out.items():
                    if (
                        "ignore" in self.instructions
                        and key in self.instructions["ignore"]
                    ):
                        continue
                    if key not in results_tensors:
                        results_tensors[key] = [None] * len(self.models)
                    results_tensors[key] = value.detach().cpu().numpy()

                model_results[i].store(Natoms, results_tensors)

        # re-order results
        with self.logger.section("postprocess"):
            out = ModelResults.mean(model_results)
            [
                self.results_logger.save(a, f"results.{n}.json")
                for n, a in enumerate(out)
            ]
        return out

    @timeit("compute_dmu_dR_deta")
    def compute_dmu_dR_deta(
        self, data: Dict[str, torch.Tensor], batch: Batch
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Compute the derivative of the dipole (mu) w.r.t. the positions (R) and lattice displacements (eta).
        The derivatives w.r.t. the positions returns the Born Effective Charges,
        while the derivatives w.r.t. the lattice displacements returns a tensor that can be related to the piezoelectric tensor.
        The conversion from this tensor to the piezoelectric one is performed in 'apply_ensemble'.
        """

        if "dipole" not in data:
            raise ValueError(
                f"The keyword 'dipole' is not in the output data of the MACE model.\nThe data provided by the model is: {list(data.keys())}"
            )
        try:
            batch["positions"]
        except:
            raise ValueError(
                f"The attribute 'positions' is not in the batch data provided to the MACE model.\nThe batch contains: {list(batch.keys())}"
            )
        dipole_components = 3
        mu = data["dipole"]  # [:,:dipole_components] # uncomment to debug
        pos = batch["positions"]
        if not isinstance(mu, torch.Tensor):
            raise ValueError(f"The dipole is not a torch.Tensor rather a {type(mu)}")
        if not isinstance(pos, torch.Tensor):
            raise ValueError(
                f"The positions are not a torch.Tensor rather a {type(pos)}"
            )

        displacement = data["displacement"]
        if displacement.requires_grad:
            res = compute_dielectric_gradients(mu, [pos, displacement])
            bec = res[0]  # (3,n_nodes,3)
            dmu_deta = res[1]  # (3,n_graphs,3,3)
        else:
            bec = compute_dielectric_gradients(mu, [pos])[0]
            dmu_deta = None

        if not isinstance(bec, torch.Tensor):
            raise ValueError(
                f"The computed Born Charges are not a torch.Tensor rather a {type(bec)}"
            )
        if tuple(bec.shape) != (dipole_components, *pos.shape):
            raise ValueError(
                f"The computed Born Charges have the wrong shape. The shape {(dipole_components,*pos.shape)} was expected but got {tuple(bec.shape)}."
            )

        if dmu_deta is not None:
            if not isinstance(dmu_deta, torch.Tensor):
                raise ValueError(
                    f"The computed piezoelectric tensor is not a torch.Tensor rather a {type(dmu_deta)}"
                )
            if tuple(dmu_deta.shape) != (dipole_components, *displacement.shape):
                raise ValueError(
                    f"The computed piezoelectric tensor has the wrong shape. The shape {(dipole_components,*displacement.shape)} was expected but got {tuple(dmu_deta.shape)}."
                )

        # Attention:
        # The tensor 'bec' has 3 dimensions.
        # Its shape is (3,*pos.shape).
        # This means that bec[0,3,2] will contain d mu_x / d R^3_z,
        # where mu_x is the x-component of the dipole and R^3_z is the z-component of the 4th (zero-indexed) atom i n the structure/batch.

        return bec, dmu_deta

    @timeit("apply_ensemble")
    def apply_ensemble(
        self,
        data: Dict[str, torch.Tensor],
        batch: Batch,
        training: bool,
        compute_bec: bool,
    ) -> Dict[str, torch.Tensor]:

        ensemble = str(self.instructions["ensemble"]).upper()
        if ensemble == "NONE":  # no ensemble (just for debugging purposes)
            pass
        else:
            extras = self.get_extras()
            if extras is None or extras == {}:
                raise ValueError("The extra information dictionary is empty.")
            Efield = np.asarray(extras["Efield"])  # in atomic units
            Efield = unit_to_user("electric-field", "v/ang", Efield)
            mu = data["dipole"]
            Efield = torch.from_numpy(Efield).to(device=self.device, dtype=mu.dtype)

            if ensemble == "E-DEBUG":  # fixed external electric field
                # This is very similar to what is done in the function 'fixed_E' in 'ipi/engine/forcefields.py'.
                if not compute_bec:
                    raise ValueError("coding error")
                bec, dmu_deta = self.compute_dmu_dR_deta(data, batch)

                interaction_energy = torch.einsum("ij,j->i", mu, Efield)
                data["energy"] -= interaction_energy
                data["forces"] += torch.einsum("ijk,i->jk", bec, Efield)

                # store to output results
                data["BEC"] = bec.moveaxis(
                    0, 1
                )  # (mu_xyz,node,R_xyz) --> (node,mu_xyz,R_xyz)

                if dmu_deta is not None:
                    cell: torch.Tensor = batch["cell"].view((-1, 3, 3))
                    volume = torch.det(cell)
                    stress_E = (
                        torch.einsum("ijkl,i->jkl", dmu_deta, Efield)
                        / volume[:, None, None]
                    )
                    data["stress"] -= stress_E

                    dmu_deta = dmu_deta.moveaxis(
                        0, 1
                    )  # (mu_xyz,graph,eta_i,eta_j) --> (graph,mu_xyz,eta_i,eta_j)
                    data["piezoelectric"] = dmu_deta2piezoelectric(
                        dmu_deta, data["dipole"], volume
                    )
                    if DEBUG:
                        # Eq. (16) of Computer Physics Communications 190 (2015) 33-50
                        cell = torch.einsum(
                            "ijk,ikl->ijl",
                            batch["cell"].view((-1, 3, 3)),
                            torch.eye(3)[None, :, :] + data["displacement"],
                        )
                        volume = torch.det(cell)
                        test = compute_dielectric_gradients(
                            data["dipole"] / volume[:, None], [data["displacement"]]
                        )[0]
                        test = test.moveaxis(
                            0, 1
                        )  # (mu_xyz,graph,eta_i,eta_j) --> (graph,mu_xyz,eta_i,eta_j)
                        assert torch.allclose(
                            test, data["piezoelectric"]
                        ), "coding error"

                        test = (
                            data["piezoelectric"]
                            + data["dipole"][:, :, None, None]
                            * torch.eye(3)[None, None, :, :]
                            / volume[:, None, None, None]
                        )
                        test = torch.einsum("ijkl,j->ikl", test, Efield)
                        assert torch.allclose(test, stress_E), "coding error"

            elif ensemble == "E":

                # Interaction energy due to the Electric Dipole Approximation
                interaction_energy = mu @ Efield
                data["energy"] -= interaction_energy

                for keyword in ["forces", "stress"]:  # "virials"
                    if data[keyword] is not None:
                        raise ValueError(f"'{keyword}' in 'data' should be None.")

                with self.logger.section("get_outputs"):
                    forces, virials, stress, hessian, edge_forces = get_outputs(
                        energy=data["energy"],
                        positions=batch["positions"],
                        cell=batch["cell"],
                        displacement=data["displacement"],
                        **self.instructions["forward_kwargs"],
                        training=training,
                    )

                to_assign = {
                    "forces": forces,
                    # "virials": virials,
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

            else:
                raise ValueError(f"Ensemble {ensemble} not implemented (yet).")

        if compute_bec and "BEC" not in data:
            bec, dmu_deta = self.compute_dmu_dR_deta(data, batch)
            # store to output results
            # (mu_xyz,node,R_xyz) --> (node,mu_xyz,R_xyz)
            data["BEC"] = bec.moveaxis(0, 1)
            if dmu_deta is not None:
                dmu_deta = dmu_deta.moveaxis(0, 1)
                cell: torch.Tensor = batch["cell"].view((-1, 3, 3))
                volume = torch.det(cell)
                data["piezoelectric"] = dmu_deta2piezoelectric(
                    dmu_deta, data["dipole"], volume
                )

        return data


# --------------------------------------- #
def dmu_deta2piezoelectric(
    dmu_deta: torch.Tensor, mu: torch.Tensor, volume: torch.Tensor
):
    """
    Convert the derivative of the dipole (mu) w.r.t. the lattice displacements (eta)
    to the piezoelectric tensor e_{ijk}.

    e_{ijk} = [ ∂μ_i / ∂η_{jk} - μ_i δ_{jk} ] / V
    """

    Nbatches = dmu_deta.shape[0]

    assert dmu_deta.shape == (
        Nbatches,
        3,
        3,
        3,
    ), f"dmu_deta must have shape (B, 3, 3, 3), got {tuple(dmu_deta.shape)}"

    assert mu.shape == (
        Nbatches,
        3,
    ), f"mu must have shape (B, 3), got {tuple(mu.shape)}"

    assert volume.shape == (
        Nbatches,
    ), f"volume must have shape (B,), got {tuple(volume.shape)}"

    delta = torch.eye(3, device=dmu_deta.device, dtype=dmu_deta.dtype)

    # Correct batched μ_i δ_jk → (B, 3, 3, 3)
    mu_delta = mu[:, :, None, None] * delta[None, None, :, :]

    e = (dmu_deta - mu_delta) / volume[:, None, None, None]

    assert e.shape == (
        Nbatches,
        3,
        3,
        3,
    ), f"output must have shape (B, 3, 3, 3), got {tuple(e.shape)}"

    return e


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
    d_dielectric_dr = d_dielectric_dr = [
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
    import argparse
    from ase.io import write

    parser = argparse.ArgumentParser(
        description="Evaluate a MACE model on structures using ExtendedMACECalculator."
    )

    parser.add_argument(
        "-m",
        "--model",
        type=str,
        required=True,
        help="Path to the trained MACE model file.",
    )
    parser.add_argument(
        "-d",
        "--device",
        type=str,
        default="cpu",
        help="Torch device (default: cpu). Example: cuda:0",
    )
    parser.add_argument(
        "-mk",
        "--mace_kwargs",
        type=str,
        default=None,
        help="JSON file with extra input arguments for the calculator.",
    )
    parser.add_argument(
        "-i",
        "--input_structures",
        type=str,
        required=True,
        help="Input file (ASE-readable).",
    )
    parser.add_argument(
        "-o",
        "--output_structures",
        type=str,
        required=True,
        help="Output file (ASE-readable).",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        default="MACE_",
        help="Prefix for saved properties.",
    )

    args = parser.parse_args()

    print(f"Loading input structures from '{args.input_structures}'...")
    structures = read(args.input_structures, index=":")
    print(f"Loaded {len(structures)} structure(s).")

    # Load extra kwargs if provided
    mace_kwargs = {}
    if args.mace_kwargs is not None:
        print(f"Loading extra MACE kwargs from '{args.mace_kwargs}'...")
        with open(args.mace_kwargs, "r") as f:
            mace_kwargs = json.load(f)
        print("Loaded extra kwargs:", mace_kwargs)

    print(
        f"Initializing ExtendedMACECalculator with model '{args.model}' on device '{args.device}'..."
    )
    calc = ExtendedMACECalculator(
        model_path=args.model, device=args.device, **mace_kwargs
    )
    print("Calculator initialized.")

    print("Evaluating structures with MACE model...")
    results = calc.compute_batched(structures)
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

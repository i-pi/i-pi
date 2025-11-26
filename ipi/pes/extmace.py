"""An (extended) interface for the [MACE](https://github.com/ACEsuit/mace) calculator"""

from .mace import MACE_driver
from .ase import ASEDriver


# --------------------------------------- #
__DRIVER_NAME__ = "extmace"
__DRIVER_CLASS__ = "Extended_MACE_driver"


class Extended_MACE_driver(MACE_driver):
    """
    MACE driver with the torch tensors exposed.
    """

    def check_parameters(self):
        """Check the arguments requuired to run the driver

        This loads the potential and atoms template in MACE
        """

        ASEDriver.check_parameters(self)  # Explicitly bypass MACE_driver

        self.ase_calculator = Extended_MACECalculator(
            model_paths=self.model,
            device=self.device,
            get_extras=lambda: self.extra,
            **self.mace_kwargs,
        )


# --------------------------------------- #
import torch
import numpy as np
from typing import List, Dict, Tuple
from mace.tools.torch_geometric.batch import Batch
from mace.calculators import MACECalculator
from mace.modules.utils import get_outputs
from ase.calculators.calculator import Calculator, all_changes
from ase.outputs import _defineprop, all_outputs
from .tools import Timer, timeit
from ipi.utils.messages import warning, verbosity
from ipi.utils.units import unit_to_user


class Extended_MACECalculator(MACECalculator):
    def __init__(
        self, instructions: dict = {}, get_extras: callable = None, *argc, **kwargs
    ):
        if get_extras is not None:
            self.get_extras = get_extras

        self.instructions = instructions
        if "forward_kwargs" not in self.instructions:
            self.instructions["forward_kwargs"] = {}
        if "ensemble" not in self.instructions:
            self.instructions["ensemble"] = "none"

        log = self.instructions.get("log", None)
        self.logger = Timer(log is not None, log)

        super().__init__(*argc, **kwargs)

    def get_extras(self) -> dict:
        return {}

    @timeit(name="preprocess()")
    def preprocess(self, atoms=None):
        """
        Preprocess the calculation: prepare the batch, result tensors, etc.
        """

        Calculator.calculate(self, atoms)

        assert not self.use_compile, "self.use_compile=True is not supported yet."

        with self.logger.section("Prepare batch"):
            batch_base: Batch = self._atoms_to_batch(atoms)

        compute_bec = False
        if "compute_BEC" in self.instructions:
            compute_bec = self.instructions["compute_BEC"]

        ensemble = str(self.instructions["ensemble"]).upper()
        if ensemble == "E-DEBUG":
            if not compute_bec:
                warning(
                    f"'compute_bec' will be switched automatically to True since you specified 'ensemble' : 'E-debug'",
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
            batch = self._clone_batch(batch_base)
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
            compute_stress = not self.use_compile
        else:
            compute_stress = False

        # -------------------#
        # Some extra parameters to the model
        if "compute_edge_forces" not in self.instructions["forward_kwargs"]:
            self.instructions["forward_kwargs"][
                "compute_edge_forces"
            ] = self.compute_atomic_stresses
        # if "compute_atomic_stresses" not in self.instructions["forward_kwargs"]:
        #     self.instructions["forward_kwargs"]["compute_atomic_stresses"] = (
        #         self.compute_atomic_stresses
        #     )
        if "compute_stress" not in self.instructions["forward_kwargs"]:
            self.instructions["forward_kwargs"]["compute_stress"] = compute_stress
        # if "compute_hessian" not in self.instructions["forward_kwargs"]:
        #     self.instructions["forward_kwargs"]["compute_hessian"] = False

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

        ret_tensors = self._create_result_tensors(
            self.model_type, self.num_models, len(atoms)
        )
        if "energy" not in ret_tensors:
            ret_tensors["energy"] = ret_tensors["energies"].clone()
        if compute_bec:
            f = ret_tensors["forces"]
            ret_tensors["BEC"] = torch.zeros(
                (len(self.models), *f.shape[1:], 3), device=self.device
            )

        return (
            batch_base,
            ret_tensors,
            node_e0,
            {
                "training": training,
                "compute_bec": compute_bec,
                "forward_kwargs": forward_kwargs,
            },
        )

    @timeit(name="calculate()", report=True)
    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        """
        Calculate properties.
        :param atoms: ase.Atoms object
        :param properties: [str], properties to be computed, used by ASE internally
        :param system_changes: [str], system changes since last calculation, used by ASE internally
        :return:
        """

        batch_base, ret_tensors, node_e0, options = self.preprocess(atoms)
        training = options["training"]
        compute_bec = options["compute_bec"]
        forward_kwargs = options["forward_kwargs"]

        # model evaluation
        with self.logger.section("Forward pass loop"):
            # loop over models in the committee
            for i, model in enumerate(self.models):
                with self.logger.section(f"model {i} -> '_clone_batch' "):
                    batch = self._clone_batch(batch_base)

                with self.logger.section(f"model {i} -> 'forward'"):
                    out = model(
                        batch.to_dict(),
                        training=training,
                        **forward_kwargs,
                    )

                # apply the external electric/dielectric field
                out = self.apply_ensemble(out, batch, training, compute_bec)

                # collect the output
                for keyword in ["energy", "forces", "stress", "dipole", "BEC"]:
                    if keyword in out:
                        ret_tensors[keyword][i] = out[keyword].detach()

                # for keyword in ["atomic_stresses", "atomic_virials"]:
                #     if keyword in out:
                #         ret_tensors.setdefault(keyword, []).append(
                #             out[keyword].detach()
                #         )

                # if "node_energy" in out:
                #     ret_tensors["node_energy"][i] = (
                #         out["node_energy"] - node_e0
                #     ).detach()

        # remove properties
        if "ignore" in self.instructions:
            for k in self.instructions["ignore"]:  # List[str]
                del ret_tensors[k]

        self.postprocess(ret_tensors)

    @timeit(name="'postprocess'")
    def postprocess(self, ret_tensors: Dict[str, torch.Tensor]):

        # process outputs
        self.results = {}
        # if self.model_type in ["MACE", "EnergyDipoleMACE"]:
        self.results["energy"] = (
            torch.mean(ret_tensors["energy"], dim=0).cpu().item()
            * self.energy_units_to_eV
        )
        self.results["free_energy"] = self.results["energy"]
        if "node_energy" in ret_tensors:
            self.results["node_energy"] = (
                torch.mean(ret_tensors["node_energy"], dim=0).cpu().numpy()
            ) * self.energy_units_to_eV
        self.results["forces"] = (
            torch.mean(ret_tensors["forces"], dim=0).cpu().numpy()
            * self.energy_units_to_eV
            / self.length_units_to_A
        )
        if self.num_models > 1:
            energies = ret_tensors["energy"].cpu().numpy() * self.energy_units_to_eV
            self.results["energy_var"] = (
                torch.var(energies, dim=0, unbiased=False).cpu().item()
                * self.energy_units_to_eV
            )
            self.results["forces_comm"] = (
                ret_tensors["forces"].cpu().numpy()
                * self.energy_units_to_eV
                / self.length_units_to_A
            )
        if ret_tensors["stress"] is not None:
            try:
                from ase.stress import full_3x3_to_voigt_6_stress
            except:
                raise ImportError("Couldn't load full_3x3_to_voigt_6_stress from ase.")
            self.results["stress"] = full_3x3_to_voigt_6_stress(
                torch.mean(ret_tensors["stress"], dim=0).cpu().numpy()
                * self.energy_units_to_eV
                / self.length_units_to_A**3
            )
            if self.num_models > 1:
                self.results["stress_var"] = full_3x3_to_voigt_6_stress(
                    torch.var(ret_tensors["stress"], dim=0, unbiased=False)
                    .cpu()
                    .numpy()
                    * self.energy_units_to_eV
                    / self.length_units_to_A**3
                )
            # if "atomic_stresses" in ret_tensors:
            #     self.results["stresses"] = (
            #         torch.mean(torch.stack(ret_tensors["atomic_stresses"]), dim=0)
            #         .cpu()
            #         .numpy()
            #         * self.energy_units_to_eV
            #         / self.length_units_to_A**3
            #     )
            # if "atomic_virials" in ret_tensors:
            #     self.results["virials"] = (
            #         torch.mean(torch.stack(ret_tensors["atomic_virials"]), dim=0)
            #         .cpu()
            #         .numpy()
            #         * self.energy_units_to_eV
            #     )
        # if self.model_type in ["DipoleMACE", "EnergyDipoleMACE"]:
        if "dipole" in ret_tensors:
            self.results["dipole"] = (
                torch.mean(ret_tensors["dipole"], dim=0).cpu().numpy()
            )
        if "BEC" in ret_tensors:
            self.results["BEC"] = torch.mean(ret_tensors["BEC"], dim=0).cpu().numpy()

    @timeit("'compute_BEC_piezo'")
    def compute_BEC_piezo(
        self, data: Dict[str, torch.Tensor], batch: Batch
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        if "dipole" not in data:
            raise ValueError(
                f"The keyword 'dipole' is not in the output data of the MACE model.\nThe data provided by the model is: {list(data.keys())}"
            )
        try:
            batch.positions
        except:
            raise ValueError(
                f"The attribute 'positions' is not in the batch data provided to the MACE model.\nThe batch contains: {list(batch.__dict__.keys())}"
            )
        dipole_components = 3
        mu = data["dipole"]  # [:,:dipole_components] # uncomment to debug
        pos = batch.positions
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
            piezo = res[1]  # (3,n_graphs,3,3)
        else:
            bec = compute_dielectric_gradients(mu, [pos])[0]
            piezo = None

        if not isinstance(bec, torch.Tensor):
            raise ValueError(
                f"The computed Born Charges are not a torch.Tensor rather a {type(bec)}"
            )
        if tuple(bec.shape) != (dipole_components, *pos.shape):
            raise ValueError(
                f"The computed Born Charges have the wrong shape. The shape {(dipole_components,*pos.shape)} was expected but got {tuple(bec.shape)}."
            )

        if piezo is not None:
            if not isinstance(piezo, torch.Tensor):
                raise ValueError(
                    f"The computed piezoelectric tensor is not a torch.Tensor rather a {type(piezo)}"
                )
            if tuple(piezo.shape) != (dipole_components, *displacement.shape):
                raise ValueError(
                    f"The computed piezoelectric tensor has the wrong shape. The shape {(dipole_components,*displacement.shape)} was expected but got {tuple(piezo.shape)}."
                )

        # Attention:
        # The tensor 'bec' has 3 dimensions.
        # Its shape is (3,*pos.shape).
        # This means that bec[0,3,2] will contain d mu_x / d R^3_z,
        # where mu_x is the x-component of the dipole and R^3_z is the z-component of the 4th (zero-indexed) atom i n the structure/batch.

        return bec, piezo

    @timeit("'apply_ensemble'")
    def apply_ensemble(
        self,
        data: Dict[str, torch.Tensor],
        batch: Batch,
        training: bool,
        compute_bec: bool,
    ) -> Dict[str, torch.Tensor]:
        # if self.instructions == {}:
        #     return data

        ensemble = str(self.instructions["ensemble"]).upper()
        if ensemble == "NONE":  # no ensemble (just for debugging purposes)
            # compute Born Effective Charges using autodiff
            pass
        else:
            extras = self.get_extras()
            if extras is None or extras == {}:
                raise ValueError("The extra information dictionary is empty.")
            Efield = np.asarray(extras["Efield"])  # in atomic units
            Efield = unit_to_user("electric-field", "v/ang", Efield)
            mu = data["dipole"][0]
            Efield = torch.from_numpy(Efield).to(device=self.device, dtype=mu.dtype)

            if ensemble == "E-DEBUG":  # fixed external electric field
                # This is very similar to what is done in the function 'fixed_E' in 'ipi/engine/forcefields.py'.
                if not compute_bec:
                    raise ValueError("coding error")
                bec, piezo = self.compute_BEC_piezo(data, batch)
                data["BEC"] = bec.moveaxis(0, 2).moveaxis(1, 2)
                if piezo is not None:
                    data["piezo"] = piezo

                interaction_energy = mu @ Efield
                data["energy"] -= interaction_energy
                data["forces"] += torch.einsum("ijk,i->jk", bec, Efield)
                if piezo is not None:
                    volume = torch.det(batch.cell)
                    data["stress"] -= (
                        torch.einsum("ijkl,i->jkl", piezo, Efield)
                        + interaction_energy / volume
                    )

            elif ensemble == "E":
                interaction_energy = mu @ Efield
                data["energy"] -= interaction_energy

                for keyword in ["forces", "stress"]:  # "virials"
                    if data[keyword] is not None:
                        raise ValueError(f"'{keyword}' in 'data' should be None.")

                with self.logger.section("'get_outputs'"):
                    forces, virials, stress, hessian, edge_forces = get_outputs(
                        energy=data["energy"],
                        positions=batch.positions,
                        cell=batch.cell,
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

                for keyword, value in to_assign.items():
                    if keyword in data and data[keyword] is not None:
                        raise ValueError(f"'{keyword}' in 'data' should be None.")
                    data[keyword] = value

            elif ensemble == "D":  # fixed dielectric displacement
                ValueError("Not implemented yet")

            else:
                raise ValueError("coding error")

        if compute_bec and "BEC" not in data:
            bec, piezo = self.compute_BEC_piezo(data, batch)
            data["BEC"] = bec.moveaxis(0, 2).moveaxis(1, 2)
            if piezo is not None:
                data["piezo"] = piezo

        return data


# --------------------------------------- #
# Function taken from https://github.com/davkovacs/mace/tree/mu_alpha
def compute_dielectric_gradients(
    dielectric: torch.Tensor, inputs: List[torch.Tensor]
) -> torch.Tensor:
    """Compute the spatial derivatives of dielectric tensor.

    Args:
        dielectric (torch.Tensor): Dielectric tensor.
        inputs (List of torch.Tensor): Atom positions.

    Returns:
        torch.Tensor: Spatial derivatives of dielectric tensor.
    """
    # dielectric = dielectric[:,0:2]
    d_dielectric_dr = d_dielectric_dr = [
        [None for _ in range(dielectric.shape[-1])] for _ in range(len(inputs))
    ]
    # [Pxyz,atoms,Rxyz]
    grad_outputs: List[torch.Tensor] = [
        torch.ones((dielectric.shape[0], 1)).to(dielectric.device)
    ]
    for i in range(dielectric.shape[-1]):
        gradients = torch.autograd.grad(
            outputs=[dielectric[:, i].unsqueeze(-1)],
            inputs=inputs,
            grad_outputs=grad_outputs,
            retain_graph=(i < dielectric.shape[-1] - 1),  # small optimization
            create_graph=False,  # small optimization
            allow_unused=False,  # small optimization
        )
        assert len(gradients) == len(inputs), "coding error"
        for j, (gradient, input) in enumerate(zip(gradients, inputs)):
            assert gradient.shape == input.shape, "coding error"
            d_dielectric_dr[j][i] = gradient.detach()
        del gradients  # cleanup
    del grad_outputs  # cleanup
    return [torch.stack(out, dim=0) for out in d_dielectric_dr]  # [Pxyz,atoms,Rxyz]

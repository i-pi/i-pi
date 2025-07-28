"""An interface for the [MACE](https://github.com/ACEsuit/mace) calculator"""

import json

from .mace import MACE_driver

# MACECalculator = None
# all_changes = None
# Calculator = None
# torch = None

__DRIVER_NAME__ = "extmace"
__DRIVER_CLASS__ = "Extended_MACE_driver"


class Extended_MACE_driver(MACE_driver):
    """
    MACE driver with the torch tensors exposed.
    """

    # def __init__(
    #     self, instructions: str = None, *args, **kwargs
    # ):
    #     # self.forward_kwargs = {}
    #     # if forward_kwargs is not None:
    #     #     with open(forward_kwargs, "r") as f:
    #     #         self.forward_kwargs = json.load(f)

    #     # self.instructions = {}
    #     # if instructions is not None:
    #     #     with open(instructions, "r") as f:
    #     #         self.instructions = json.load(f)
    #     # if "forward_kwargs" not in self.instructions:
    #     #     self.instructions["forward_kwargs"] = {}

    #     super().__init__(*args, **kwargs)

    def check_parameters(self):
        """Check the arguments requuired to run the driver

        This loads the potential and atoms template in MACE
        """

        super().check_parameters()

        self.ase_calculator = Extended_MACECalculator(
            model_paths=self.model,
            device=self.device,
            # instructions=self.instructions,
            # forward_kwargs=self.forward_kwargs,
            get_extras=lambda: self.extra,
            **self.mace_kwargs
        )

OK = True
try:
    import torch
    from typing import List, Optional
    from mace.tools.torch_geometric.batch import Batch
    from mace.calculators import MACECalculator
    # from mace.tools.torch_geometric.dataloader import DataLoader
    from ase.calculators.calculator import Calculator, all_changes
except:
    OK = False
    
if OK:

    class Extended_MACECalculator(MACECalculator):
        
        def __init__(
            self,
            instructions: dict = {},
            get_extras: callable = None,
            *argc,
            **kwargs
        ):
            if get_extras is not None:
                self.get_extras = get_extras
            
            self.instructions = instructions            
            if "forward_kwargs" not in self.instructions:
                self.instructions["forward_kwargs"] = {}
                
            super().__init__(*argc, **kwargs)
            
        def get_extras(self) -> dict:
            return {}

        def apply_ensemble(self, data: dict) -> dict:
            if self.instructions == {}:
                return data
            if "ensemble" not in self.instructions:
                raise ValueError(
                    "You need to specify an ensemble (among 'none', 'E', and 'D') in the instructions file."
                )
            ensemble = str(self.instructions["ensemble"]).upper()
            if ensemble == "NONE":  # no ensemble (just for debugging purposes)
                pass
            elif ensemble == "E":  # fixed external electric field
                pass
            elif ensemble == "D":  # fixed dielectric displacement
                pass
            return data

        def calculate(self, atoms=None, properties=None, system_changes=all_changes):
            """
            Calculate properties.
            :param atoms: ase.Atoms object
            :param properties: [str], properties to be computed, used by ASE internally
            :param system_changes: [str], system changes since last calculation, used by ASE internally
            :return:
            """
            # call to base-class to set atoms attribute
            Calculator.calculate(self, atoms)

            batch_base:Batch = self._atoms_to_batch(atoms)
            assert isinstance(batch_base, Batch), "hello"
            
            #-------------------#
            # Some boolean flags
            compute_bec = "compute_BEC" in self.instructions and self.instructions["compute_BEC"] == True
            
            # Attention:
            # if we want to compute the Born Charges we need to call 'torch.autograd.grad' on the dipoles w.r.t. the positions.
            # However, since the forces are always computed, MACE always calls 'torch.autograd.grad' on the energy w.r.t. the positions.
            # This happens in 'compute_forces' in '/mace/modules/utils.py'.
            # If 'training' == False, in that function the computational graph will be destroy and the Born Charges can not be computed afterwards.
            # For this reason, we set 'training' == True so that the computational graph is preserved and we can call 'torch.autograd.grad' in 'compute_dielectric_gradients'.
            # If you don't believe me, please have a look at the keyword 'retain_graph' in '/mace/modules/utils.py' in the function 'compute_forces'.
            training = self.use_compile or compute_bec

            if self.model_type in ["MACE", "EnergyDipoleMACE"]:
                batch = self._clone_batch(batch_base)
                node_heads = batch["head"][batch["batch"]]
                num_atoms_arange = torch.arange(batch["positions"].shape[0])
                node_e0 = self.models[0].atomic_energies_fn(batch["node_attrs"])[
                    num_atoms_arange, node_heads
                ]
                compute_stress = not self.use_compile
            else:
                compute_stress = False

            ret_tensors = self._create_result_tensors(
                self.model_type, self.num_models, len(atoms)
            )
            for i, model in enumerate(self.models):
                batch = self._clone_batch(batch_base)
                out = model(
                    batch.to_dict(),
                    compute_stress=compute_stress,
                    training=training, # yes, this is correct
                    compute_edge_forces=self.compute_atomic_stresses,
                    compute_atomic_stresses=self.compute_atomic_stresses,
                    **self.instructions["forward_kwargs"]
                )

                # compute Born Effective Charges using autodiff
                if compute_bec:
                    out["BEC"] = self.compute_BEC(out,batch)

                # apply the external electric/dielectric field
                out = self.apply_ensemble(out)

                if self.model_type in ["MACE", "EnergyDipoleMACE"]:
                    ret_tensors["energies"][i] = out["energy"].detach()
                    ret_tensors["node_energy"][i] = (
                        out["node_energy"] - node_e0
                    ).detach()
                    ret_tensors["forces"][i] = out["forces"].detach()
                    if out["stress"] is not None:
                        ret_tensors["stress"][i] = out["stress"].detach()
                if self.model_type in ["DipoleMACE", "EnergyDipoleMACE"]:
                    ret_tensors["dipole"][i] = out["dipole"].detach()
                if self.model_type in ["MACE"]:
                    if out["atomic_stresses"] is not None:
                        ret_tensors.setdefault("atomic_stresses", []).append(
                            out["atomic_stresses"].detach()
                        )
                    if out["atomic_virials"] is not None:
                        ret_tensors.setdefault("atomic_virials", []).append(
                            out["atomic_virials"].detach()
                        )

            del out

            # remove properties
            if "ignore" in self.instructions:
                for k in self.instructions["ignore"]:  # List[str]
                    del ret_tensors[k]

            self.results = {}
            if self.model_type in ["MACE", "EnergyDipoleMACE"]:
                self.results["energy"] = (
                    torch.mean(ret_tensors["energies"], dim=0).cpu().item()
                    * self.energy_units_to_eV
                )
                self.results["free_energy"] = self.results["energy"]
                self.results["node_energy"] = (
                    torch.mean(ret_tensors["node_energy"], dim=0).cpu().numpy()
                )
                self.results["forces"] = (
                    torch.mean(ret_tensors["forces"], dim=0).cpu().numpy()
                    * self.energy_units_to_eV
                    / self.length_units_to_A
                )
                if self.num_models > 1:
                    self.results["energies"] = (
                        ret_tensors["energies"].cpu().numpy() * self.energy_units_to_eV
                    )
                    self.results["energy_var"] = (
                        torch.var(ret_tensors["energies"], dim=0, unbiased=False)
                        .cpu()
                        .item()
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
                        raise ImportError(
                            "Couldn't load full_3x3_to_voigt_6_stress from ase."
                        )
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
                if "atomic_stresses" in ret_tensors:
                    self.results["stresses"] = (
                        torch.mean(torch.stack(ret_tensors["atomic_stresses"]), dim=0)
                        .cpu()
                        .numpy()
                        * self.energy_units_to_eV
                        / self.length_units_to_A**3
                    )
                if "atomic_virials" in ret_tensors:
                    self.results["virials"] = (
                        torch.mean(torch.stack(ret_tensors["atomic_virials"]), dim=0)
                        .cpu()
                        .numpy()
                        * self.energy_units_to_eV
                    )
            if self.model_type in ["DipoleMACE", "EnergyDipoleMACE"]:
                self.results["dipole"] = (
                    torch.mean(ret_tensors["dipole"], dim=0).cpu().numpy()
                )
                if self.num_models > 1:
                    self.results["dipole_var"] = (
                        torch.var(ret_tensors["dipole"], dim=0, unbiased=False)
                        .cpu()
                        .numpy()
                    )

        @staticmethod
        def compute_BEC(data: dict,batch:Batch):
            # still to be implemented
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
            mu = data["dipole"]
            pos = batch.positions
            if not isinstance(mu,torch.Tensor):
                raise ValueError(f"The dipole is not a torch.Tensor rather a {type(mu)}")
            if not isinstance(pos,torch.Tensor):
                raise ValueError(f"The positions are not a torch.Tensor rather a {type(pos)}")
            bec = compute_dielectric_gradients(mu,pos)
            if not isinstance(bec,torch.Tensor):
                raise ValueError(f"The computed Born Charges are not a torch.Tensor rather a {type(bec)}")
            if tuple(bec.shape) != (3,*pos.shape):
                raise ValueError(f"The computed Born Charges have the wrong shape. The shape {(*pos.shape,3)} was expected but got {tuple(bec.shape)}.")
            data["BEC"] = bec
            # Attention:
            # The tensor 'bec' has 3 dimensions.
            # Its shape is (3,*pos.shape).
            # This means that bec[0,3,2] will contain d mu_x / d R^3_z,
            # where mu_x is the x-component of the dipole and R^3_z is the z-component of the 4th (zero-indexed) atom i n the structure/batch.
            return data

#---------------------------------------#
# Function taken from https://github.com/davkovacs/mace/tree/mu_alpha
def compute_dielectric_gradients(dielectric: torch.Tensor, positions: torch.Tensor) -> torch.Tensor:
    """Compute the spatial derivatives of dielectric tensor.

    Args:
        dielectric (torch.Tensor): Dielectric tensor.
        positions (torch.Tensor): Atom positions.

    Returns:
        torch.Tensor: Spatial derivatives of dielectric tensor.
    """
    # dielectric = dielectric[:,0:2]
    d_dielectric_dr = [None]*dielectric.shape[-1]
    for i in range(dielectric.shape[-1]):
        grad_outputs: List[Optional[torch.Tensor]] = [
            torch.ones((dielectric.shape[0], 1)).to(dielectric.device)
        ]
        gradient = torch.autograd.grad(
            outputs=[dielectric[:, i].unsqueeze(-1)],
            inputs=[positions],
            grad_outputs=grad_outputs,
            retain_graph=True,
            create_graph=True,
            allow_unused=True,
        )[0]
        d_dielectric_dr[i] = gradient
    # output = torch.stack(d_dielectric_dr, dim=0)
    # if gradient is None:
    #     return torch.zeros((positions.shape[0], dielectric.shape[-1], 3))
    return torch.stack(d_dielectric_dr, dim=0) # [Pxyz,atoms,Rxyz]
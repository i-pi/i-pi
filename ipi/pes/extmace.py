"""An (extended) interface for the [MACE](https://github.com/ACEsuit/mace) calculator"""

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

    def check_parameters(self):
        """Check the arguments requuired to run the driver

        This loads the potential and atoms template in MACE
        """

        super().check_parameters()

        self.ase_calculator = Extended_MACECalculator(
            model_paths=self.model,
            device=self.device,
            get_extras=lambda: self.extra,
            **self.mace_kwargs,
        )


OK = True
try:
    import time
    from datetime import datetime
    import torch
    import numpy as np
    from typing import List, Optional, Dict
    from mace.tools.torch_geometric.batch import Batch
    from mace.calculators import MACECalculator
    from mace.modules.utils import get_outputs
    from ase.calculators.calculator import Calculator, all_changes

    LOG_ENABLED = False
except:
    OK = False

if OK:

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
            if "log" not in self.instructions:
                self.instructions["log"] = False

            global LOG_ENABLED
            LOG_ENABLED = self.instructions["log"]

            super().__init__(*argc, **kwargs)

        def get_extras(self) -> dict:
            return {}

        def calculate(self, atoms=None, properties=None, system_changes=all_changes):
            """
            Calculate properties.
            :param atoms: ase.Atoms object
            :param properties: [str], properties to be computed, used by ASE internally
            :param system_changes: [str], system changes since last calculation, used by ASE internally
            :return:
            """

            # time single operations
            with Timer("calculate()", LOG_ENABLED):
                # call to base-class to set atoms attribute
                with Timer("Base calculate", LOG_ENABLED):
                    Calculator.calculate(self, atoms)

                with Timer("Prepare batch", LOG_ENABLED):
                    batch_base: Batch = self._atoms_to_batch(atoms)

                compute_bec = False
                if "compute_BEC" in self.instructions:
                    compute_bec = self.instructions["compute_BEC"]

                if self.instructions["ensemble"] == "E-debug":
                    compute_bec = True

                if compute_bec:
                    from ase.outputs import _defineprop, all_outputs

                    if "BEC" not in all_outputs:
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
                    node_e0 = self.models[0].atomic_energies_fn(batch["node_attrs"])[
                        num_atoms_arange, node_heads
                    ]
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
                    self.instructions["forward_kwargs"][
                        "compute_stress"
                    ] = compute_stress
                # if "compute_hessian" not in self.instructions["forward_kwargs"]:
                #     self.instructions["forward_kwargs"]["compute_hessian"] = False

                forward_kwargs = self.instructions["forward_kwargs"].copy()

                if self.instructions["ensemble"] == "E":
                    # disable all autograd flags to turn them on again in 'self.apply_ensemble'
                    forward_kwargs["compute_force"] = False
                    forward_kwargs["compute_stress"] = False
                    forward_kwargs["compute_virials"] = False
                    # forward_kwargs["compute_hessian"] = False
                    forward_kwargs["compute_edge_forces"] = False
                    # forward_kwargs["compute_atomic_stresses"] = False

                ret_tensors = self._create_result_tensors(
                    self.model_type, self.num_models, len(atoms)
                )
                if compute_bec:
                    f = ret_tensors["forces"]
                    ret_tensors["BEC"] = torch.zeros(
                        (len(self.models), *f.shape[1:], 3), device=self.device
                    )
                    # torch.zeros((3,*f.shape),dtype=f.dtype,device=f.device)

                with Timer("Forward pass loop", LOG_ENABLED):
                    for i, model in enumerate(self.models):
                        with Timer(f"model {i} -> '_clone_batch' ", LOG_ENABLED):
                            batch = self._clone_batch(batch_base)

                        with Timer(f"model {i} -> 'forward'", LOG_ENABLED):
                            out = model(
                                batch.to_dict(),
                                training=training,
                                **forward_kwargs,
                                compute_displacement=True,
                            )

                        # apply the external electric/dielectric field
                        with Timer(f"model {i} -> 'apply_ensemble'", LOG_ENABLED):
                            out = self.apply_ensemble(out, batch, training, compute_bec)

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
                        if "BEC" in out:
                            ret_tensors["BEC"][i] = out["BEC"].detach()

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
                            ret_tensors["energies"].cpu().numpy()
                            * self.energy_units_to_eV
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
                            torch.mean(
                                torch.stack(ret_tensors["atomic_stresses"]), dim=0
                            )
                            .cpu()
                            .numpy()
                            * self.energy_units_to_eV
                            / self.length_units_to_A**3
                        )
                    if "atomic_virials" in ret_tensors:
                        self.results["virials"] = (
                            torch.mean(
                                torch.stack(ret_tensors["atomic_virials"]), dim=0
                            )
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
                if "BEC" in ret_tensors:
                    self.results["BEC"] = (
                        torch.mean(ret_tensors["BEC"], dim=0).cpu().numpy()
                    )

            if LOG_ENABLED:
                Timer.report()
            return

        @staticmethod
        def compute_BEC(data: Dict[str, torch.Tensor], batch: Batch) -> torch.Tensor:
            with Timer("'compute_BEC'", LOG_ENABLED):
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
                if not isinstance(mu, torch.Tensor):
                    raise ValueError(
                        f"The dipole is not a torch.Tensor rather a {type(mu)}"
                    )
                if not isinstance(pos, torch.Tensor):
                    raise ValueError(
                        f"The positions are not a torch.Tensor rather a {type(pos)}"
                    )
                bec = compute_dielectric_gradients(mu, pos)
                if not isinstance(bec, torch.Tensor):
                    raise ValueError(
                        f"The computed Born Charges are not a torch.Tensor rather a {type(bec)}"
                    )
                if tuple(bec.shape) != (3, *pos.shape):
                    raise ValueError(
                        f"The computed Born Charges have the wrong shape. The shape {(*pos.shape,3)} was expected but got {tuple(bec.shape)}."
                    )
                # Attention:
                # The tensor 'bec' has 3 dimensions.
                # Its shape is (3,*pos.shape).
                # This means that bec[0,3,2] will contain d mu_x / d R^3_z,
                # where mu_x is the x-component of the dipole and R^3_z is the z-component of the 4th (zero-indexed) atom i n the structure/batch.

                # we need to reshape 'bec' so that the first two axis will become d.o.f. once flattened
                # and the last axis will the the component of the dipole
                bec = bec.moveaxis(0, 2).moveaxis(1, 2)
                # bec.sum(dim=0) should be a (3,3) matrix filled with zeros.
                return bec

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
                Efield = np.asarray(extras["Efield"])
                mu = data["dipole"][0]
                Efield = torch.from_numpy(Efield).to(device=self.device, dtype=mu.dtype)

                if ensemble == "E-debug":  # fixed external electric field
                    # This is very similar to what is done in the function 'fixed_E' in 'ipi/engine/forcefields.py'.
                    if not compute_bec:
                        raise ValueError("coding error")
                    data["BEC"] = self.compute_BEC(data, batch)

                    Z = data["BEC"]
                    data["energy"] -= mu @ Efield
                    data["forces"] += torch.einsum("ijk,i->jk", Z, Efield)

                elif ensemble == "E":
                    data["energy"] -= mu @ Efield

                    if data["forces"] is not None:
                        raise ValueError("'forces' in 'data' should be None.")

                    with Timer("'get_outputs'", LOG_ENABLED):
                        forces, virials, stress, hessian, edge_forces = get_outputs(
                            energy=data["energy"],
                            positions=batch.positions,
                            cell=batch.cell,
                            displacement=data["displacement"],
                            **self.instructions["forward_kwargs"],
                            training=training,
                        )

                    if forces is not None:
                        data["forces"] = forces
                    if virials is not None:
                        data["virials"] = virials
                    if stress is not None:
                        data["stress"] = stress
                    if hessian is not None:
                        data["hessian"] = hessian
                    if edge_forces is not None:
                        data["edge_forces"] = edge_forces

                elif ensemble == "D":  # fixed dielectric displacement
                    ValueError("Not implemented yet")

                else:
                    raise ValueError("coding error")

            if compute_bec and "BEC" not in data:
                data["BEC"] = self.compute_BEC(data, batch)

            return data

    # ---------------------------------------#
    # Function taken from https://github.com/davkovacs/mace/tree/mu_alpha
    def compute_dielectric_gradients(
        dielectric: torch.Tensor, positions: torch.Tensor
    ) -> torch.Tensor:
        """Compute the spatial derivatives of dielectric tensor.

        Args:
            dielectric (torch.Tensor): Dielectric tensor.
            positions (torch.Tensor): Atom positions.

        Returns:
            torch.Tensor: Spatial derivatives of dielectric tensor.
        """
        # dielectric = dielectric[:,0:2]
        d_dielectric_dr = [None] * dielectric.shape[-1]
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
        return torch.stack(d_dielectric_dr, dim=0)  # [Pxyz,atoms,Rxyz]

    class Timer:
        _stack = []
        _records = []

        def __init__(self, name, log=True):
            self.name = name
            self.log = log
            self.level = len(Timer._stack)

        def __enter__(self):
            self.start = time.time()
            Timer._stack.append(self)
            return self

        def __exit__(self, exc_type, exc_val, exc_tb):
            self.end = time.time()
            self.elapsed = self.end - self.start
            Timer._stack.pop()
            if self.log:
                Timer._records.append((self.start, self.level, self.name, self.elapsed))

        @staticmethod
        def report():
            for t, level, name, elapsed in sorted(Timer._records):
                timestamp = datetime.fromtimestamp(t).strftime("%Y-%m-%d %H:%M:%S")
                indent = "    " * level
                print(f"{timestamp}{indent} {name}: {elapsed:.4f} s")

        @staticmethod
        def reset():
            Timer._records.clear()

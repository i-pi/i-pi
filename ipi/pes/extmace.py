"""An interface for the [MACE](https://github.com/ACEsuit/mace) calculator"""

import json

from .mace import MACE_driver

MACECalculator = None
all_changes = None
Calculator = None
torch = None

__DRIVER_NAME__ = "extmace"
__DRIVER_CLASS__ = "Extended_MACE_driver"


class Extended_MACE_driver(MACE_driver):
    """
    MACE driver with the torch tensors exposed.
    """

    def __init__(
        self, forward_kwargs: str = None, instructions: str = None, *args, **kwargs
    ):
        self.forward_kwargs = {}
        if forward_kwargs is not None:
            with open(forward_kwargs, "r") as f:
                self.forward_kwargs = json.load(f)

        self.instructions = {}
        if instructions is not None:
            with open(instructions, "r") as f:
                self.instructions = json.load(f)

        super().__init__(*args, **kwargs)

    def check_parameters(self):
        """Check the arguments requuired to run the driver

        This loads the potential and atoms template in MACE
        """

        super().check_parameters()

        self.ase_calculator = Extended_MACECalculator(
            model_paths=self.model,
            device=self.device,
            instructions=self.instructions,
            forward_kwargs=self.forward_kwargs,
            get_extras=lambda: self.extra,
            **self.mace_kwargs
        )


try:
    import torch
except:
    pass
    # raise ImportError("Couldn't load torch")

try:
    from mace.calculators import MACECalculator
except:
    pass
    # raise ImportError("Couldn't load mace bindings")

try:
    from ase.calculators.calculator import Calculator, all_changes
except:
    pass
    # raise ImportError("Couldn't load ase bindings")

if MACECalculator is not None:

    class Extended_MACECalculator(MACECalculator):
        def get_extras(self) -> dict:
            return {}

        def __init__(
            self,
            instructions: dict = {},
            forward_kwargs: dict = {},
            get_extras: callable = None,
            *argc,
            **kwargs
        ):
            if get_extras is not None:
                self.get_extras = get_extras
            self.forward_kwargs = forward_kwargs
            self.instructions = instructions
            super().__init__(*argc, **kwargs)

        def apply_ensemble(self, data: dict) -> dict:
            if self.instructions == {}:
                return data
            extras = self.get_extras()
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

            batch_base = self._atoms_to_batch(atoms)

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
                    training=self.use_compile,
                    compute_edge_forces=self.compute_atomic_stresses,
                    compute_atomic_stresses=self.compute_atomic_stresses,
                    **self.forward_kwargs
                )

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

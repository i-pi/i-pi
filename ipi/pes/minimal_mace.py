"""To be written"""

import numpy as np
from ipi.utils.messages import verbosity, warning

try:
    from .dummy import Dummy_driver
except:
    from dummy import Dummy_driver

try:
    from ase import Atoms
except ImportError as e:
    message = (
        "could not find 'ase.Atoms': make sure that 'ase' is installed."
    )
    warning(f"{message}: {e}")
    raise ImportError(message) from e

try:
    from ase.calculators.calculator import all_changes, all_properties
except ImportError as e:
    message = (
        "could not find 'all_changes' and 'all_properties' in 'ase': make sure that ase is installed."
    )
    warning(f"{message}: {e}")
    raise ImportError(message) from e

try:
    from mace.calculators import MACECalculator
except ImportError as e:
    message = (
        "could not find 'MACECalculator': make sure that 'MACE' is installed."
    )
    warning(f"{message}: {e}")
    raise ImportError(message) from e

__DRIVER_NAME__ = "minimal_mace"
__DRIVER_CLASS__ = "Minimal_MACE_driver"

# This driver inherits from MACECalculator so that it can overwrite the `calculate` method

class ExposedMACECalculator(MACECalculator):
    pass

class Minimal_MACE_driver(Dummy_driver):
    """
    To be written
    """

    def __init__(self, template, *args, **kwargs):
        warning(
            "THIS PES HAS NOT BEEN TESTED FOLLOWING CONVERSION TO THE NEW PES API.",
            verbosity.low,
        )

        try:
            from ase.outputs import _defineprop, all_outputs

            # avoid duplicate
            # it complains with a committee of ffdirect MACE models
            if "node_energy" not in all_outputs:
                _defineprop("node_energy", dtype=float, shape=("natoms",))
        except ImportError:
            raise ValueError("Could not find or import the ASE module")

        # set the device and dtype for torch
        # torch.set_default_device(device)
        # torch.set_default_dtype(getattr(torch, dtype))
        
        # it loads only once mode
        # self.model:torch.nn.Module = torch.load(model, map_location=device)
        # self.model.to(device)
        # self.model.eval()
        
        # initialize
        ASEDriver().__init__(self,template, *args, **kwargs)
        MACECalculator.__init__(self,*args, **kwargs)

    def check_parameters(self):
        super().check_parameters()
        self.ase_calculator = self
        
   

# pylint: disable=dangerous-default-value
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
        )
        if self.model_type in ["MACE", "EnergyDipoleMACE"]:
            ret_tensors["energies"][i] = out["energy"].detach()
            ret_tensors["node_energy"][i] = (out["node_energy"] - node_e0).detach()
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

        if out["stress"] is not None:
            self.results["stress"] = full_3x3_to_voigt_6_stress(
                torch.mean(ret_tensors["stress"], dim=0).cpu().numpy()
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
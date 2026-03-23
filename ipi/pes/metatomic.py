"""
Interface with metatomic ML models (https://docs.metatensor.org/metatomic/), that can be
used to perform calculations with arbitrary machine learning potentials.
"""

import json
import warnings
import concurrent.futures
import numpy as np

from ipi.utils.units import unit_to_internal, unit_to_user
from .dummy import Dummy_driver

from ipi.utils.messages import warning, info, verbosity

torch = None
mts = None
mta = None
vesin_metatomic = None
ase_io = None

__DRIVER_NAME__ = "metatomic"
__DRIVER_CLASS__ = "MetatomicDriver"

from time import time


class MetatomicDriver(Dummy_driver):
    """
    Driver for ``metatomic`` MLIPs.

    The driver requires the path to a torchscript model, and a template file that
    describes the chemical makeup of the structure. Requires the metatomic-torch
    library.

    Command-line:
        i-pi-py_driver -m metatomic -o template=template.xyz,model=model.json [...]

    Parameters:
        :param template: string, filename of an ASE-readable structure file to
            initialize atomic number and types
        :param model: string, filename of the torchscript model file
        :param device: string, optional, device the model should be run on
            can specify GPU ID, and specify multiple devices, which activates
            multiprocessing ["cpu" | "cuda", "cuda:0", "cuda:0;cuda:1", ...]
        :param extensions: string, optional, path to the compiled extensions for the
            model
        :param check_consistency: bool, optional, whether to perform various consistency
            checks when evaluating the model
        :param energy_ensemble: bool, optional, whether to compute an ensemble of
            energies for uncertainty quantification
        :param force_virial_ensemble: bool, optional, whether to compute an ensemble of
            forces and virials for uncertainty quantification (warning: this can be
            computationally expensive)
        :param non_conservative: bool, optional, whether to use non-conservative forces
            and stresses
        :param energy_variant: string, optional, which energy variant to pick for
            a multi-head model
        :param non_conservative_variant: string, optional, which energy variant to pick for
            a multi-head model (pick a different version for the direct predictions
        :param uncertainty_variant: string, optional, which energy variant to pick for
            a multi-head model (pick a different version than energy_variant for uncertainty and ensemble)
        :param uncertainty_threshold: float, optional, whether the model should
            estimate uncertainty and warn whenever it exceeds a selected threshold
            (in eV/atom)
    """

    def __init__(
        self,
        template,
        model,
        device="cpu",
        extensions="",
        check_consistency=False,
        energy_ensemble=False,
        force_virial_ensemble=False,
        non_conservative=False,
        energy_variant=None,
        non_conservative_variant=None,
        uncertainty_variant=None,
        uncertainty_threshold=0.0,
        *args,
        **kwargs,
    ):
        # Initialize dependencies
        global torch, mta, mts, vesin_metatomic, ase_io
        if torch is None or mta is None or mts is None or vesin_metatomic is None:
            try:
                import torch
                import metatensor.torch as mts
                import metatomic.torch as mta
                import vesin.metatomic as vesin_metatomic
                import ase.io as ase_io
            except ImportError as e:
                message = (
                    "could not find the metatomic driver dependencies, "
                    "make sure they are installed with "
                    "`python -m pip install --upgrade metatomic-torch vesin`"
                )
                warning(f"{message}: {e}")
                raise ImportError(message) from e

        self.model = model
        self.device = device
        self.extensions = extensions
        self.check_consistency = check_consistency
        self.energy_ensemble = energy_ensemble
        self.force_virial_ensemble = force_virial_ensemble
        self.non_conservative = non_conservative
        self.template = template

        # Determine suffixes for different model outputs
        self.energy_suffix = "" if energy_variant is None else f"/{energy_variant}"
        self.non_conservative_suffix = (
            self.energy_suffix
            if non_conservative_variant is None
            else f"/{non_conservative_variant}"
        )
        self.uncertainty_suffix = (
            self.energy_suffix
            if uncertainty_variant is None
            else f"/{uncertainty_variant}"
        )

        self.uncertainty_threshold = unit_to_internal(
            "energy", "electronvolt", uncertainty_threshold
        )
        super().__init__(*args, **kwargs)

        info(f"Model arguments:\n{args}\n{kwargs}", self.verbose)

    def check_parameters(self):
        """Check the arguments required to run the driver

        This loads the potential and atoms template and validates configuration.
        """
        self._check_metatomic_version()
        self._check_consistency_options()

        super().check_parameters()

        self._load_models_and_device()
        self._load_template()
        self._register_model_outputs()

        # Show the model metadata to the users
        info(f"Metatomic model information:\n{self.model.metadata()}", self.verbose)

    def _check_metatomic_version(self):
        """Ensures the installed metatomic version is compatible."""
        metatomic_major, metatomic_minor, *_ = mta.__version__.split(".")
        metatomic_major = int(metatomic_major)
        metatomic_minor = int(metatomic_minor)

        if metatomic_major != 0 or metatomic_minor != 1:
            raise ImportError(
                "this code is only compatible with metatomic-torch == v0.1, "
                f"found version v{mta.__version__} at '{mta.__file__}'"
            )

    def _check_consistency_options(self):
        """Checks for conflicting options."""
        if self.force_virial_ensemble and not self.energy_ensemble:
            raise ValueError(
                "force_virial_ensemble can only be set to True if "
                "energy_ensemble is also True"
            )

        if self.force_virial_ensemble and self.non_conservative:
            raise ValueError(
                "`force_virial_ensemble` and `non_conservative` cannot be set to True "
                "at the same time"
            )

    def _load_models_and_device(self):
        """Loads the model(s) onto the requested device(s)."""
        model_path = self.model

        if ";" in self.device:
            self.devices = [d.strip() for d in self.device.split(";")]
        else:
            self.devices = [self.device]

        self.models = []
        for device in self.devices:
            model = mta.load_atomistic_model(
                model_path, extensions_directory=self.extensions
            )
            model = model.to(device)
            self.models.append(model)

        # ensure self.device and self.modelare initialized (to the first device)
        # to ensure compatibility with single-device methods
        self.model = self.models[0]
        self.device = self.devices[0]

        self._dtype = getattr(torch, self.model.capabilities().dtype)

        # Initialize executor for multi-device execution
        if hasattr(self, "executor") and self.executor is not None:
            self.executor.shutdown()

        if len(self.devices) > 1:
            self.executor = concurrent.futures.ThreadPoolExecutor(
                max_workers=len(self.devices)
            )
        else:
            self.executor = None

    def _load_template(self):
        """Loads the atomic types from the template file."""
        atoms = ase_io.read(self.template)
        self._types = torch.from_numpy(atoms.numbers).to(dtype=torch.int32)

    def _register_model_outputs(self):
        """Registers the required outputs for the model evaluation."""
        outputs = {
            f"energy{self.energy_suffix}": mta.ModelOutput(
                quantity="energy",
                unit="eV",
                per_atom=False,
            )
        }

        if self.non_conservative:
            if "non_conservative_forces" not in self.model.capabilities().outputs:
                raise ValueError(
                    "Non-conservative evaluation was requested, but "
                    "this model does not support non-conservative forces. "
                )
            else:
                outputs[f"non_conservative_forces{self.non_conservative_suffix}"] = (
                    mta.ModelOutput(
                        quantity="force",
                        unit="eV/Angstrom",
                        per_atom=True,
                    )
                )
            if "non_conservative_stress" not in self.model.capabilities().outputs:
                warnings.warn(
                    "Non-conservative evaluation was requested, but "
                    "this model does not support non-conservative stresses. "
                    "Setting them to `NaN`; make sure your simulation does not require "
                    "them."
                )
            else:
                outputs[f"non_conservative_stress{self.non_conservative_suffix}"] = (
                    mta.ModelOutput(
                        quantity="pressure",
                        unit="eV/Angstrom^3",
                        per_atom=False,
                    )
                )

        if self.energy_ensemble or self.uncertainty_threshold > 0.0:
            outputs[f"energy_uncertainty{self.uncertainty_suffix}"] = mta.ModelOutput(
                quantity="energy",
                unit="eV",
                per_atom=True,
            )
            outputs[f"energy_ensemble{self.uncertainty_suffix}"] = mta.ModelOutput(
                quantity="energy",
                unit="eV",
                per_atom=False,
            )

        self.evaluation_options = mta.ModelEvaluationOptions(
            length_unit="Angstrom",
            outputs=outputs,
        )

    def _prepare_system(self, cell, pos, model=None, device=None):
        """Prepares a metatomic System object from cell and positions."""
        if model is None:
            model = self.model
        if device is None:
            device = self.device

        # Convert units to Angstrom
        positions = unit_to_user("length", "angstrom", pos)
        cell = unit_to_user("length", "angstrom", cell.T)

        positions = torch.from_numpy(positions).to(dtype=self._dtype)
        cell = torch.from_numpy(cell).to(dtype=self._dtype)
        pbc = torch.norm(cell, dim=1) != 0.0
        strain = None

        if not self.non_conservative:
            # Enable gradients for virial calculation
            positions.requires_grad_(True)
            strain = torch.eye(3, requires_grad=True, dtype=self._dtype)
            positions = positions @ strain
            cell = cell @ strain

        system = mta.System(self._types, positions, cell, pbc)
        system = system.to(device)

        # Compute neighbor lists using vesin
        vesin_metatomic.compute_requested_neighbors(
            system, system_length_unit="Angstrom", model=model
        )

        return system, strain

    def _process_standard_outputs(self, outputs, systems, strains, device):
        """Extracts and computes energy, forces, and virials."""
        num_systems = len(systems)
        energy_tensor = outputs[f"energy{self.energy_suffix}"].block().values

        if self.non_conservative:
            # For non-conservative, forces are direct outputs
            forces_tensor = (
                outputs[f"non_conservative_forces{self.non_conservative_suffix}"]
                .block()
                .values
            )
            if "non_conservative_stress" in outputs:
                stresses_tensor = (
                    outputs[f"non_conservative_stress{self.non_conservative_suffix}"]
                    .block()
                    .values
                )
            else:
                stresses_tensor = torch.full(
                    (num_systems, 3, 3),
                    torch.nan,
                    device=device,
                    dtype=self._dtype,
                )
            # Virial from stress
            virials_tensor = torch.stack(
                [
                    -stress_tensor * torch.abs(torch.det(system.cell))
                    for stress_tensor, system in zip(stresses_tensor, systems)
                ]
            )
        else:
            # For conservative, compute gradients (forces and virials)
            gradients = torch.autograd.grad(
                energy_tensor,
                [system.positions for system in systems] + strains,
                grad_outputs=-torch.ones_like(energy_tensor),
            )
            # First half are forces, second half are virials
            forces_tensors, virial_tensors = (
                gradients[: len(systems)],
                gradients[len(systems) :],
            )
            forces_tensor = torch.stack(forces_tensors)
            virials_tensor = torch.stack(virial_tensors)

        # Convert units back to internal
        energies = unit_to_internal(
            "energy",
            "electronvolt",
            energy_tensor.detach()
            .to(device="cpu", dtype=torch.float64)
            .numpy()
            .reshape(num_systems, 1),
        )
        forces = unit_to_internal(
            "force",
            "ev/ang",
            forces_tensor.detach()
            .reshape(num_systems, -1, 3)
            .to(device="cpu", dtype=torch.float64)
            .numpy(),
        )
        virials = unit_to_internal(
            "energy",
            "electronvolt",
            virials_tensor.detach()
            .reshape(num_systems, 3, 3)
            .to(device="cpu", dtype=torch.float64)
            .numpy(),
        )
        # Symmetrize virial
        virials = (virials + virials.swapaxes(1, 2)) * 0.5

        return energies, forces, virials

    def _process_uncertainty(self, outputs, extras_dicts):
        """Extracts uncertainty information if requested."""
        if self.uncertainty_threshold <= 0.0:
            return

        num_systems = len(extras_dicts)
        energy_uncertainty_tensor = (
            outputs[f"energy_uncertainty{self.uncertainty_suffix}"].block().values
        )

        energy_uncertainty = unit_to_internal(
            "energy",
            "electronvolt",
            energy_uncertainty_tensor.detach()
            .to(device="cpu", dtype=torch.float64)
            .numpy(),
        ).reshape(num_systems, -1)

        for extras_dict, energy_uq in zip(extras_dicts, energy_uncertainty):
            extras_dict["energy_uncertainty"] = list(energy_uq)
            if np.max(energy_uq) > self.uncertainty_threshold:
                warning(
                    f"the estimated atomic energy uncertainty {unit_to_user('energy','electronvolt',np.max(energy_uq)):.4f} eV/atom "
                    f"exceeds the selected threshold of {unit_to_user('energy','electronvolt',self.uncertainty_threshold):.4f} eV/atom",
                    verbosity.medium,
                )

    def _compute_force_virial_ensemble(self, systems, strains, model, device):
        """Computes force and virial ensembles using Jacobian calculation."""

        # Function wrapper for jacobian calculation
        def _compute_ensemble(*positions_and_strain_tuple):
            positions_list = positions_and_strain_tuple[: len(systems)]
            strain_list = positions_and_strain_tuple[len(systems) :]
            new_systems = [
                mta.System(
                    system.types,
                    positions @ strain,
                    system.cell @ strain,
                    system.pbc,
                )
                for system, positions, strain in zip(
                    systems, positions_list, strain_list
                )
            ]
            for options in model.requested_neighbor_lists():
                # Re-attach pre-computed neighbor lists
                for system, new_system in zip(systems, new_systems):
                    neighbors = mts.detach_block(system.get_neighbor_list(options))
                    mta.register_autograd_neighbors(
                        new_system,
                        neighbors,
                        check_consistency=self.check_consistency,
                    )
                    new_system.add_neighbor_list(options, neighbors)

            return (
                model(
                    new_systems,
                    self.evaluation_options,
                    check_consistency=self.check_consistency,
                )[f"energy_ensemble{self.uncertainty_suffix}"]
                .block()
                .values.sum(0)
            )

        # Compute Jacobian
        # Note: vectorize=False is required to avoid mysterious runtime errors
        # (potential PyTorch/Metatomic bug workaround)
        jacobians = torch.autograd.functional.jacobian(
            _compute_ensemble,
            tuple(
                [system.positions for system in systems]
                + [s.to(device) for s in strains]
            ),
            vectorize=False,
        )

        force_ensemble_tensors = [-j for j in jacobians[: len(systems)]]
        virial_ensemble_tensors = [-j for j in jacobians[len(systems) :]]

        return force_ensemble_tensors, virial_ensemble_tensors

    def _process_ensemble(self, outputs, systems, strains, model, device, extras_dicts):
        """Extracts ensemble energy and optionally force/virial ensembles."""
        if not self.energy_ensemble:
            return

        energy_ensembles_tensor = (
            outputs[f"energy_ensemble{self.uncertainty_suffix}"].block().values
        )
        energy_ensembles = unit_to_internal(
            "energy",
            "electronvolt",
            energy_ensembles_tensor.detach()
            .to(device="cpu", dtype=torch.float64)
            .numpy(),
        )
        for extras_dict, energy_ensemble in zip(extras_dicts, energy_ensembles):
            extras_dict["committee_pot"] = list(energy_ensemble)

        if self.force_virial_ensemble:
            force_ensemble_tensors, virial_ensemble_tensors = (
                self._compute_force_virial_ensemble(systems, strains, model, device)
            )

            force_ensembles = [
                unit_to_internal(
                    "force",
                    "ev/ang",
                    force_ensemble_tensor.detach()
                    .to(device="cpu", dtype=torch.float64)
                    .numpy(),
                )
                for force_ensemble_tensor in force_ensemble_tensors
            ]
            virial_ensembles = [
                unit_to_internal(
                    "energy",
                    "electronvolt",
                    virial_ensemble_tensor.detach()
                    .to(device="cpu", dtype=torch.float64)
                    .numpy(),
                )
                for virial_ensemble_tensor in virial_ensemble_tensors
            ]

            # Symmetrize virial ensembles
            virial_ensembles = [
                (virial_ensemble + virial_ensemble.swapaxes(1, 2)) / 2.0
                for virial_ensemble in virial_ensembles
            ]

            for extras_dict, force_ensemble, virial_ensemble in zip(
                extras_dicts, force_ensembles, virial_ensembles
            ):
                extras_dict["committee_force"] = list(force_ensemble.flatten())
                extras_dict["committee_virial"] = list(virial_ensemble.flatten())

    def _process_outputs(self, outputs, systems, strains, model=None):
        """Processes raw model outputs into i-PI compatible formats."""
        if model is None:
            model = self.model

        device = systems[0].positions.device
        num_systems = len(systems)

        # 1. Standard Energy, Force, Virial
        energies, forces, virials = self._process_standard_outputs(
            outputs, systems, strains, device
        )

        extras_dicts = [{} for _ in range(num_systems)]

        # 2. Uncertainty
        self._process_uncertainty(outputs, extras_dicts)

        # 3. Ensembles
        self._process_ensemble(outputs, systems, strains, model, device, extras_dicts)

        extras_list = [json.dumps(extras_dict) for extras_dict in extras_dicts]

        energies = energies.flatten()
        forces = forces.reshape((len(systems), -1, 3))
        virials = virials.reshape((len(systems), 3, 3))

        # Remove net force in the non-conservative case
        if self.non_conservative:
            forces = forces - forces.mean(1, keepdims=True)

        return [
            (e, f, v, extras)
            for (e, f, v, extras) in zip(energies, forces, virials, extras_list)
        ]

    def compute_structure(self, cell, pos):
        """Computes properties for a single structure."""
        system, strain = self._prepare_system(cell, pos)

        outputs = self.model(
            [system],
            self.evaluation_options,
            check_consistency=self.check_consistency,
        )

        return self._process_outputs(outputs, [system], [strain])[0]

    def _compute_chunk(self, device_idx, cell_chunk, pos_chunk):
        """Helper to compute a batch chunk on a specific device."""
        model = self.models[device_idx]
        device = self.devices[device_idx]

        sys_batch = []
        strain_batch = []

        # Prepare systems on the target device
        for c, p in zip(cell_chunk, pos_chunk):
            sy, st = self._prepare_system(c, p, model=model, device=device)
            sys_batch.append(sy)
            strain_batch.append(st)

        if torch.cuda.is_available() and "cuda" in str(device):
            torch.cuda.synchronize(device=device)

        # Run model
        outputs = model(
            sys_batch,
            self.evaluation_options,
            check_consistency=self.check_consistency,
        )

        # Process outputs
        return self._process_outputs(outputs, sys_batch, strain_batch, model=model)

    def compute(self, cell, pos):
        """Calls the model evaluation, taking care of both serial and batched execution.

        The implementation handles both batched evaluations (many structures are sent at
        once for evaluation by one instance) and concurrent evaluation (the batch is
        split across multiple devices, each running its own instance of the model).
        """

        if isinstance(cell, list):
            # Batched mode
            if not isinstance(pos, list) or len(cell) != len(pos):
                raise ValueError(
                    "Both position and cell should be given as lists to run in batched mode"
                )

            n_devices = len(self.devices)
            n_samples = len(cell)
            pre_time = -time()

            if n_devices == 1:
                processed_outputs = self._compute_chunk(0, cell, pos)
            else:
                # Multi-device batch splitting
                batch_size = (n_samples + n_devices - 1) // n_devices
                futures = []

                for i in range(n_devices):
                    start = i * batch_size
                    end = min(start + batch_size, n_samples)
                    if start >= end:
                        break

                    cell_chunk = cell[start:end]
                    pos_chunk = pos[start:end]
                    futures.append(
                        self.executor.submit(
                            self._compute_chunk, i, cell_chunk, pos_chunk
                        )
                    )

                processed_outputs = []
                for f in futures:
                    processed_outputs.extend(f.result())

            pre_time += time()

            info(
                "Total compute time %e" % (pre_time),
                verbosity.high,
            )

            return processed_outputs

        else:
            # Single structure mode
            return self.compute_structure(cell, pos)

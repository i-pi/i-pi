"""
Interface with metatensor
(https://lab-cosmo.github.io/metatensor/latest/atomistic/index.html), that can
be used to perform calculations based on different types of machine learning
potentials
"""

import json
import warnings

from ipi.utils.units import unit_to_internal, unit_to_user
from .dummy import Dummy_driver

from ipi.utils.messages import warning, info
import ase.io

torch = None
mts = None
mta = None
vesin_torch_metatensor = None

__DRIVER_NAME__ = "metatensor"
__DRIVER_CLASS__ = "MetatensorDriver"


class MetatensorDriver(Dummy_driver):
    """
    Driver for `metatensor` MLIPs
    The driver requires specification of a torchscript model,
    and a template file that describes the chemical makeup
    of the structure. Requires the metatensor-torch library

    Command-line:
        i-pi-py_driver -m metatensor -o template=template.xyz,model=model.json [...]

    Parameters:
        :param template: string, filename of an ASE-readable structure file
            to initialize atomic number and types
        :param model: string, filename of the torchscript model file
        :param device: string, optional, ["cpu" | "cuda"]
        :param extensions: string, optional, path to the compiled extensions for the
            model
        :param check_consistency: bool, optional, whether to perform various consistency
            checks when evaluating the model
        :param energy_ensemble: bool, optional, whether to compute an ensemble of
            energies for uncertainty quantification
        :param force_virial_ensemble: bool, optional, whether to compute an ensemble
            of forces and virials for uncertainty quantification (warning: this can be
            computationally expensive)
        :param non_conservative: bool, optional, whether to use non-conservative forces
            and stresses
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
        *args,
        **kwargs,
    ):
        global torch, mts, mta, vesin_torch_metatensor
        if (
            torch is None
            or mts is None
            or mta is None
            or vesin_torch_metatensor is None
        ):
            try:
                import torch
                import metatensor.torch as mts
                import metatensor.torch.atomistic as mta
                import vesin.torch.metatensor as vesin_torch_metatensor
            except ImportError as e:
                message = (
                    "could not find the metatensor driver dependencies, "
                    "make sure they are installed with "
                    "`python -m pip install metatensor[torch] vesin[torch]`"
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
        super().__init__(*args, **kwargs)

        info(f"Model arguments:\n{args}\n{kwargs}", self.verbose)

    def check_parameters(self):
        """Check the arguments required to run the driver

        This loads the potential and atoms template in metatensor
        """

        metatensor_major, metatensor_minor, *_ = mts.__version__.split(".")
        metatensor_major = int(metatensor_major)
        metatensor_minor = int(metatensor_minor)

        if metatensor_major != 0 or metatensor_minor < 6:
            raise ImportError(
                "this code is only compatible with metatensor-torch >= v0.6, "
                f"found version v{mts.__version__} "
                f"at '{mts.__file__}'"
            )

        if self.force_virial_ensemble and not self.energy_ensemble:
            raise ValueError(
                "force_virial_ensemble can only be set to True if "
                "energy_ensemble is also True"
            )

        # these two are not compatible, since force and stress ensembles are calculated
        # by autograd on the memebers of the energy ensemble
        if self.force_virial_ensemble and self.non_conservative:
            raise ValueError(
                "`force_virial_ensemble` and `non_conservative` cannot be set to True "
                "at the same time"
            )

        super().check_parameters()

        # Load the model
        model_path = self.model
        self.model = mta.load_atomistic_model(
            model_path, extensions_directory=self.extensions
        )
        self.model = self.model.to(self.device)
        self._dtype = getattr(torch, self.model.capabilities().dtype)

        # read the template and extract the corresponding atomic types
        atoms = ase.io.read(self.template)
        self._types = torch.from_numpy(atoms.numbers).to(dtype=torch.int32)

        # Register the requested outputs
        outputs = {"energy": mta.ModelOutput(quantity="energy", unit="eV")}
        if self.non_conservative:
            if "non_conservative_forces" not in self.model.capabilities().outputs:
                raise ValueError(
                    "Non-conservative evaluation was requested, but "
                    "this model does not support non-conservative forces. "
                )
            else:
                outputs["non_conservative_forces"] = mta.ModelOutput(
                    quantity="force", unit="eV/Angstrom", per_atom=True
                )
            if "non_conservative_stress" not in self.model.capabilities().outputs:
                warnings.warn(
                    "Non-conservative evaluation was requested, but "
                    "this model does not support non-conservative stresses. "
                    "Setting them to `nan`; make sure your simulation does not require "
                    "them."
                )
            else:
                outputs["non_conservative_stress"] = mta.ModelOutput(
                    quantity="pressure", unit="eV/Angstrom^3"
                )
        if self.energy_ensemble:
            outputs["energy_ensemble"] = mta.ModelOutput(quantity="energy", unit="eV")
        self.evaluation_options = mta.ModelEvaluationOptions(
            length_unit="Angstrom",
            outputs=outputs,
        )

        # Show the model metadata to the users
        info(f"Metatomic model data:\n{self.model.metadata()}", self.verbose)

    def __call__(self, cell, pos):
        """Get energies, forces and virials from the atomistic model."""
        positions = unit_to_user("length", "angstrom", pos)
        cell = unit_to_user("length", "angstrom", cell.T)

        positions = torch.from_numpy(positions).to(dtype=self._dtype)
        cell = torch.from_numpy(cell).to(dtype=self._dtype)
        pbc = torch.norm(cell, dim=1) != 0.0

        if not self.non_conservative:
            positions.requires_grad_(True)
            # this is to compute the virial (which is related to the derivative with
            # respect to the strain)
            strain = torch.eye(3, requires_grad=True, dtype=self._dtype)
            positions = positions @ strain
            cell = cell @ strain

        system = mta.System(self._types, positions, cell, pbc)
        # compute the requires neighbor lists with vesin
        vesin_torch_metatensor.compute_requested_neighbors(
            system, system_length_unit="A", model=self.model
        )
        system = system.to(self.device)

        outputs = self.model(
            [system],
            self.evaluation_options,
            check_consistency=self.check_consistency,
        )
        energy_tensor = outputs["energy"].block().values

        if self.non_conservative:
            forces_tensor = outputs["non_conservative_forces"].block().values
            if "non_conservative_stress" in outputs:
                stress_tensor = outputs["non_conservative_stress"].block().values
            else:
                stress_tensor = torch.full(
                    (3, 3), torch.nan, device=self.device, dtype=self._dtype
                )
            virial_tensor = -stress_tensor * torch.abs(torch.det(cell))
        else:
            forces_tensor, virial_tensor = torch.autograd.grad(
                energy_tensor,
                (positions, strain),
                grad_outputs=-torch.ones_like(energy_tensor),
            )

        energy = unit_to_internal(
            "energy",
            "electronvolt",
            energy_tensor.detach().to(device="cpu", dtype=torch.float64).numpy().item(),
        )
        forces = unit_to_internal(
            "force",
            "ev/ang",
            forces_tensor.detach()
            .reshape(3, -1)
            .to(device="cpu", dtype=torch.float64)
            .numpy(),
        )
        virial = unit_to_internal(
            "energy",
            "electronvolt",
            virial_tensor.detach()
            .reshape(3, 3)
            .to(device="cpu", dtype=torch.float64)
            .numpy(),
        )
        # symmetrize the virial for rotationally unconstrained models
        virial = (virial + virial.T) / 2.0

        extras_dict = {}

        if self.energy_ensemble:
            energy_ensemble_tensor = (
                outputs["energy_ensemble"].block().values.squeeze(0)
            )
            energy_ensemble = unit_to_internal(
                "energy",
                "electronvolt",
                energy_ensemble_tensor.detach()
                .to(device="cpu", dtype=torch.float64)
                .numpy(),
            )
            extras_dict["committee_pot"] = list(energy_ensemble)

            # the ensemble of forces and virials require a more expensive calculation
            # so they are controlled by a separate option
            if self.force_virial_ensemble:
                # this function definition is necessary to use
                # torch.autograd.functional.jacobian, which is vectorized
                def _compute_ensemble(positions, strain):
                    new_system = mta.System(
                        self._types, positions @ strain, cell @ strain, pbc
                    )
                    for options in self.model.requested_neighbor_lists():
                        # get the pre-computed NL, and re-attach it in the compute graph
                        neighbors = mts.detach_block(system.get_neighbor_list(options))
                        mta.register_autograd_neighbors(
                            new_system,
                            neighbors,
                            check_consistency=self.check_consistency,
                        )
                        new_system.add_neighbor_list(options, neighbors)

                    return (
                        self.model(
                            [new_system],
                            self.evaluation_options,
                            check_consistency=self.check_consistency,
                        )["energy_ensemble"]
                        .block()
                        .values.squeeze(0)
                    )

                minus_force_ensemble_tensor, minus_virial_ensemble_tensor = (
                    torch.autograd.functional.jacobian(
                        _compute_ensemble, (positions, strain), vectorize=True
                    )
                )
                force_ensemble_tensor = -minus_force_ensemble_tensor
                virial_ensemble_tensor = -minus_virial_ensemble_tensor
                force_ensemble = unit_to_internal(
                    "force",
                    "ev/ang",
                    force_ensemble_tensor.detach()
                    .to(device="cpu", dtype=torch.float64)
                    .numpy(),
                )
                virial_ensemble = unit_to_internal(
                    "energy",
                    "electronvolt",
                    virial_ensemble_tensor.detach()
                    .to(device="cpu", dtype=torch.float64)
                    .numpy(),
                )
                # symmetrize the virial for rotationally unconstrained models
                virial_ensemble = (
                    virial_ensemble + virial_ensemble.swapaxes(1, 2)
                ) / 2.0
                extras_dict["committee_force"] = list(force_ensemble.flatten())
                extras_dict["committee_virial"] = list(virial_ensemble.flatten())

        extras = json.dumps(extras_dict)
        return energy, forces, virial, extras

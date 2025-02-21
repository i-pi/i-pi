"""
Interface with metatensor
(https://lab-cosmo.github.io/metatensor/latest/atomistic/index.html), that can
be used to perform calculations based on different types of machine learning
potentials
"""

import json

from .dummy import Dummy_driver

from ipi.utils.messages import warning

torch = None
mtt = None
mta = None

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
        *args,
        **kwargs,
    ):
        global torch, mtt, mta
        if torch is None and mtt is None and mta is None:
            try:
                import torch
                import metatensor.torch as mtt
                import metatensor.torch.atomistic as mta
            except ImportError as e:
                warning(f"Could not find or import metatensor.torch: {e}")

        self.model = model
        self.device = device
        self.extensions = extensions
        self.check_consistency = check_consistency
        self.energy_ensemble = energy_ensemble
        self.force_virial_ensemble = force_virial_ensemble
        print("ARGS ", kwargs)
        super().__init__(template, *args, **kwargs)

    def check_parameters(self):
        """Check the arguments required to run the driver

        This loads the potential and atoms template in metatensor
        """

        metatensor_major, metatensor_minor, *_ = mtt.__version__.split(".")
        metatensor_major = int(metatensor_major)
        metatensor_minor = int(metatensor_minor)

        if metatensor_major != 0 or metatensor_minor < 6:
            raise ImportError(
                "this code is only compatible with metatensor-torch >= v0.6, "
                f"found version v{mtt.__version__} "
                f"at '{mtt.__file__}'"
            )

        if self.force_virial_ensemble and not self.energy_ensemble:
            raise ValueError(
                "force_virial_ensemble can only be set to True if "
                "energy_ensemble is also True"
            )

        super().check_parameters()

        # Load the model
        model_path = self.model
        self.model = mta.load_atomistic_model(
            model_path, extensions_directory=self.extensions
        )
        self.model = self.model.to(self.device)

        # Register the requested outputs
        outputs = {"energy": mta.ModelOutput(quantity="energy", unit="Hartree")}
        if self.energy_ensemble:
            outputs["energy_ensemble"] = mta.ModelOutput(
                quantity="energy", unit="Hartree"
            )
        self.evaluation_options = mta.ModelEvaluationOptions(
            length_unit="Bohr",
            outputs=outputs,
        )

        # Show the model metadata to the users
        print(self.model.metadata())

    def __call__(self, cell, pos):
        """Get energies, forces and virials from the atomistic model."""

        ase_atoms = self.template_ase.copy()
        ase_atoms.positions = pos
        ase_atoms.cell = cell

        types, positions, cell, pbc = mta.ase_calculator._ase_to_torch_data(
            atoms=ase_atoms, dtype=torch.float64, device=self.device
        )
        positions.requires_grad_(True)
        strain = torch.eye(
            3, requires_grad=True, device=self._device, dtype=self._dtype
        )
        positions = positions @ strain
        cell = cell @ strain
        system = mta.System(types, positions, cell, pbc)

        for options in self.model.requested_neighbor_lists():
            neighbors = mta.ase_calculator._compute_ase_neighbors(
                ase_atoms, options, dtype=self._dtype, device=self._device
            )
            mta.register_autograd_neighbors(
                system,
                neighbors,
                check_consistency=self.check_consistency,
            )
            system.add_neighbor_list(options, neighbors)

        outputs = self._model(
            [system],
            self.evaluation_options,
            check_consistency=self.check_consistency,
        )
        energy_tensor = outputs["energy"].block().values
        forces_tensor, virial_tensor = torch.autograd.grad(
            energy_tensor,
            (positions, strain),
            grad_outputs=-torch.ones_like(energy_tensor),
        )

        energy = energy_tensor.detach().cpu().numpy()
        forces = forces_tensor.detach().cpu().numpy()
        virial = virial_tensor.detach().cpu().numpy()

        extras_dict = {}

        if self.energy_ensemble:
            energy_ensemble_tensor = (
                outputs["energy_ensemble"].block().values.squeeze(0)
            )
            energy_ensemble = energy_ensemble_tensor.detach().cpu().numpy()
            extras_dict["committee_pot"] = list(energy_ensemble)

            if self.force_virial_ensemble:
                # this function definition is necessary to use
                # torch.autograd.functional.jacobian, which is vectorized
                def _compute_ensemble(positions, strain):
                    return (
                        self._model(
                            [mta.System(types, positions @ strain, cell @ strain, pbc)],
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
                force_ensemble = force_ensemble_tensor.detach().cpu().numpy()
                virial_ensemble = virial_ensemble_tensor.detach().cpu().numpy()
                extras_dict["committee_force"] = list(force_ensemble.flatten())
                extras_dict["committee_virial"] = list(virial_ensemble.flatten())

        extras = json.dumps(extras_dict)
        return energy, forces, virial, extras

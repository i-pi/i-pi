"""
Interface with metatomic-core models (https://docs.metatensor.org/metatomic/), i.e.
any model implementing the metatomic C API: TorchScript models loaded through the
``metatomic-torch`` plugin, JAX/NumPy models defined in Python, or native plugins.

Unlike the ``metatomic`` driver (TorchScript only, through ``metatomic.torch``),
this driver only depends on the backend-independent ``metatomic`` (metatomic-core)
Python package, ``metatensor`` and ``vesin``.
"""

import importlib.util
import pathlib
import sys

import numpy as np

from ipi.utils.messages import info, warning
from ipi.utils.units import unit_to_internal, unit_to_user
from .dummy import Dummy_driver

mta = None
metatensor = None
vesin = None

__DRIVER_NAME__ = "metatomic_core"
__DRIVER_CLASS__ = "MetatomicCoreDriver"


def _as_bool(value):
    return value if isinstance(value, bool) else str(value).lower() in ("1", "true")


class MetatomicCoreDriver(Dummy_driver):
    """
    Driver for models implementing the metatomic C API (metatomic-core).

    Command-line:
        i-pi-py_driver -m metatomic_core -o template=template.xyz,model=<model> [...]

    Parameters:
        :param template: string, ASE-readable structure file giving the atomic types
            and the periodic boundary conditions
        :param model: string, the model to run. Either ``path/to/file.py:name``, where
            ``name`` is a ``metatomic.BaseModel`` subclass (instantiated without
            arguments), a function returning a model, or a model instance; or any
            string accepted by ``metatomic.load_model`` (e.g. a ``.pt`` file for the
            metatomic-torch plugin, or the model name of a native plugin)
        :param plugins: string, optional, ``;``-separated list of plugin shared
            libraries to load before loading the model. Bare names are looked up in
            metatomic's ``libexec/metatomic`` directory, and ``torch`` stands for the
            metatomic-torch plugin. The torch plugin is loaded automatically for
            ``.pt`` models.
        :param plugin_name: string, optional, the plugin ``metatomic.load_model``
            should use
        :param option_<key>: optional, string option ``key`` given to
            ``metatomic.load_model`` (e.g. ``option_cutoff=7.0``)
        :param energy_variant: string, optional, which ``energy/<variant>`` output
            to use for a multi-head model
        :param check_consistency: bool, optional, whether metatomic should check the
            data given to and returned by the model
        :param enable_x64: bool, optional (default true), enable 64-bit floats in JAX
            before the model is constructed (JAX models cast their parameters to
            ``float64`` on construction)
        :param skin: float, optional, neighbor list skin in Angstrom (default 2.0)
    """

    def __init__(
        self,
        template,
        model,
        plugins=None,
        plugin_name=None,
        energy_variant=None,
        check_consistency=False,
        enable_x64=True,
        skin=2.0,
        *args,
        **kwargs,
    ):
        global mta, metatensor, vesin
        try:
            import metatomic as mta
            import metatensor
            import vesin
        except ImportError as e:
            message = (
                "could not find the metatomic_core driver dependencies, make sure "
                "they are installed with `python -m pip install metatomic vesin ase`"
            )
            warning(f"{message}: {e}")
            raise ImportError(message) from e

        self.template = template
        self.model_spec = str(model)
        self.plugins = [p for p in str(plugins or "").split(";") if p.strip()]
        self.plugin_name = plugin_name
        self.model_options = {
            key[len("option_") :]: str(kwargs.pop(key))
            for key in list(kwargs)
            if key.startswith("option_")
        }
        self.energy_variant = energy_variant
        self.check_consistency = _as_bool(check_consistency)
        self.enable_x64 = _as_bool(enable_x64)
        self.skin = float(skin)
        super().__init__(*args, **kwargs)

    def check_parameters(self):
        import ase.io

        atoms = ase.io.read(self.template)
        self._types = atoms.numbers.astype(np.int32)
        self._pbc = np.array(atoms.pbc, dtype=bool)

        # outputs of native plugins use the C++ metatensor array type, older
        # metatensor releases need to be told how to read them from Python
        for origin in (
            "metatensor::SimpleDataArray",
            "metatensor_torch::TorchDataArray",
        ):
            try:
                metatensor.register_external_data_wrapper(
                    origin, metatensor.ExternalCpuArray
                )
            except Exception:
                pass

        self.model = self._load_model()
        capabilities = self.model.capabilities()
        self._dtype = np.dtype(capabilities.dtype)
        self._energy = self._energy_request(capabilities)
        self._pair_lists = self.model.requested_pair_lists([self._energy])
        inputs = self.model.requested_inputs([self._energy])
        if len(inputs) != 0:
            raise ValueError(
                "the metatomic_core driver can not provide the additional inputs "
                f"requested by this model: {[q.name for q in inputs]}"
            )
        length = mta.unit_conversion_factor(capabilities.length_unit, "angstrom")
        self._neighbors = [
            vesin.NeighborList(
                cutoff=options.cutoff * length,
                full_list=options.full_list,
                skin=self.skin,
            )
            for options in self._pair_lists
        ]
        unknown = set(self._types.tolist()) - set(capabilities.atomic_types)
        if unknown:
            raise ValueError(f"the model does not support atomic types {unknown}")
        info(f"metatomic-core model:\n{self.model.metadata()}", self.verbose)
        info(f"computing {self._energy}", self.verbose)

    def _load_model(self):
        if self.enable_x64:
            try:  # must happen before any JAX array is created by the model
                import jax

                jax.config.update("jax_enable_x64", True)
            except ImportError:
                pass

        path, _, name = self.model_spec.rpartition(":")
        if path.endswith(".py"):
            path = pathlib.Path(path).resolve()
            sys.path.insert(0, str(path.parent))
            spec = importlib.util.spec_from_file_location(path.stem, path)
            module = importlib.util.module_from_spec(spec)
            sys.modules[path.stem] = module
            spec.loader.exec_module(module)
            model = getattr(module, name)
            if not isinstance(model, mta.BaseModel):
                model = model()
            if not isinstance(model, mta.BaseModel):
                raise TypeError(f"{self.model_spec} is not a metatomic.BaseModel")
            return model

        plugins = list(self.plugins)
        if self.model_spec.endswith(".pt") and "torch" not in plugins:
            plugins.append("torch")
        for plugin in plugins:
            mta.load_plugin(self._plugin_path(plugin.strip()))
        return mta.load_model(
            self.model_spec,
            options=self.model_options or None,
            plugin_name=self.plugin_name,
        )

    @staticmethod
    def _plugin_path(plugin):
        if plugin == "torch":
            import metatomic_torch

            return metatomic_torch.utils.plugin_path
        if pathlib.Path(plugin).exists():
            return plugin
        libexec = pathlib.Path(mta.__file__).parent / "libexec" / "metatomic"
        candidates = [libexec / plugin, libexec / f"{plugin}.so"]
        candidates += [libexec / f"metatomic-{plugin}-plugin.so"]
        for candidate in candidates:
            if candidate.exists():
                return str(candidate)
        raise ValueError(f"could not find metatomic plugin '{plugin}'")

    def _energy_request(self, capabilities):
        name = (
            "energy" if self.energy_variant is None else f"energy/{self.energy_variant}"
        )
        candidates = [q for q in capabilities.outputs if q.name == name]
        if not candidates and self.energy_variant is None:
            # a model with a single, variant-only energy output
            candidates = [
                q for q in capabilities.outputs if q.name.startswith("energy/")
            ]
            if len({q.name for q in candidates}) > 1:
                raise ValueError(
                    "the model has several energy variants, select one with "
                    f"energy_variant: {sorted({q.name for q in candidates})}"
                )
        if not candidates:
            raise ValueError(f"the model does not provide a '{name}' output")
        # prefer per-system energies, summed per-atom energies otherwise
        declared = sorted(candidates, key=lambda q: q.sample_kind != "system")[0]
        if "positions" not in declared.gradients:
            raise ValueError(f"'{declared.name}' has no gradients w.r.t. positions")
        gradients = [g for g in ("positions", "strain") if g in declared.gradients]
        if "strain" not in gradients:
            warning("the model has no strain gradients: the virial will be zero")
        return mta.Quantity(
            name=declared.name,
            unit="eV",
            sample_kind=declared.sample_kind,
            gradients=gradients,
        )

    def _system(self, cell, pos):
        positions = unit_to_user("length", "angstrom", pos).reshape(-1, 3)
        # i-PI stores cell vectors as columns, metatomic as rows
        cell = unit_to_user("length", "angstrom", cell.T).copy()
        cell[~self._pbc] = 0.0
        if len(positions) != len(self._types):
            raise ValueError(
                f"got {len(positions)} atoms, the template has {len(self._types)}"
            )
        system = mta.System(
            "angstrom",
            self._types.copy(),
            np.ascontiguousarray(positions, dtype=self._dtype),
            np.ascontiguousarray(cell, dtype=self._dtype),
            self._pbc.copy(),
        )
        for options, neighbors in zip(self._pair_lists, self._neighbors):
            i, j, shifts, vectors = neighbors.compute(
                points=positions,
                box=cell,
                periodic=self._pbc,
                quantities="ijSD",
                copy=True,
            )
            samples = metatensor.Labels(
                [
                    "first_atom",
                    "second_atom",
                    "cell_shift_a",
                    "cell_shift_b",
                    "cell_shift_c",
                ],
                np.column_stack([i, j, shifts]).astype(np.int32).reshape(-1, 5),
            )
            pairs = metatensor.TensorBlock(
                np.ascontiguousarray(vectors, dtype=self._dtype).reshape(-1, 3, 1),
                samples,
                [metatensor.Labels.range("xyz", 3)],
                metatensor.Labels.range("distance", 1),
            )
            system.add_pairs(options, pairs)
        return system

    def compute(self, cell, pos):
        """Evaluates one structure, or a list of structures in a single model call"""
        batched = isinstance(cell, list)
        cells, positions = (cell, pos) if batched else ([cell], [pos])
        if len(cells) != len(positions):
            raise ValueError("cell and pos lists should have the same length")

        systems = [self._system(c, p) for c, p in zip(cells, positions)]
        (energy,) = mta.execute_model(
            self.model,
            systems,
            [self._energy],
            check_consistency=self.check_consistency,
        )
        block = energy.block().to(arrays="numpy")
        n_systems = len(systems)

        energies = np.zeros(n_systems)
        np.add.at(energies, block.samples.column("system"), block.values[:, 0])

        gradient = block.gradient("positions")
        forces = [np.zeros((len(self._types), 3)) for _ in systems]
        for s, a, value in zip(
            gradient.samples.column("system"),
            gradient.samples.column("atom"),
            np.asarray(gradient.values)[..., 0],
        ):
            forces[s][a] -= value

        virials = np.zeros((n_systems, 3, 3))
        if "strain" in block.gradients_list():
            gradient = block.gradient("strain")
            parent = block.samples.column("system")[gradient.samples.column("sample")]
            np.add.at(virials, parent, -np.asarray(gradient.values)[..., 0])

        results = [
            (
                unit_to_internal("energy", "electronvolt", e),
                unit_to_internal("force", "ev/ang", f).reshape(np.shape(p)),
                unit_to_internal("energy", "electronvolt", v.T),
                "",
            )
            for e, f, v, p in zip(energies, forces, virials, positions)
        ]
        return results if batched else results[0]

    def compute_structure(self, cell, pos):
        return self.compute(cell, pos)

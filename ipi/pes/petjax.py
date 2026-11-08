"""Interface with pet-jax (https://github.com/lab-cosmo/pet-jax), a JAX
implementation of the PET models, such as PET-MAD"""

import itertools

import numpy as np

from ipi.pes._ase import ASEDriver
from ipi.utils.array_backend import backend
from ipi.utils.messages import info, verbosity
from ipi.utils.units import unit_to_internal, unit_to_user

__DRIVER_NAME__ = "petjax"
__DRIVER_CLASS__ = "PETJAX_driver"

jax = None
jnp = None


class PETJAX_driver(ASEDriver):
    """Driver for PET models (e.g. PET-MAD) evaluated with pet-jax, from a
    checkpoint converted with `petjax-convert`.

    With `neighbors="vesin"` (the default), the structures are evaluated by the
    ASE calculator of pet-jax, which builds the neighbor lists with vesin on
    the host, keeps them while the atoms move less than `skin`, and runs the
    model with JAX. This works with any array backend of i-PI and any system
    size, but the steps are not compiled with `--jit`.

    With `neighbors="jax"`, the neighbor lists are computed with JAX, by
    checking all the pairs of atoms within the cell images that the cutoff
    reaches, and the model takes and returns the arrays of the backend (numpy
    or jax). The evaluation is a function of the positions and the cell only,
    so with `--backend jax --jit` it is compiled together with the rest of the
    step, and the structures of a batch (e.g. the beads) are evaluated in one
    call. The number of pairs and of selected neighbors per atom are fixed at
    the start, from the template, with the margins `pair_margin` and
    `extra_neighbors`: the simulation stops with an error if they are
    exceeded. The cost grows with the square of the number of atoms, so this
    is meant for up to a few thousand atoms.

    Command-line:
        i-pi-py_driver -m petjax -o template=init.xyz,checkpoint=pet-mad-xs [...]

    Parameters:
        :param template: string, filename of an ASE-readable structure file,
            which gives the atomic numbers, the periodicity and, for
            `neighbors="jax"`, the sizes of the neighbor lists
        :param checkpoint: string, a pet-jax checkpoint directory
        :param dtype: string, "float32" (default) or "float64"
        :param non_conservative: bool, whether to use the forces and stress
            predicted directly by the model (faster, but not conservative)
        :param neighbors: string, "vesin" (default) or "jax"
        :param skin: float, Verlet skin of the vesin neighbor lists, in angstrom
        :param pair_margin: float, ratio between the number of pairs in the
            neighbor lists computed with JAX and that of the template
        :param extra_neighbors: int, neighbors per atom selected by the model
            beyond those of the template, for `neighbors="jax"`
    """

    def __init__(
        self,
        template,
        checkpoint,
        dtype="float32",
        non_conservative=False,
        neighbors="vesin",
        skin=0.5,
        pair_margin=1.25,
        extra_neighbors=8,
        *args,
        **kwargs,
    ):
        if neighbors not in ("vesin", "jax"):
            raise ValueError(f"neighbors must be 'vesin' or 'jax', not '{neighbors}'")
        self.checkpoint = checkpoint
        self.dtype = dtype
        self.non_conservative = non_conservative
        self.neighbors = neighbors
        self.skin = skin
        self.pair_margin = pair_margin
        self.extra_neighbors = extra_neighbors
        super().__init__(template, *args, **kwargs)

    def check_parameters(self):
        """Loads the template and the model."""

        global jax, jnp

        super().check_parameters()
        try:
            import jax
            import jax.numpy as jnp
            import petjax
        except ImportError as err:
            raise ImportError(
                "could not import pet-jax, which the petjax driver needs: install it "
                "with `pip install pet-jax` (see https://github.com/lab-cosmo/pet-jax)"
            ) from err

        # the stress is only computed for periodic structures
        self.pbc = np.asarray(self.template_ase.pbc, dtype=bool)
        if not self.pbc.any() and "stress" in self.capabilities:
            self.capabilities.remove("stress")
        params, metadata = petjax.load_checkpoint(self.checkpoint)
        model = petjax.UPET(
            **metadata["config"],
            direct_forces=self.non_conservative,
            direct_stress=self.non_conservative and self.pbc.any(),
        )
        if self.neighbors == "vesin":
            # one calculator, with its own neighbor list, for each structure
            # of a batch (e.g. each bead), created when first needed
            self._calculators = []
            self._new_calculator = lambda: petjax.UPETCalculator(
                model,
                params,
                metadata,
                default_dtype=self.dtype,
                stress=bool(self.pbc.any()),
                skin=self.skin,
            )
        else:
            self._setup_jax(petjax, model, params, metadata)

    def _setup_jax(self, petjax, model, params, metadata):
        """Sets up the evaluation with neighbor lists computed with JAX, with
        sizes taken from the template."""

        from petjax.utils import cast_floats

        if backend() == "torch":
            raise ValueError(
                "neighbors='jax' needs the numpy or the jax backend of i-PI: use "
                "neighbors='vesin' with the torch backend"
            )
        self.backend_arrays = True
        self.traceable = True

        # the same global settings as the ASE calculator of pet-jax
        if self.dtype == "float64":
            jax.config.update("jax_enable_x64", True)
        jax.config.update("jax_default_matmul_precision", "high")
        self._dtype = jnp.float64 if self.dtype == "float64" else jnp.float32

        self._params = cast_floats(params, self._dtype)
        self._predict = petjax.get_predict_fn(model, stress=bool(self.pbc.any()))
        self._cutoff = float(model.cutoff)

        atoms = self.template_ase
        natoms = len(atoms)
        numbers = atoms.get_atomic_numbers()
        try:
            shift = sum(float(metadata["shifts"][int(z)]) for z in numbers)
        except KeyError as err:
            raise ValueError(
                f"The model was not trained on the atomic number {err.args[0]}"
            ) from err
        self._shift = unit_to_internal("energy", "electronvolt", shift)

        # the neighbor lists of the template give the sizes of those
        # computed in the steps (the last pair is a padding sentinel)
        structure = petjax.to_structure(
            atoms, self._cutoff, skin=0.0, bucket_strategy=None
        )
        npairs = int(structure["pair_mask"].sum())
        selection = petjax.select_edges(
            structure,
            model.num_neighbors_adaptive,
            self._cutoff,
            model.cutoff_width_adaptive,
            method=model.adaptive_cutoff_method,
        )
        k_sel = max(int(selection.counts.max()), 1) + self.extra_neighbors
        # an even number of tokens (the neighbors and the central atom) keeps
        # the attention on XLA's fused path
        self._k_sel = k_sel + k_sel % 2
        self._capacity = int(np.ceil(npairs * self.pair_margin)) + 1

        # the cell images that the cutoff reaches from atoms inside the cell
        cell = np.array(atoms.get_cell()[:], dtype=float)
        if not self.pbc.any():
            cell = np.eye(3)
        nimages = [
            int(np.ceil(self._cutoff / s)) if p else 0
            for s, p in zip(_plane_spacings(cell), self.pbc)
        ]
        shifts = np.array(
            list(itertools.product(*[range(-n, n + 1) for n in nimages])), dtype=int
        )
        index = {tuple(s): k for k, s in enumerate(shifts)}
        self._nimages = np.array(nimages)
        self._shifts = shifts
        self._negative = np.array([index[tuple(-s)] for s in shifts])
        self._zero = index[(0, 0, 0)]
        self._natoms = natoms
        self._numbers = np.concatenate([numbers, [0]])
        self._evaluate = jax.jit(self._evaluate_batch)
        info(
            f" @petjax: neighbor lists computed with JAX: {len(shifts)} cell images, "
            f"up to {self._capacity - 1} pairs ({npairs} in the template) and "
            f"{self._k_sel} selected neighbors per atom",
            verbosity.low,
        )

    def compute(self, cell, pos):
        """Evaluates one structure, or lists of them as one batch."""

        if self.neighbors == "vesin":
            if isinstance(cell, list):
                return [self._compute_with(k, *s) for k, s in enumerate(zip(cell, pos))]
            return self._compute_with(0, cell, pos)
        batched = isinstance(cell, list)
        cells = jnp.stack(cell) if batched else jnp.asarray(cell)[None]
        positions = jnp.stack(pos) if batched else jnp.asarray(pos)[None]
        energies, forces, virials = self._evaluate(cells, positions)
        if backend() == "numpy":
            energies, forces, virials = (
                np.asarray(x, dtype=float) for x in (energies, forces, virials)
            )
        # the composition shifts are large and constant: added in the dtype of
        # the backend (float64), not that of the model
        energies = energies + self._shift
        results = [
            (energies[b], forces[b], virials[b], "") for b in range(len(energies))
        ]
        return results if batched else results[0]

    def _compute_with(self, index, cell, pos):
        """Evaluates a structure with the calculator of the position `index`
        in a batch, whose neighbor list follows that structure."""

        while len(self._calculators) <= index:
            self._calculators.append(self._new_calculator())
        self.ase_calculator = self._calculators[index]
        return self.compute_structure(cell, pos)

    def _evaluate_batch(self, cells, positions):
        """The energies (without the composition shifts), forces and virials
        of a batch of structures, from their cells (i-PI's cell vectors as
        columns) and positions, in atomic units."""

        dtype = positions.dtype
        length = unit_to_user("length", "angstrom", 1.0)
        cell = (jnp.swapaxes(cells, 1, 2) * length).astype(self._dtype)
        pos = (positions.reshape(len(positions), -1, 3) * length).astype(self._dtype)
        energy, forces, virial, overflow = jax.vmap(self._evaluate_one)(cell, pos)
        jax.debug.callback(self._check_overflow, *[x.any() for x in overflow])

        energy = unit_to_internal("energy", "electronvolt", energy.astype(dtype))
        forces = unit_to_internal("force", "ev/ang", forces.astype(dtype))
        virial = unit_to_internal("energy", "electronvolt", virial.astype(dtype))
        virial = (virial + jnp.swapaxes(virial, 1, 2)) * 0.5
        return energy, forces.reshape(len(forces), -1), virial

    def _evaluate_one(self, cell, positions):
        """The energy, forces and virial of a structure in eV and angstrom,
        from its cell (vectors as rows) and positions, and whether the sizes
        of the neighbor lists are exceeded."""

        natoms, nimages = self._natoms, len(self._shifts)
        periodic = jnp.asarray(self.pbc)
        if not self.pbc.any():
            cell = jnp.eye(3, dtype=cell.dtype)

        # the pairs are found between the positions wrapped into the cell
        wrap = jnp.where(periodic, jnp.floor(positions @ jnp.linalg.inv(cell)), 0)
        wrapped = positions - wrap @ cell
        images = jnp.asarray(self._shifts, dtype=cell.dtype) @ cell
        d = wrapped[None, :, None, :] - wrapped[:, None, None, :] + images
        itself = jnp.eye(natoms, dtype=bool)[:, :, None] & (
            jnp.arange(nimages) == self._zero
        )
        within = (jnp.sum(d * d, axis=-1) < self._cutoff**2) & ~itself
        within = jax.lax.stop_gradient(within).reshape(-1)

        # the pairs, sorted by center (i, j, image), and their reverse pairs
        npairs = jnp.sum(within)
        capacity = self._capacity
        (index,) = jnp.nonzero(within, size=capacity, fill_value=0)
        i, j, k = (
            index // (natoms * nimages),
            (index // nimages) % natoms,
            index % nimages,
        )
        rank = jnp.cumsum(within) - 1
        reverse = rank[(j * natoms + i) * nimages + jnp.asarray(self._negative)[k]]
        valid = jnp.arange(capacity) < npairs
        # the image of j relative to i for the unwrapped positions
        shifts = jnp.asarray(self._shifts)[k] + (wrap[i] - wrap[j]).astype(int)

        dummy = natoms
        structure = {
            "positions": jnp.concatenate(
                [positions, jnp.zeros((1, 3), positions.dtype)]
            ),
            "cell": cell,
            "atomic_numbers": jnp.asarray(self._numbers),
            "atom_mask": jnp.arange(natoms + 1) < natoms,
            "centers": jnp.where(valid, i, dummy),
            "others": jnp.where(valid, j, dummy),
            "cell_shifts": jnp.where(valid[:, None], shifts, 0).astype(cell.dtype),
            "reverse": jnp.where(valid, reverse, jnp.arange(capacity)),
            "pair_mask": valid,
            "k_sel_sizer": jnp.zeros(self._k_sel, dtype=bool),
        }
        result = self._predict(self._params, structure)

        forces = result["forces"][:natoms]
        if self.non_conservative:
            # the net force of the direct forces is removed, as pet-jax does
            forces = forces - forces.mean(axis=0)
        if self.pbc.any():
            virial = -result["stress"]
        else:
            virial = jnp.zeros((3, 3), dtype=cell.dtype)
        # the cell images must cover the cutoff after the cell changes
        spacings = jnp.abs(jnp.linalg.det(cell)) / jnp.linalg.norm(
            jnp.cross(jnp.roll(cell, -1, axis=0), jnp.roll(cell, -2, axis=0)), axis=1
        )
        too_small = jnp.any(periodic & (self._cutoff > self._nimages * spacings))
        overflow = (npairs > capacity - 1, result["overflow"], too_small)
        return result["energy"], forces, virial, overflow

    def _check_overflow(self, pairs, neighbors, images):
        """Stops the simulation if the neighbor lists computed with JAX
        exceed the sizes fixed at the start."""

        if pairs:
            raise RuntimeError(
                "petjax: there are more pairs within the cutoff than the neighbor "
                "lists can hold: increase pair_margin"
            )
        if neighbors:
            raise RuntimeError(
                "petjax: the model selects more neighbors than the neighbor lists "
                "can hold: increase extra_neighbors"
            )
        if images:
            raise RuntimeError(
                "petjax: the cell has shrunk below the cutoff of the neighbor lists "
                "set from the template: use a template with a smaller cell"
            )


def _plane_spacings(cell):
    """The distances between the lattice planes of a cell (vectors as rows)."""

    volume = abs(np.linalg.det(cell))
    return volume / np.linalg.norm(np.cross(cell[[1, 2, 0]], cell[[2, 0, 1]]), axis=1)

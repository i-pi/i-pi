# PET-MAD through pet-jax

Runs short path-integral simulations of 32 water molecules with the universal
[PET-MAD](https://github.com/lab-cosmo/upet) machine-learning potential
(version 1.6, size "xs"), evaluated with
[pet-jax](https://github.com/lab-cosmo/pet-jax), a JAX implementation of the
PET models, through the `petjax` driver of i-PI.

## Installation

Install pet-jax in the same environment as i-PI:

```bash
python -m pip install pet-jax
```

or, for the latest version, `python -m pip install "git+https://github.com/lab-cosmo/pet-jax"`.
Install a CUDA-enabled version of JAX (e.g. `python -m pip install "jax[cuda12]"`)
to run on a GPU.

`input-jit.xml` needs a version of pet-jax that computes the neighbor lists
with [vesin-jax](https://github.com/Luthaf/vesin) (0.6.2 or later), currently
that of the `vesin-jax` branch. vesin-jax is not on PyPI yet, and is built from
the vesin repository when pet-jax is installed, which needs `cmake` and a C++
compiler. It supports jax 0.10 and 0.11: pip downgrades a newer jax. The
versions that are known to work together are in `requirements/petjax.txt`, at
the root of the i-PI repository:

```bash
python -m pip install -r requirements/petjax.txt
```

| package | needed for | version |
|---|---|---|
| pet-jax | the `petjax` driver | `vesin-jax` branch for `neighbors:jax` |
| jax, jaxlib | pet-jax, `--backend jax` | 0.10 or 0.11 with vesin-jax |
| flax, ase, vesin, marathon-train | pet-jax | installed with pet-jax |
| vesin-jax | `neighbors:jax` | 0.6.2, from the vesin repository |
| torch, metatomic-torch, metatrain | converting the checkpoints (`pet-jax[convert]`) | |

The model is not included in this folder. pet-jax reads its own checkpoints,
converted from those of metatrain: the conversion needs some further packages
(torch, metatomic and metatrain), which are not needed to run the model.
`getmodel.sh` downloads PET-MAD-xs v1.6.0 from Hugging Face and converts it
into the `pet-mad-xs/` folder:

```bash
python -m pip install "pet-jax[convert]"
bash getmodel.sh
```

## Running the examples

`input.xml` uses the default `neighbors:vesin` option of the driver: the
neighbor lists are computed by vesin, and kept while the atoms move less than
a Verlet skin, as in the ASE calculator of pet-jax. This works with any array
backend of i-PI, and for systems of any size.

```bash
i-pi input.xml
```

`input-jit.xml` uses `neighbors:jax`: the neighbor lists are computed by
vesin-jax inside JAX, with the same Verlet skin, so that evaluating the model
only takes the positions and the cell. With the jax backend and `--jit`, i-PI
then compiles the whole MD step, including the evaluation of the model on the
four beads (as one batch), into a single call:

```bash
i-pi --backend jax --device gpu --jit input-jit.xml
```

(use `--device cpu` without a GPU). The log reports `step compiled` after the
first few steps. The sizes of the neighbor lists (the number of pairs of each
atom, and of the neighbors that the model selects) are fixed at the start, from
the template structure, with some margin (the `pair_margin` and
`extra_neighbors` parameters): i-PI stops with an error if a simulation
exceeds them. vesin-jax needs the 64-bit mode of JAX, which the driver enables;
the model still runs in single precision, unless `dtype:float64` is set.

The other parameters of the driver are listed in its documentation
(`ipi/pes/petjax.py`): `non_conservative:True` uses the forces and stress
predicted directly by the model, which are faster but not conservative, and
`dtype:float64` runs the model in double precision.

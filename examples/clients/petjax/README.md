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

`input-jit.xml` uses `neighbors:jax`: the neighbor lists are computed with
JAX, by checking all the pairs of atoms (with the periodic images that the
cutoff reaches), so that evaluating the model only takes the positions and
the cell. With the jax backend and `--jit`, i-PI then compiles the whole MD
step, including the evaluation of the model on the four beads (as one batch),
into a single call:

```bash
i-pi --backend jax --device gpu --jit input-jit.xml
```

(use `--device cpu` without a GPU). The log reports `step compiled` after the
first few steps. The sizes of the neighbor lists are fixed at the start, from
the template structure, with some margin (the `pair_margin` and
`extra_neighbors` parameters): i-PI stops with an error if a simulation
exceeds them. As the cost of finding the pairs grows with the square of the
number of atoms, this option is meant for systems of up to a few thousand
atoms.

The other parameters of the driver are listed in its documentation
(`ipi/pes/petjax.py`): `non_conservative:True` uses the forces and stress
predicted directly by the model, which are faster but not conservative, and
`dtype:float64` runs the model in double precision.

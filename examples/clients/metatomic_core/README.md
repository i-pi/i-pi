# metatomic-core example

NVE dynamics of 108 argon atoms with the `metatomic_core` driver, which runs any
model implementing the [metatomic](https://docs.metatensor.org/metatomic) C API
through the backend-independent `metatomic` Python package: TorchScript models
(through the metatomic-torch plugin), JAX or NumPy models written in Python, and
native plugins.

## Installation

```bash
pip install metatomic vesin ase
# for TorchScript (.pt) models, also: pip install metatomic-torch
```

## Running the example

```bash
i-pi input.xml & sleep 2

# a Lennard-Jones model shipped as a native plugin (metatomic-lj-plugin.so)
i-pi-py_driver -a metatomic-core -u -m metatomic_core -o template=argon.xyz,model=metatomic-lj-model,plugins=lj,option_sigma=3.4,option_epsilon=0.0104,option_cutoff=7.0,option_atomic_type=18,option_length_unit=angstrom,option_energy_unit=eV

# or the same potential as a metatomic.NumpyModel defined in a Python file
i-pi-py_driver -a metatomic-core -u -m metatomic_core -o template=argon.xyz,model=lj_model.py:model
```

Other kinds of models are selected with the `model` option:

- `model=file.py:Name`, where `Name` is a `metatomic.BaseModel` instance, a
  subclass (instantiated without arguments) or a function returning a model;
- `model=model.pt`, a TorchScript model, loaded with the metatomic-torch plugin;
- `model=<name>`, anything else `metatomic.load_model` accepts, after loading
  the plugins given in `plugins=a.so;b.so` (bare names are looked up in
  metatomic's `libexec/metatomic` directory). `option_<key>=<value>` options are
  given to the plugin, and `plugin_name` selects which plugin should load the
  model.

The other options are `energy_variant` (pick `energy/<variant>` for multi-head
models), `check_consistency`, `enable_x64` (enable 64-bit floats in JAX before
building the model, default true) and `skin` (neighbor list skin, in Angstrom).
The model is evaluated on CPU, and returns the energy, the forces and (when the
model provides strain gradients) the virial. The periodic boundary conditions
are taken from the template file.

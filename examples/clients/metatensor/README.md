# Metatensor example

Runs an example of a metatensor model, implementing a basic Lennard-Jones model
with the parameters for Nickel.

This driver can be used with any model following the [metatensor atomistic
models](https://lab-cosmo.github.io/metatensor/latest/atomistic/index.html)
interface. This can be used for classical force-fields, but it more intended for
machine learning potentials.

## Installation

The code is compatible with metatensor-torch v0.5, which you can install with

```bash
# install all metatensor packages simultaneously
pip install "metatensor[torch]"

# install packages individually, with explicit control over the installed versions
pip install "metatensor-torch ==0.5.*" "metatensor-operations ==0.2.*"
```

## Running the example

```bash
i-pi input.xml & sleep 1 
i-pi-py_driver -a metatensor -u -m metatensor -o nickel.xyz,nickel-lj.pt

# with all the optional parameters:
i-pi-py_driver -a metatensor -u -m metatensor -o nickel.xyz,nickel-lj.pt,device=cpu,extensions=some-extensions-dir/,check_consistency=True
```

The options (after `-o`) are as follow:

- the path to a template file, we will use it to get the types for all atoms in
  the system
- the path to the model file
- `device` controls which torch device to use to run the model
- `extensions` is the path to a directory containing TorchScript extensions. If
  the model requires such extensions, we will try to load them from this
  directory first
- `check_consistency` controls whether we should run some extra internal
  consistency checks about data given to the model and data returned by the
  model.

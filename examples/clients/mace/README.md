# MACE Example

Runs an example of the pretrained [MACE model](https://github.com/ACEsuit/mace-mp/releases/download/mace_mp_0c/mace-small-density-agnesi-stress.model) for simple molecular dynamics of liquid water at 300 K.

## Installation

To be able to run MACE, check out the installation instructions at [https://github.com/ACEsuit/mace/tree/main?tab=readme-ov-file#installation](https://github.com/ACEsuit/mace/tree/main?tab=readme-ov-file#installation).

## Running the Example

To download the pretrained [MACE model](https://github.com/ACEsuit/mace-mp/releases/download/mace_mp_0c/mace-small-density-agnesi-stress.model):

```bash
bash getmodel.sh
```

Then, run the following to start the i-PI simulation:

```bash
i-pi input.xml & sleep 10
i-pi-py_driver -a driver -u -m mace -o init.xyz,mace.model
```

### Options (after `-o`)

- Path to a template file, used to retrieve the types for all atoms in the system.
- Path to the model file.


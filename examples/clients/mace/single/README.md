# MACE Example

Runs a ten-step molecular dynamics simulation of a water molecule at 300 K
using the pretrained
[MACE-MP-0b2 small model](https://huggingface.co/mace-foundations/mace-mp-0/blob/main/mace-mp-0b2-small.model).
The model is MIT-licensed and approximately 68 MB.

## Installation

Install MACE in the same Python environment as i-PI:

```bash
python -m pip install mace-torch
```

See the [MACE installation instructions](https://github.com/ACEsuit/mace#installation)
for other installation methods and accelerator-specific options.

## Running the Example

Download the pretrained model as `mace.model`:

```bash
bash getmodel.sh
```

Then start the simulation from this directory:

```bash
i-pi input.xml
```

The `<ffdirect>` force field loads MACE directly inside the i-PI process, so no
separate `i-pi-py_driver` process is required:

```xml
<ffdirect name='driver' pbc='False'>
  <pes>mace</pes>
  <parameters>{template:init.xyz,model:mace.model,device:cpu}</parameters>
</ffdirect>
```

The parameters select:

- `template:init.xyz`, which supplies the atomic species expected by MACE;
- `model:mace.model`, the model downloaded by `getmodel.sh`; and
- `device:cpu`, which runs the calculation on the CPU.

For an example that evaluates multiple ring-polymer beads in one MACE batch,
see [`../mace-batched`](../mace-batched/README.md).

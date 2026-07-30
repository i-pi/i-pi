# BaTiO3 with a local MACE-POLAR-1 checkpoint

This example runs MACE-POLAR-1-M through i-PI's in-process `mace` driver. The
checkpoint is stored in this directory and passed to the driver like a normal
MACE `.model` file; no socket client or custom PES is needed.

MACE-POLAR-1 is a molecular foundation model trained on OMol25. This periodic
BaTiO3 cell is therefore a wiring/smoke-test example, not a validated BaTiO3
potential. In particular, the total dipole reported by MACE-POLAR is not
well-defined under periodic boundary conditions and is deliberately not
written by `input.xml`.

## Install

Use one Python environment for MACE and i-PI. MACE-POLAR support requires a
recent MACE version and its long-range graph dependency:

```bash
python -m pip install -U 'mace-torch>=0.3.16'
python -m pip install 'git+https://github.com/WillBaldwin0/graph_electrostatics.git'
python -m pip install -e ../../../..
```

The MACE-POLAR checkpoints use the Academic Software License (ASL). Review and
accept its terms before downloading a checkpoint.

## Download and run

From this directory:

```bash
./download_mace_polar.sh
i-pi input.xml
```

The download script writes `MACE-POLAR-1-M.model` here. `input.xml` then loads
it with:

```xml
<ffdirect name='shell' pbc='True'>
  <pes>mace</pes>
  <parameters>{
    template:start.extxyz,
    model:MACE-POLAR-1-M.model,
    device:cpu,
    mace_kwargs:mace_kwargs.json
  }</parameters>
</ffdirect>
```

`mace_kwargs.json` sets `model_type` to `PolarMACE`, which is required even
though the checkpoint is loaded from an ordinary local path. It uses float32
for MD and avoids copying unused charge-density tensors from the accelerator.
Use `device:cuda` for a CUDA GPU. For higher-accuracy static checks, change
`default_dtype` to `float64`.

When `instructions.compute_BEC` is enabled, i-PI uses the strain-corrected
(`proper`) dipole by default to calculate Born effective charges and the
piezoelectric tensor. To use the raw dipole returned by the model instead, set
`use_proper_dipole` to `false`, either in `mace_kwargs.json` or directly in the
force-field parameters:

```xml
<parameters>{...,use_proper_dipole:false}</parameters>
```

`start.extxyz` supplies the global electronic state expected by MACE-POLAR:
neutral charge, singlet spin multiplicity, and zero external field. Change
these metadata when simulating a different charge, spin state, or static
field.

## Selecting MACE outputs

MACE-POLAR can return many optional tensors in addition to the energy, forces,
stress, and dipole. The exact set can vary between MACE versions and
checkpoints. Transferring every tensor from a GPU is unnecessary when i-PI does
not use it, so `mace_kwargs.json` ignores the optional outputs encountered by
this example.

At the first force evaluation, i-PI prints all outputs produced by the model.
Registered outputs are grouped into per-atom and per-structure properties;
unregistered outputs are printed with their raw tensor shapes. If several
unsupported outputs remain, i-PI reports all of them in one error rather than
stopping at the first one. For example:

```text
Unknown model properties for a batch with natoms=[5]:
- 'atomic_embedding': raw shape (5, 16); inferred per-atom shape ["natoms", 16]
- 'graph_tensor': raw shape (1, 3, 3); inferred per-structure shape [3, 3]
- 'spin_density': raw shape (5, 2, 4); inferred per-atom shape ["natoms", 2, 4]
```

The same error provides two consolidated JSON snippets. If the properties are
not needed, merge every reported name into the existing ignore list:

```json
{
    "instructions": {
        "ignore": [
            "atomic_embedding",
            "graph_tensor",
            "spin_density"
        ]
    }
}
```

To retain an output in the i-PI extras instead, register its shape at the top
level of `mace_kwargs.json` and do not put that property in `ignore`:

```json
{
    "model_type": "PolarMACE",
    "default_dtype": "float32",
    "ase_like_properties": {
        "atomic_embedding": ["natoms", 16],
        "graph_tensor": [3, 3],
        "spin_density": ["natoms", 2, 4]
    },
    "instructions": {
        "ignore": ["other_unused_property"]
    }
}
```

`"natoms"` marks a per-atom leading dimension. A shape without `"natoms"` is
per-structure; i-PI removes the leading batch dimension shown in the raw
tensor. These shapes are inferred from dimensions alone, so verify that the
suggestion matches the physical meaning of the property before retaining it.

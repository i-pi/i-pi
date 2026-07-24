# Batched MACE Example

Runs an example of the pretrained [MACE-MP-0b2 small model](https://huggingface.co/mace-foundations/mace-mp-0/blob/main/mace-mp-0b2-small.model) for path-integral molecular dynamics of a water molecule at 300 K. The four ring-polymer beads are collected into one in-process MACE batch. The model is MIT-licensed and approximately 68 MB.

## Installation

To be able to run MACE, check out the installation instructions at [https://github.com/ACEsuit/mace/tree/main?tab=readme-ov-file#installation](https://github.com/ACEsuit/mace/tree/main?tab=readme-ov-file#installation).

## Running the Example

To download the pretrained MACE-MP-0b2 small model as `mace.model`:

```bash
bash getmodel.sh
```

Then, start the i-PI simulation from this directory:

```bash
i-pi input.xml
```

The `<ffdirect>` force field loads MACE directly inside the i-PI process. No separate `i-pi-py_driver` process is needed.

## How batching is configured

The relevant part of [`input.xml`](input.xml) is:

```xml
<ffdirect name='driver' pbc='False'>
  <pes>mace</pes>
  <parameters>{template:init.xyz,model:mace.model,device:cpu,mace_kwargs:mace_kwargs.json}</parameters>
  <batch_size>4</batch_size>
</ffdirect>
```

These settings mean:

- `template:init.xyz` supplies the atomic species expected by MACE. This three-atom extended XYZ configuration is copied from `examples/init_files/h2o.extxyz`; its positions and periodic cell are expressed in angstrom.
- `model:mace.model` selects the model downloaded by `getmodel.sh`. The downloaded model is not part of the repository.
- `device:cpu` runs the example on the CPU. Change it to `device:cuda` to use a CUDA GPU.
- `mace_kwargs:mace_kwargs.json` loads the additional MACE settings described below.
- `<batch_size>4</batch_size>` tells i-PI to collect four force requests before calling MACE.

The system also uses `nbeads='4'`, so the four ring-polymer beads fill one complete batch at every force evaluation. `BatchedMACE` evaluates all four structures together unless a smaller internal batch size is explicitly configured.

## Avoiding unnecessary GPU-to-CPU transfers

MACE can produce optional properties that i-PI does not use. Copying those tensors from a GPU to the CPU costs time. The driver therefore prints a summary once at startup showing:

- every property available after the MACE evaluation; and
- the subset copied to the CPU for i-PI.

The committed [`mace_kwargs.json`](mace_kwargs.json) contains:

```json
{
  "instructions": {
    "ignore": ["displacement"]
  }
}
```

Here, MACE still creates `displacement` because it is needed internally to calculate stress, but the driver does not copy it to the CPU as a returned result.

To suppress more optional results, add their names to the `ignore` list after checking the startup summary. Never ignore `energy`, `forces`, or `stress`, because i-PI requires them.

## Recorded output

[`example.txt`](example.txt) contains the terminal output from a successful ten-step CPU run. In particular:

- `Launching batch evaluation, 4 / 4` confirms that all four beads were collected into one call.
- `MACE output summary (printed once)` confirms that the output report is only emitted once.
- The CPU-output lists show that `forces`, `energy`, and `stress` were returned, while `displacement` was not copied.

The reported timings and optional acceleration warnings depend on the machine and installed packages.

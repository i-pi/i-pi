# MACE example

This example runs a ten-step molecular dynamics simulation of a water molecule
using the pretrained MACE-MP-0b2 small model.

## Run

Install MACE in the same Python environment as i-PI:

```bash
python -m pip install -r ../../../requirements/mace.txt
```

Then download the model and run the basic in-process example:

```bash
bash getmodel.sh
i-pi input.xml
```

## Execution-path comparison

The [`comparison`](comparison) folder compares three equivalent NVE
trajectories: ASE `MACECalculator` over `ffsocket`, the bundled
`i-pi-py_driver` over `ffsocket`, and native `ffdirect`. It also checks the
trajectories and plots end-to-end timings.

```bash
cd comparison
./run.sh
```

The comparison runner downloads `../mace.model` automatically if necessary.

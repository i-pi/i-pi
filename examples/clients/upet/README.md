# PET-MAD through the metatomic driver

Runs short simulations of 32 water molecules with the universal
[PET-MAD](https://github.com/lab-cosmo/upet) machine-learning potential
(version 1.6, size "xs"), using the metatomic driver
(see also `examples/clients/metatomic`).

## Installation

Install the compatible upet/metatomic dependency set in the same environment as i-PI:

```bash
python -m pip install -r ../../../requirements/upet.txt
```

Then export the model to a TorchScript file, `pet-mad-xs-v1.6.0.pt`:

```bash
bash getmodel.sh
```

## Running the examples

`input.xml` runs path-integral MD with four beads, using a `<ffdirect>` force
field that evaluates all beads as a single batch. The batch is split between two
model instances running in parallel threads (`device:cpu;cpu`; on GPUs you would
use e.g. `device:cuda:0;cuda:1`). The per-atom energy uncertainty is written to
`simulation.uq_*`, and a warning is printed if it exceeds `uncertainty_threshold`.

```bash
i-pi input.xml
```

`input-socket.xml` runs the same simulation over a UNIX socket. All beads are
sent to `i-pi-py_driver` in one batched request. The driver also computes an
ensemble of energies, which is written to `socket.committee_pot_*`.

```bash
i-pi input-socket.xml & sleep 5
i-pi-py_driver -a upet -u -m metatomic \
    -o template=init.xyz,model=pet-mad-xs-v1.6.0.pt,energy_ensemble=True
```

`input-rotations.xml` uses the fast, non-conservative forces predicted directly by
PET-MAD. These are neither conservative nor exactly rotationally equivariant, so
they are wrapped in `<ffrotations>`. At every step it evaluates a random rotation of
the structure and its inversion as a single batch, and averages the two predictions.
A strong Langevin thermostat absorbs the resulting noise.

```bash
i-pi input-rotations.xml
```

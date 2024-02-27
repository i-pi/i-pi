metatensor - i-PI example
========================

Runs an example of a metatensor model, implementing a basic Lennard-Jones model
with the parameters for Nickel.

This driver can be used with any model following the [metatensor atomistic
models](https://lab-cosmo.github.io/metatensor/latest/atomistic/index.html)
interface. This can be used for classical force-fields, but it more intended for
machine learning potentials.


```bash
i-pi input.xml; sleep 1 &
i-pi-py_driver -a metatensor -u -m metatensor -o nickel.xyz,nickel-lj.pt
```

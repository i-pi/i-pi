metatensor - i-PI example
========================

Runs a _very_ contrieved example of a metatensor model. 
The model is just a hard-coded Einstein crystal model for a
chunk of diamond. No periodic boundary conditions, no ability
to work with a different number of atoms or just a differently
oriented sample. Really, this is just to show how to run the
driver. Given that metatensor model takes a torchsript file
for a model, the very same machinery can be used for actual
ML potentials, or torchscript-implemented models. 

```bash
i-pi input.xml; sleep 1 &
i-pi-py_driver -a metatensor -u -m metatensor -o initial.xyz,harmonic-model.pt
```

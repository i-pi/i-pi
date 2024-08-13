Metadynamics calculations using PLUMED
======================================

There are three examples of different variations-on-a-theme of the use of `<ffplumed>`
to perform metadynamics calculations with i-PI

`pimd_metadynamics` : metadynamics using a centroid bias for a PIMD calculation
`wte_metadynamics`  : well-tempered ensemble. requires a workaround to communicate the energy as a CV
`remd_metadynamics` : replica exchange combined with metadynamics

These examples can be run using the Zundel PES implemented in the i-pi driver,

```bash
i-pi-driver -u -h zundel -m zundel
```

However, they require having PLUMED libraries and Python bindings in the appropriate system paths.

Note also that the example in `pimd_metadynamics` provides a demonstration of how to retrieve quantities
computed PLUMED-side (such as CV values) back into i-PI as `extras` that can be associated with 
individual structures and printed out alongside the other i-PI outputs. This functionality requires
a recent version of PLUMED, at least 2.10.

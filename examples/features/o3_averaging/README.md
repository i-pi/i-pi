O3 rotational averaging. 
========================

Examples of the averaging module to handle non-equivariant potentials.
The `<ffrotations>` socket allows to evaluate a potential multiple times,
using random rotations or a grid of rotations, to reduce the impact of the
lack of rotational equivariance on the quality of the dynamics. 
The examples use the `noo3_h2o` potential, that assumes the structure is 
a sequence of H2O molecules and adds a orientation-dependent potential to each
water molecule.


`noo3`: no correction is applied. The structure of the liquid will be badly
    disrupted.
`random`: the structure is randomly oriented at each evaluation. 
    This effectively introduces a noise in the forces, that must be contrasted
    by an aggressive thermostat. There is a large drift in the conserved quantity
    but the mean potential is OK and water molecules do not orient artificially
`grid`: the potential is averaged over a grid of rotations. 
    This reduces the orientational error, but comes at the cost of many additional
    energy evaluations per step.

Conservative ring-polymer mean-field examples
=============================================

`NVE-f/DW-NVE-bafab` and `NVE-f/DW-NVE-fbabf` demonstrate the bath-induced
mean-field correction in a one-dimensional double well. Both explicitly use
`bath_mode=none`: they contain no stochastic or dissipative bath propagation.
The two historical splitting names are retained for comparison.

Run `bash run.sh` from either example with i-PI and its Python driver installed.
The driver supplies the exact coupling `g(x)=x*SD(x)` and its derivative.
`friction_lambda.dat` supplies `[omega, Lambda(omega)]` in atomic units and
sets the correction strength. The retained `fric_SD.dat` contains J(omega);
the new table is its conversion Lambda=J/omega, with the zero-frequency value
obtained from the exponential form. Omit the driver's legacy eta0 argument: nonzero values are rejected. A zero
value is only a compatibility placeholder and does not disable friction.

The correction is `U_mf = 0.5 sum_n alpha_n F_n^2`, with alpha integrated
from Lambda. The current output property is `friction_energy_mf`.

For a complete RPMD-GLE bath, supply externally fitted Ap matrices for every
normal mode and use an appropriate non-Markovian splitting. See the electronic
friction documentation for the mode ordering, normalization, restart and
backend requirements. These examples alone are not RPMD-GLE TCF benchmarks.

Classical molecular dynamics. 
=============================

Examples of classical molecular dynamics simulations. 
The NPzT simulations are barostatted with an anisotropic cell barostat, and must be run with LAMMPS,
as the built-in driver only supports cubic cells for the TIP4P/f potential

`nvt_ensemble`: global thermostatted classical NVT simulations
`NPzT_ensemble_BZP`: classical NVT simulations with the BZP barostat with only the z-component cell allowed to fluctuate
`NPzT_ensemble_MTTK`: classical NVT simulations with the MTTK barostat with only the z-component cell allowed to fluctuate

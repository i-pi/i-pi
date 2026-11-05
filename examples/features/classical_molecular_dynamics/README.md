Classical molecular dynamics
============================

Examples of classical molecular dynamics simulations of water, at constant volume, pressure or stress.
The ice examples get the forces from LAMMPS (`in.lmp`), as the built-in q-TIP4P/F driver only supports cubic cells.

`nvt_ensemble_liquid_water`: liquid water with a global (SVR) thermostat

`NPT_ensemble_SCR`: the stochastic cell-rescaling barostat for a harmonic model, which runs without an external code

`NPT_ensemble_SCR_water`: liquid water at constant pressure with the stochastic cell-rescaling barostat

`NPT_ensemble_BZP_ice`: ice at constant pressure with the isotropic (BZP) barostat, which only changes the volume

`NPT_ensemble_MTTK_ice`: ice at constant pressure with the flexible (MTTK) barostat, which changes the lengths and
angles of the cell

`NST_ensemble_RGB_ice`: ice at constant stress with the anisotropic (RGB) barostat

`NPzT_ensemble_RGB_ice`: ice at constant stress with the anisotropic barostat, with only the z component of the
cell allowed to fluctuate

`NVsigmaT_ensemble_MTTK_ice`: ice in the (N, V, $\sigma_a = 0$, T) ensemble with the flexible barostat, where the
cell fluctuates at constant volume

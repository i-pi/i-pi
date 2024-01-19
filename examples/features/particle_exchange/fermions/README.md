Trapped cold non-interacting fermions
=================================================
An example of how to obtain fermionic statistics from a bosonic simulation, here, of 3 non-interacting fermions in an anisotropic harmonic trap.
The fermions are treated as bosons throughout the simulation. An additional property, 'fermionic_sign', is recorded. Then, in the analysis, the value of the estimator should be multiplied by the fermionic sign in each step, and the total result should be divided the average fermionic sign. See Hirshberg et al.'s doi:10.1063/5.0008720.

When the average sign is close to zero, the results become harder to converge. This is the reason for the choice of a higher temperature compared to other examples of trapped bosons. P=12 beads suffice for convergence in this temperature. The analytical result for the total energy is 0.00091881 Ha, higher than the result for distinguishable particles (0.00078146 Ha) or bosons (0.0006917 Ha) in the same setting.

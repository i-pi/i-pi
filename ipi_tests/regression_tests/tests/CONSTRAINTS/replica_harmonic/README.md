Reference data use the exact continuous-time harmonic solution, not a recorded
i-PI trajectory. The atom has mass 2, the physical spring constant is 0.2,
and the two-replica spring frequency is 2 k_B T / hbar at 10000 K.

With replica 0 fixed, replica 1 oscillates about
`2 m omega_P^2 q0 / (0.2 + 2 m omega_P^2)` with frequency
`sqrt(0.2/m + 2 omega_P^2)` and initial momentum 0.3.

With the centroid fixed, both coordinates oscillate about 0.75 with frequency
`sqrt(0.2/m + 4 omega_P^2)`. Initial momenta after centroid projection are
-0.15 and +0.15. The centroid position and momentum are constant.

The absolute tolerance of 1e-7 allows the second-order splitting error over
10 steps of 0.01 atomic time units; unit tests check fixed coordinates much
more tightly and also check longer trajectories, thermal statistics, and restart.

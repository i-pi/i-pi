# Fixed-replica and fixed-centroid PIMD

Run `i-pi input.xml`. This example uses the built-in harmonic forcefield,
so no external driver is needed.

`<fixbeads>[0]</fixbeads>` inside `<motion mode="dynamics">` fixes all
atomic coordinates of replica 0 at their initial values. Replica indices are
zero-based. Multiple replicas can be fixed, for example `[0, 4]`.
The other replicas continue to move, including their spring motion.

To hold each atomic centroid fixed instead, replace the `fixbeads` line with
`<fixcentroid>True</fixcentroid>`. This replaces the former `nvt-cc` dynamics
mode: use ordinary `nvt` dynamics with the motion option. The
`constrained_centroid_liquid_water` example also uses this interface.

Both options support NVE and NVT dynamics and either OBABO or BAOAB splitting
(with multiple time steps). They cannot be enabled together, and automatically
disable `fixcom`. They can be combined with `fixatoms` or `fixatoms_dof`;
shared constraints are counted once for thermostat degrees of freedom.
Constraints and current positions survive XML checkpoint/restart.
Within `MultiMotion`, all motions must use the same fixed-bead set;
different sets (including fixed beads in only some motions) are not implemented.

Fixed replicas require ordinary closed, distinguishable-particle Trotter
springs, RPMD dynamical masses, the exact spring propagator, and (for NVT)
a Cartesian Langevin thermostat. Fixed centroids retain the thermostat
options of the former `nvt-cc` mode and do not support bosonic exchange. Unsupported
combinations raise an error before integration.

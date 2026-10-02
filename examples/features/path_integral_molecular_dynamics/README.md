Path integral molecular dynamics. 
=================================================
`standard_pimd`: path integral molecular dynamics using the PILE thermostat and MTS for propagating the ring polymer modes.
`cayley_pimd`: path integral molecular dynamics using the PILE thermostat and Cayley integrator for propagating the ring polymer modes.
`constrained_centroid`: constant-temperature PIMD, with the centroid held fixed in the initial position. Useful to compute the centroid potential of mean force.
`piglet`: path integral molecular dynamics using the PIGLET thermostat
`pimd+mts`:  path integral molecular dynamics with MTS algorithm to integrate short and long ranged forces with different timesteps. Includes option to print out the slow and/or fast force components.
`rpc`:  path integral molecular dynamics with RPC algorithm to use different number of replicas for short and long ranged forces.
scpimd:  path integral molecular dynamics with the Suzuki-Chin splitting.
`eco_pimd`: path integral molecular dynamics with economised ("eco") ring-polymer frequencies, that accelerate convergence with the number of beads by fitting the spring frequencies to reproduce the radii of gyration of harmonic oscillators up to a maximum physical frequency.
`standard_constant_pressure`: path integral molecular dynamics using the PILE thermostat and MTS for propagating the ring polymer modes in the NPT ensemble.

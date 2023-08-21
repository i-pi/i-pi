Path integral molecular dynamics. 
=================================================
standard pimd: path integral molecular dynamics using the PILE thermostat and MTS for propagating the ring polymer modes.
cayley pimd: path integral molecular dynamics using the PILE thermostat and Cayley integrator for propagating the ring polymer modes.
piglet: path integral molecular dynamics using the PIGLET thermostat
pimd+mts:  path integral molecular dynamics with MTS algorithm to integrate short and long ranged forces with different timesteps.
rpc:  path integral molecular dynamics with RPC algorithm to use different number of replicas for short and long ranged forces.
scpimd:  path integral molecular dynamics with the Suzuki-Chin splitting.
replica_exchange: path integral molecular dynamics with replica-exchange acorss high-temperatures and pressures.
standard_constant_pressure: path integral molecular dynamics using the PILE thermostat and MTS for propagating the ring polymer modes in the NPT ensemble.
two_temperature_PILE: an example of how one could run PACMD with a two-temperature PILE protocol. Settings cannot prevent "temperature leakage", do not use for production.

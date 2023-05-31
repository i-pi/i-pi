Weighted baseline model
=======================

This runs a simulation in which one potential is meant as a "safe baseline"
and a committee of potentials is used as a "correction". The uncertainty in the
correction is estimated from the spread of the potentials, as in [doi.org/10.1021/acs.jctc.8b00959].

Here we take as a simple example that does not require specialized ML software
a simulation of water in which the baseline is taken to be the monomer "bound" term,
and the LJ and long-range parts is the correction. The difference between the two
members of the committee is just the kspace cutoff. 

The baseline error is taken to be rather large (it is!) so that it will basically
never switch off the correction and the simulation is stable. 

The committee model also computes useful quantities for post-processing. 
This example also showcases how to use ring polymer contraction on these quantities
by using the `force_extras` keyword.

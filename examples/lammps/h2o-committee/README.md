Thermodynamic uncertainty estimation
====================================

This runs a simulation in which several potentials are computed, and 
the uncertainty in the potential correction is estimated from the spread of the
members of the comittee, as in [doi.org/10.1021/acs.jctc.8b00959].

Here we take as a simple example that does not require specialized ML software
a simulation of water in which the members differ by the kspace cutoff and the
details of the LJ potential.
The simulation prints out files that are needed to compute how the error in the 
potential propagates on thermodynamic properties, specifically on the g(r).

This also showcases the "active learning" feature, in which structures
with a large error are printed out so they can be further analyzed.

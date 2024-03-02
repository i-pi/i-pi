Weighted baseline model
=======================

This runs a simulation in which one potential is meant as a "safe baseline"
and a committee of potentials is used as a "correction". The uncertainty in the
correction is estimated from the spread of the potentials, as in [Imbalzano 2021](http://doi.org/10.1021/acs.jctc.8b00959).

Here we take as a simple example that does not require specialized ML software
a simulation of water in which the members differ by the application of a small 
LJ-like perturbation.
The simulation prints out files that are needed to compute how the error in the 
potential propagates on thermodynamic properties, specifically on the g(r).

The committee model also computes useful quantities for post-processing. 
This example also showcases how to use ring polymer contraction on these quantities
by using the `force_extras` keyword.

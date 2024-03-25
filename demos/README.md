i-PI demonstrations
===================

These examples provide more sophisticated and realistic simulations setups, that 
demonstrate the usage of i-PI in a number of actual application scenarios. 
This README provide an overview of the content of the various folders.


p-H2: a tutorial
----------------

This is an introductory tutorial demonstrating classical and quantum simulations
of para-hydrogen close to a phase transition. The example is discussed as a tutorial
as part of the main documentations of i-PI. 


Alchemical isotope exchanges
----------------------------

This example demonstrates the use of the alchemical isotope exchange method, which simplifies
sampling equilibrium isotope fractionation between different molecular sites or different
thermodynamic phases of matter.


Kinetic Monte Carlo for Al6XXX
------------------------------

This is an experimental (and somewhat hacky) implementation of kinetic monte carlo, 
specifically designed for the study of vacancy diffusion in aluminum alloys of the 6000 series (Al alloyed with ~1% Si and Mg).
It demonstrates the use of i-PI with a ML potential, through LAMMPS patched with N2P2. Generalizing the implementation would require some coding, but most of the heavy lifting is already done. 

2D-IR-Raman
-----------

This example shows how one can use i-PI to compute 2D spectrum corresponding to a hybrid IR-Raman pulse sequence. The example contains additional scripts that run i-PI in a slightly different way than in the original executable, as well as processing scripts for computing the appropriate spectroscopic response function.

Solid-liquid free energy with uncertainty quantification
--------------------------------------------------------

This example shows how to perform an interface-pinning simulation of a liquid/solid coexistence supercell. 
Pinning is achieved using an order parameter computed by Plumed, and the model is a committee of ML potentials.
Statistical reweighting is used to compute the error on the difference in chemical potential between
water and ice at the set temperature. Repeating the simulation at multiple temperatures can be used to determine
the melting point.

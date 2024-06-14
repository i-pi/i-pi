Welcome to i-PI's documentation!
================================

i-PI is a interface for advanced molecular simulations written in
Python 3. The main goal is to decouple the problem of evolving the 
ionic positions and the problem of computing the inter-atomic
forces. i-PI was initially developed for Path Integral Molecular
Dynamics (PIMD) simulations, and contains probably the most
comprehensive array of PIMD techniques. Since v2.0, however, it also
contains several general-purpose methods for molecular simulations,
ranging from phonon calculators to replica exchange molecular dynamics.

.. figure:: ../figures/ipi-scheme.*
   :width: 90.0%

   Schematic representation of the functioning of i-PI.

The implementation is based on a client-server paradigm, where i-PI acts
as the server and deals with the propagation of the nuclear dynamics,
whereas the calculation of the potential energy, forces and the
potential energy part of the pressure virial is delegated to one or more
instances of an external code, acting as clients. The implementation of i-PI is
efficient enough that, by tuning the socket parameters and avoiding 
excessive I/O it can be used to run empirical forcefields or machine learning
interatomic potentials with tens of thousands of atoms with only a small overhead. 

This documentation is structured as follows:

.. toctree::
   :maxdepth: 2

   introduction
   features
   getting-started
   units
   input-files
   input-tags
   output-files
   output-tags
   distributed
   tutorials
   faq
   troubleshooting
   contributing
   bibliography


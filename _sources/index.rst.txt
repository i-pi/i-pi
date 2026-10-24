Welcome to i-PI's documentation!
================================

i-PI is a force engine written in Python 3 with the goal of performing standard and advanced molecular simulations.
The implementation is based on a client-server paradigm, where i-PI acts
as the server and deals with the propagation of the nuclear dynamics,
whereas the calculation of the potential energy, forces and the
potential energy part of the pressure virial is delegated to one or more
instances of an external code, acting as clients.
Thus, i-PI effectively decouples the problem of evolving the 
ionic positions and the problem of computing the system-specific properties

i-PI was initially developed for Path Integral Molecular
Dynamics (PIMD) simulations, and contains probably the most
comprehensive array of PIMD techniques. However, it has evolved over time to become
a general-purpose code with methods ranging from phonon calculators to replica exchange molecular dynamics.


The implementation of i-PI is efficient enough that, by tuning the socket parameters and avoiding 
excessive I/O, it can be used to run simulations using empirical forcefields or machine learning
interatomic potentials with tens of thousands of atoms with only a small overhead. 

This documentation is structured as follows:

.. toctree::
   :maxdepth: 2

   introduction
   getting-started
   features
   units
   input-files
   input-tags
   output-files
   output-tags
   python
   distributed
   tutorials
   onlinereso
   faq
   troubleshooting
   contributing
   bibliography


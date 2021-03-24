Welcome to i-PI's documentation!
================================

i-PI is a interface for advanced molecular simulations written in
Python, designed to be used together with an *ab initio* evaluation of
the interactions between the atoms. The main goal is to decouple the
problem of evolving the ionic positions to sample the appropriate
thermodynamic ensemble and the problem of computing the inter-atomic
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
instances of an external code, acting as clients. Since the main focus
is on performing *ab initio* PIMD – where the cost of the force
evaluation is overwhelming relative to the ionic dynamics – clarity has
been privileged over speed. Still, the implementation of i-PI is
efficient enough that it can be used with empirical forcefields to
perform simple benchmarks and preparatory simulations.

The documentation including the python help strings of the keywords is
best accessed through the latex-generated PDF file that can be compiled from the code
repository for the most recent version or downloaded from here (static
version updated periodically): `manual-pdf`_

.. _manual-pdf: ./_static/manual.pdf

This documentation is structured as follows:

.. toctree::
   :maxdepth: 2

   introduction
   getting-started
   user-guide
   tutorials
   input-reference
   troubleshooting
   faq
   contributing
   bibliography

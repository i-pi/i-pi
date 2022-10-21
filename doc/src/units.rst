Internal units
====================

All the units used internally by i-PI are atomic units, as given below.
By default, both input and output data are given in atomic units, but in
most cases the default units can be overridden if one wishes so. For
details on how to do this, see `3.1.1.1 <#inputunits>`__ and
`3.2.1 <#propertyfile>`__.

.. container:: center

   =========== ============= ================
   Unit        Name          S.I. Value
   =========== ============= ================
   Length      Bohr radius   5.2917721e-11 m
   Time        N.A.          2.4188843e-17 s
   Mass        Electron mass 9.1093819e-31 kg
   Temperature Hartree       315774.66 K
   Energy      Hartree       4.3597438e-18 J
   Pressure    N.A.          2.9421912e13 Pa
   =========== ============= ================

Regarding the specification of these units in the i-PI input files, the
user is able to specify units both in the i-PI input file or in the
structure file. If the structure file is of the .xyz format, the units
specifications should be present in the comment line. Examples of such
inputs can be found in ``examples/pes-regtest/io-units/`` . The code
then behaves in the following way, depending on the user’s choice:

-  If no units are specified in the input, i-PI tries to guess from the
   comment line of the structure file and it nothing is present, assumes
   atomic units, except for the cases discussed below.

-  If units are specified in both input and structure file and they
   match, conversion happens just once. If they do not match, an error
   is raised and i-PI stops.


There are a few exceptions:

  - The total simulation time, specified by the flag <total_time>, has to be provided  in seconds.
  - The PDB standard defines units to be in angstrom, and so default units for PDB are angstrom.  
  - When using ASE format for I/O, the units are dictated by the ASE defaults. 
    So for example, the ASE extended xyz format  expects positions in angstroms.
    (See for example ``examples/ase-io/``)


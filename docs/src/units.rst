Internal units and conventions
==============================

Everything is a.u.
------------------

All the units used internally by i-PI are atomic units, as given below.
*All* of them, even the temperature.
By default, both input and output data are given in atomic units, but in
most cases the default units can be overridden if one wishes so. For
details on how to do this, see :ref:`inputunits` and
:ref:`propertyfile`.

.. container:: center

   =========== ============= ================
   Unit        Name          S.I. Value
   =========== ============= ================
   Length      Bohr radius   5.2917721e-11 m
   Time        Atomic units  2.4188843e-17 s
   Mass        Electron mass 9.1093819e-31 kg
   Temperature Hartree       315774.66 K
   Energy      Hartree       4.3597438e-18 J
   Pressure    Atomic units  2.9421912e13 Pa
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


There are a few exceptions, concerning quantities that are not used to describe the atomic-scale system, or I/O that is constrained by existing standards:

  - Quantities that define run-time execution options, such as the simulation wall-clock time limit (specified by the flag <total_time>) or the latency of the forcefield socket (specified by the parameter <latency>), have to be provided in seconds.
  - The PDB standard defines units to be in angstrom, and so default units for PDB are angstrom.  
  - When using ASE format for I/O, the units are dictated by the ASE defaults. 
    So for example, the ASE extended xyz format  expects positions in angstroms.
    (See for example ``examples/ase-io/``)


Lattice parameters
------------------

The unit cell is defined internally by i-PI with the lattice vectors stored 
in the columns - i.e. as a :math:`3\times 3` array in which ``h[i,j]`` contains the i-th 
Cartesian component of the j-th lattice vector.
Furthermore, the cell is constrained to be oriented with the first vector
along the :math:`x` axis, the second vector within the :math:`xy` plane, and
the third in an arbitrary position (so that the cell is stored internally as an 
upper-triangular matrix). 
This is also reflected in how the cell parameters should be provided in 
the input file: when given in an array format, the correct format is e.g.
``<h units='angstrom'> [10, 1, 2, 0, 9, -1, 0, 0, 11] </h>``. 


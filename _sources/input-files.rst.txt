.. _inputfiles:

Input files
===========

.. _ifilestructure:

Input file format and structure
-------------------------------

In order to give the clearest layout, xml formatting was chosen as the
basis for the main input file. An xml file consists of a set of
hierarchically nested tags. There are three parts to an xml tag. Each
tag is identified by a tag name, which specifies the class or variable
that is being initialized. Between the opening and closing tags there
may be some data, which may or may not contain other tags. This is used
to specify the contents of a class object, or the value of a variable.
Finally tags can have attributes, which are used for metadata, i.e. data
used to specify how the tag should be interpreted. As an example, a
‘mode’ attribute can be used to select between different thermostatting
algorithms, specifying how the options of the thermostat class should be
interpreted.

A xml tag has the following syntax:

The syntax for the different types of tag data is given below:

.. container:: center

   ========== ==========================================
   Data type  Syntax
   ========== ==========================================
   Boolean    <tag>True</tag> or <tag>False</tag>
   Float      <tag>11.111</tag> or <tag>1.1111e+1</tag>
   Integer    <tag>12345</tag>
   String     <tag>string_data</tag>
   Tuple      <tag> (int1, int2, …)</tag>
   Array      <tag> [ entry1, entry2, …] </tag>
   Dictionary <tag>{name1: data1, name2: data2, …}</tag>
   ========== ==========================================

Note that arrays are always given as one-dimensional lists. In cases
where a multi-dimensional array must be entered, one can use the ‘shape’
attribute, that determines how the list will be reshaped into a
multi-dimensional array. For example, the bead positions are represented
in the code as an array of shape (number of beads, 3*number of atoms).
If we take a system with 20 atoms and 8 beads, then this can be
represented in the xml input file as:

.. code-block::

   <beads nbeads=’8’ natoms=’20’> <q shape=‘(8,60)’>[ q11x, q11y, q11z,
   q12x, q12y, ... ]</q> ... </beads>

If ‘shape’ is not specified, a 1D array will be assumed.

The code uses the hierarchical nature of the xml format to help read the
data; if a particular object is held within a parent object in the code,
then the tag for that object will be within the appropriate parent tags.
This is used to make the structure of the simulation clear.

For example, the system that is being studied is partly defined by the
thermodynamic ensemble that should be sampled, which in turn may be
partly defined by the pressure, and so on. To make this dependence clear
in the code the global simulation object which holds all the data
contains an ensemble object, which contains a pressure variable.

Therefore the input file is specified by having a :ref:`simulation` tag, containing an
:ref:`ensemble` tag, which itself contains a :ref:`pressure` tag, which will contain a float
value corresponding to the external pressure. In this manner, the class
structure can be constructed iteratively.

For example, suppose we want to generate a *NPT* ensemble at an external
pressure of :math:`10^{-7}` atomic pressure units. This would be
specified by the following input file:

.. code-block::

   <system> <ensemble> <pressure> 1e-7 </pressure> ... </ensemble> ...
   </system>

To help detect any user error the recognized tag names, data types and
acceptable options are all specified in the code in a specialized input
class for each class of object. A full list of all the available tags
and a brief description of their function is given at 
`input tags <input-tags.html>`_.

.. _inputunits:

Overriding default units
~~~~~~~~~~~~~~~~~~~~~~~~

Many of the input parameters, such as the pressure in the above example,
can be specified in more than one unit. Indeed, often the atomic unit is
inconvenient to use, and we would prefer something else. Let us take the
above example, but instead take an external pressure of 3 MPa. Instead
of converting this to the atomic unit of pressure, it is possible to use
pascals directly using:

.. code-block::

   <system> <ensemble> <pressure units=‘pascal’> 3e6 </pressure> ...
   </ensemble> ... </system>

The code can also understand S.I. prefixes, so this can be simplified
further using:

.. code-block::

   <system> <ensemble> <pressure units=‘megapascal’> 3 </pressure> ...
   </ensemble> ... </system>

A full list of which units are defined for which dimensions can be found
in the units.py module.

Initialization section
----------------------

The input file can contain a :ref:`initializer` tag, which contains a number of fields that
determine the starting values of the various quantities that define the
state of the simulation – atomic positions, cell parameters, velocities,
…. These fields (:ref:`positions`, :ref:`velocities`, :ref:`cell`, 
:ref:`masses`, :ref:`labels`, :ref:`file` ) specify how the values should be obtained:
either from a manually-input list or from an external file.

.. _configfile:

Configuration files
~~~~~~~~~~~~~~~~~~~

Instead of initializing the atom positions manually, the starting
configuration can be specified through a separate data file. The name of
the configuration file is specified within one of the possible fields of
an :ref:`initializer` tag. The file format is specified with the “mode” attribute. The
currently accepted file formats are:

-  pdb

-  xyz

-  chk

the last of which will be described in the next section.

Depending on the field name, the values read from the external file will
be used to initialize one component of the simulation or another (e.g.
the positions or the velocities). The :ref:`initfile` tag can be used as a shortcut to
initialize the atom positions, labels, masses and possibly the cell
parameters at the same time. For instance,

.. code-block::

   <initialize nbeads="8"> <file mode="pdb"> init.pdb </file>
   </initialize>

is equivalent to

.. code-block::

   <initialize nbeads="8"> <positions mode="pdb"> init.pdb </positions>
   <labels mode="pdb"> init.pdb </labels> <masses mode="pdb"> init.pdb
   </masses> <cell mode="pdb"> init.pdb </cell> </initialize>

In practice, the using the :ref:`initfile` tag will only read the information that can
be inferred from the given file type, so for an ‘xyz’ file that does not contain a cell, the cell
parameters will not be initialized.

Initialization from checkpoint files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

i-PI gives the option to output the entire state of the simulation at a
particular timestep as an xml input file, called a checkpoint file (see
:ref:`checkpoints` for details). As well as being a valid input for
i-PI, a checkpoint can also be used inside an :ref:`initializer` tag to specify the
configuration of the system, discarding other parameters of the
simulation such as the current time step or the chosen ensemble. Input
from a checkpoint is selected by using “chk” as the value of the “mode”
attribute. As for the configuration file, a checkpoint file can be used
to initialize either one or many variables depending on which tag name
is used.

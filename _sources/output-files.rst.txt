.. _outputfiles:

Output files
============

i-PI uses a very flexible mechanism to specify how and how often atomic
configurations and physical properties should be output. Within the :ref:`outputs` tag
of the xml input file the user can specify multiple tags, each one of
which will correspond to a particular output file. Each file is managed
separately by the code, so what is output to a particular file and how
often can be adjusted for different files independently.

For example, some of the possible output properties require more than
one force evaluation per time step to calculate, and so can considerably
increase the computational cost of a simulation unless they are computed
once every several time steps. On the other hand, for properties such as
the conserved energy quantity it is easy, and often useful, to output
them every time step as they are simple to compute and do not take long
to output to file.

There are three types of output file that can be specified; property
files for system level properties, trajectory files for atom/bead level
properties, and checkpoint files which save the state of the system and
so can be used to restart the simulation from a particular point. For a
brief overview of the format of each of these types of files, and some
of their more common uses, see :ref:`part1`. To give a more in depth
explanation of each of these files, they will now be considered in turn.

.. _propertyfile:

Properties
----------

This is the output file for all the system and simulation level
properties, such as the total energy and the time elapsed. It is
designed to track a small number of important properties throughout a
simulation run, and as such has been formatted to be used as input for
plotting programs such as gnuplot.

The file starts with a header, which describes the properties being
written in the different columns and their output units. This is
followed by the actual data. Each line corresponds to one instant of the
simulation. The file is fixed formatted, with two blank characters at
the start of each row, then the data in the same order as the header
row. By default, each column is 16 characters wide and every float is
written in exponential format with 8 digits after the decimal point.

For example, if we had asked for the current time step, the total
simulation time in picoseconds, and the potential energy in
electronvolt, then the properties output file would look something like:

.. code-block::

   # column 1 –> step : The current simulation time step. # column 2 –>
   timepicosecond : The elapsed simulation time. # column 3 –>
   potential{electronvolt} : The physical system potential energy.
   0.00000000e+00 0.00000000e+00 -1.32860475e+04 1.00000000e+00
   1.00000000e-03 -1.32865789e+04 ...

The properties that are output are determined by the :ref:`properties` tag in the xml
input file. The format of this tag is:

.. code-block::

   <properties stride= filename= flush= shape=> [ prop1name{units}(arg1;
   ... ), prop2name{units}(...), ... ] </properties>

e.g.

The attributes have the following meanings:

stride
   The number of steps between each output to file

filename
   The name of the output file

flush
   The number of output lines between buffer flushes

shape
   The number of properties in the list.

The tag data is an array of strings, each of which contains three
different parts:

-  The property name, which describes which type of property is to be
   output. This is a mandatory part of the string.

-  The units that the property will be output in. These are specified
   between curly brackets. If this is not specified, then the property
   will be output in atomic units. Note that some properties can only be
   output in atomic units.

-  The arguments to be passed to the function. These are specified
   between standard brackets, with each argument separated by a
   semi-colon. These may or may not be mandatory depending on the
   property, as some arguments have well defined default values. The
   arguments can be specified by either of two different syntaxes,
   (name1=arg1; …) or (arg1; …).

   The first syntax uses keyword arguments. The above example would set
   the variable with the name “name1” the value “arg1”. The second
   syntax uses positional arguments. This syntax relies on the arguments
   being specified in the correct order, as defined in the relevant
   function in the property.py module, since the user has not specified
   which variable to assign the value to.

   The two syntaxes may be mixed, but positional arguments must be
   specified first otherwise undefined behaviour will result. If no
   arguments are specified, then the defaults as defined in the
   properties.py module will be used.

See also the documentation of the :ref:`properties` tag, and
the full :ref:`property_list`.

.. _trajectories:

Trajectory files
----------------

These are the output files for atomic or bead level properties, such as
the bead positions. In contrast to properties files, they output data
for all atomic degrees of freedom, in a format that can be read by
visualization packages such as VMD.

Multiple trajectory files can be specified, each described by a separate :ref:`trajectory`
tag within the :ref:`outputs` section of the input file. The allowable file formats for
the trajectory output files are the same as for the configuration input
files, given in :ref:`configfile`.

These tags have the format:

.. code-block::

    <trajectory stride=`' filename=`' format=`' cell_units=`' flush=`' bead=`'>
         traj_name{units}(arg1;...)
    </trajectory>

This is very similar to the :ref:`properties` tag, except that it has the additional 
tags “format” and “cell_units”, and only one *traj_name* quantity can be
specified per file. ‘format’ specifies the format of the output file,
and ‘cell_units’ specifies the units in which the cell dimensions are
output. Depending on the quantity being output, the trajectory may
consist of just one file per time step (e.g. the position of the
centroid) or of several files, one for each bead, whose name will be
automatically determined by appending the bead index to the specified
“filename” attribute (e.g. the beads position). In the latter case it is
also possible to output the quantity computed for a single bead by
specifying its (zero-based) index in the “bead” attribute.

See also the :ref:`trajectory` tag, and the full
:ref:`trajectory_list`. 

.. _checkpoints:

Checkpoint files
----------------

As well as the above output files, the state of the system at a
particular time step can also be saved to file. These checkpoint files
can later be used as input files, with all the information required to
restore the state of the system to the point at which the file was
created.

This is specified by the :ref:`checkpoint` tag which has the syntax:

.. code-block::

   <checkpoint stride= filename= overwrite=> step </checkpoint>

Again, this is similar to the :ref:`properties` and 
:ref:`trajectory` tags, but instead of having a value
which specifies what to output, the value simply gives a number to
identify the current checkpoint file. There is also one additional
attribute, “overwrite”, which specifies whether each new checkpoint file
overwrites the old one, or whether all checkpoint files are kept. If
they are kept, they will be written not to the file “filename”, but
instead an index based on the value of “step” will be appended to it to
distinguish between different files.

If the ‘step’ parameter is not specified, the following syntax can also
be used:

.. code-block::

   <checkpoint stride= filename= overwrite=/>

Soft exit and RESTART
~~~~~~~~~~~~~~~~~~~~~

As well as outputting checkpoint files during a simulation run, i-PI
also creates a checkpoint automatically at the end of the simulation,
with file name “RESTART”. In the same way as the checkpoint files
discussed above, it contains the full state of the simulation. It can be
used to seamlessly restart the simulation if the user decides that a
longer run is needed to gather sufficient statistics, or if i-PI is
terminated before the desired number of steps have been completed.

i-PI will try to generate a RESTART file when it terminates, either
because *total_time* has elapsed, or because it received a (soft) kill
signal by the operating system. A soft exit can also be forced by
creating an empty file named “EXIT” in the directory in which i-PI is
running.

An important point to note is that since each time step is split into
several parts, it is only at the end of each step that all the variables
are consistent with each other in such a way that the simulation can be
restarted from them without changing the dynamics. Thus if a soft exit
call is made during a step, then the restart file that is created must
correspond to the state of the system *before* that step began. To this
end, the state of the system is saved at the start of every step.

In order to restart i-PI from a file named “RESTART" one simply has to
run

.. code-block::

   > python i-pi RESTART


Reading output files
--------------------

It can be useful to parse the output files of i-PI into a format that can
be readily manipulated in a custom Python script. To this end, i-PI provides
a few utilities. `ipi.read_output` that can be used to parse a property output 
file, that returns data blocks as a dictionary of numpy array, and additional
information from the header (such as units, and the description of each
property) as a separate dictionary

.. code-block::

   from ipi import read_output
   data, info = read_output("simulation.out")

Trajectory files can be read with `ipi.read_trajectory`. This reads the 
trajectory output into a list of `ase.Atoms` objects (hence this functionality
has a dependency on `ase`), converting positions and cell to angstrom, and 
moving other properties to arrays and converting 
them to ASE units (e.g. forces are in eV/Å). `extras` output files (that
contain either numerical data, or raw strings, returned by the driver code
in addition to energy and forces) can be processed by using `format='extras'`
and an option: the call will then return a dictionary with an entry having the
name of the type of extra (if present) and either a list of the raw strings, 
or a numpy array with the data. A second dictionary entry contains the list
of step numbers. 

.. code-block::

   from ipi import read_trajectory
   data = read_trajectory("simulation.dipoles", format="extras")


User guide
==========

Input files
-----------

.. _ifilestructure:

Input file format and structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

Therefore the input file is specified by having a tag, containing an
tag, which itself contains a “pressure” tag, which will contain a float
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
and a brief description of their function is given in
chapter `4 <#hierarchy>`__.

.. _inputunits:

Overriding default units
^^^^^^^^^^^^^^^^^^^^^^^^

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
~~~~~~~~~~~~~~~~~~~~~~

The input file can contain a tag, which contains a number of fields that
determine the starting values of the various quantities that define the
state of the simulation – atomic positions, cell parameters, velocities,
…. These fields (, , , , , ) specify how the values should be obtained:
either from a manually-input list or from an external file.

.. _configfile:

Configuration files
^^^^^^^^^^^^^^^^^^^

Instead of initializing the atom positions manually, the starting
configuration can be specified through a separate data file. The name of
the configuration file is specified within one of the possible fields of
an tag. The file format is specified with the “mode” attribute. The
currently accepted file formats are:

-  pdb

-  xyz

-  chk

the last of which will be described in the next section.

Depending on the field name, the values read from the external file will
be used to initialize one component of the simulation or another (e.g.
the positions or the velocities). The tag can be used as a shortcut to
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

In practice, the using the tag will only read the information that can
be inferred from the given file type, so for an ‘xyz’ file, the cell
parameters will not be initialized.

Initialization from checkpoint files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

i-PI gives the option to output the entire state of the simulation at a
particular timestep as an xml input file, called a checkpoint file (see
`3.2.3 <#checkpoint>`__ for details). As well as being a valid input for
i-PI, a checkpoint can also be used inside an tag to specify the
configuration of the system, discarding other parameters of the
simulation such as the current time step or the chosen ensemble. Input
from a checkpoint is selected by using “chk” as the value of the “mode”
attribute. As for the configuration file, a checkpoint file can be used
to initialize either one or many variables depending on which tag name
is used.

.. _outputfiles:

Output files
------------

i-PI uses a very flexible mechanism to specify how and how often atomic
configurations and physical properties should be output. Within the tag
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
of their more common uses, see `5.1 <#part1>`__. To give a more in depth
explanation of each of these files, they will now be considered in turn.

.. _propertyfile:

Properties
~~~~~~~~~~

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
   potentialelectronvolt : The physical system potential energy.
   0.00000000e+00 0.00000000e+00 -1.32860475e+04 1.00000000e+00
   1.00000000e-03 -1.32865789e+04 ...

The properties that are output are determined by the tag in the xml
input file. The format of this tag is:

.. code-block::

   <properties stride= filename= flush= shape=> [ prop1nameunits(arg1;
   ... ), prop2name...(...), ... ] </properties>

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

The different available properties are:

.. _trajectories:

Trajectory files
~~~~~~~~~~~~~~~~

These are the output files for atomic or bead level properties, such as
the bead positions. In contrast to properties files, they output data
for all atomic degrees of freedom, in a format that can be read by
visualization packages such as VMD.

Multiple trajectory files can be specified, each described by a separate
tag within the section of the input file. The allowable file formats for
the trajectory output files are the same as for the configuration input
files, given in `3.1.2.1 <#configfile>`__.

These tags have the format:

This is very similar to the tag, except that it has the additional tags
“format” and “cell_units”, and only one *traj_name* quantity can be
specified per file. ‘format’ specifies the format of the output file,
and ‘cell_units’ specifies the units in which the cell dimensions are
output. Depending on the quantity being output, the trajectory may
consist of just one file per time step (e.g. the position of the
centroid) or of several files, one for each bead, whose name will be
automatically determined by appending the bead index to the specified
“filename” attribute (e.g. the beads position). In the latter case it is
also possible to output the quantity computed for a single bead by
specifying its (zero-based) index in the “bead” attribute.

The quantities that can be output in trajectory files are:

.. _checkpoint:

Checkpoint files
~~~~~~~~~~~~~~~~

As well as the above output files, the state of the system at a
particular time step can also be saved to file. These checkpoint files
can later be used as input files, with all the information required to
restore the state of the system to the point at which the file was
created.

This is specified by the tag which has the syntax:

.. code-block::

   <checkpoint stride= filename= overwrite=> step </checkpoint>

Again, this is similar to the and tags, but instead of having a value
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
^^^^^^^^^^^^^^^^^^^^^

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

.. _distrib:

Distributed execution
---------------------

.. _communication-protocol-1:

Communication protocol
~~~~~~~~~~~~~~~~~~~~~~

i-PI is based on a clear-cut separation between the evolution of the
nuclear coordinates and the evaluation of energy and forces, which is
delegated to an external program. The two parts are kept as independent
as possible, to minimize the client-side implementation burden, and to
make sure that the server will be compatible with any empirical or *ab
initio* code that can compute inter-atomic forces for a given
configuration.

Once a communication channel has been established between the client and
the server (see `3.3.3 <#sockets>`__), the two parties exchange minimal
information: i-PI sends the atomic positions and the cell parameters to
the client, which computes energy, forces and virial and returns them to
the server.

The exchange of information is regulated by a simple communication
protocol. The server polls the status of the client, and when the client
signals that is ready to compute forces i-PI sends the atomic positions
to it. When the client responds to the status query by signalling that
the force evaluation is finished, i-PI will prepare to receive the
results of the calculation. If at any stage the client does not respond
to a query, the server will wait and try again until a prescribed
timeout period has elapsed, then consider the client to be stuck,
disconnect from it and reassign its force evaluation task to another
active instance. The server assumes that 4-byte integers, 8-byte floats
and 1-byte characters are used. The typical communication flow is as
follows:

#. a header string “**STATUS**” is sent by the server to the client that
   has connected to it;

#. a header string is then returned, giving the status of the client
   code. Recognized messages are:

   “NEEDINIT”:
      if the client code needs any initialising data, it can be sent
      here. The server code will then send a header string “INIT”,
      followed by an integer corresponding to the bead index, another
      integer giving the number of bits in the initialization string,
      and finally the initialization string itself.

   “READY”:
      sent if the client code is ready to calculate the forces. The
      server socket will then send a string “POSDATA”, then nine floats
      for the cell vector matrix, then another nine floats for the
      inverse matrix. The server socket will then send one integer
      giving the number of atoms, then the position data as 3 floats for
      each atom giving the 3 cartesian components of its position.

   “HAVEDATA”:
      is sent if the client has finished computing the potential and
      forces. The server socket then sends a string “GETFORCE”, and the
      client socket returns “FORCEREADY”. The potential is then returned
      as a float, the number of atoms as an integer, then the force data
      as 3 floats per atom in the same way as the positions, and the
      virial as 9 floats in the same way as the cell vector matrix.
      Finally, the client may return an arbitrary string containing
      additional data that have been obtained by the electronic
      structure calculation (atomic charges, dipole moment, …). The
      client first returns an integer specifying the number of
      characters, and then the string, which will be output verbatim if
      this “extra” information is requested in the output section (see
      `3.2.2 <#trajectories>`__). The string can be formatted in the
      JSON format, in which case i-PI can extract and process individual
      fields, that can be printed separately to different files.

#. The server socket waits until the force data for each replica of the
   system has been calculated and returned, then the MD can be
   propagated for one more time step, and new force requests will be
   dispatched.

Parallelization
~~~~~~~~~~~~~~~

As mentioned before, one of the primary advantages of using this type of
data transfer is that it allows multiple clients to connect to an i-PI
server, so that different replicas of the system can be assigned to
different client codes and their forces computed in parallel. In the
case of *ab initio* force evaluation, this is a trivial level of
parallelism, since the cost of the force calculation is overwhelming
relative to the overhead involved in exchanging coordinates and forces.
Note that even if the parallelization over the replicas is trivial,
often one does not obtain perfect scaling, due to the fact that some of
the atomic configurations might require more steps to reach
self-consistency, and the wall-clock time per step is determined by the
slowest replica.

i-PI maintains a list of active clients, and distributes the forces
evaluations among those available. This means that, if desired, one can
run an :math:`n`-bead calculation using only :math:`m<n` clients, as the
server takes care of sending multiple replicas to each client per MD
step. To avoid having clients idling for a substantial amount of time,
:math:`m` should be a divisor of :math:`n`. The main advantage of this
approach, compared to one that rigidly assigns one instance of the
client to each bead, is that if each client is run as an independent job
in a queue (see `2.3.3 <#hpc>`__), i-PI can start performing PIMD as
soon as a single job has started, and can carry on advancing the
simulation even if one of the clients becomes unresponsive.

Especially for *ab initio* calculations, there is an advantage in
running with :math:`m=n`. i-PI will always try to send the coordinates
for one path integral replica to the client that computed it at the
previous step: this reduces the change in the particle positions between
force evaluations, so that the charge density/wavefunction from the
previous step is a better starting guess and self-consistency can be
achieved faster. Also, receiving coordinates that represent a continuous
trajectory makes it possible to use extrapolation strategies that might
be available in the client code.

Obviously, most electronic-structure client codes provide a further
level of parallelisation, based on OpenMP and/or MPI. This is fully
compatible with i-PI, as it does not matter how the client does the
calculation since only the forces, potential and virial are sent to the
server, and the communication is typically performed by the master
process of the client.

Sockets
~~~~~~~

The communication between the i-PI server and the client code that
evaluates forces is implemented through sockets. A socket is a data
transfer device that is designed for internet communication, so it
supports both multiple client connections to the same server and two-way
communication. This makes sockets ideal for use in i-PI, where each
calculation may require multiple instances of the client code. A socket
interface can actually function in two different modes.

UNIX-domain sockets are a mechanism for local, inter-process
communication. They are fast, and best suited when one wants to run i-PI
with empirical potentials, and the latency of the communication with the
client becomes a significant overhead for the calculation. UNIX-domain
sockets create a special file in the local file system, that serves as a
rendezvous point between server and clients, and are uniquely identified
by the name of the file itself, that can be specified in the “address”
tag of in the xml input file and in the input of the client.

Unfortunately, UNIX sockets do not allow one to run i-PI and the clients
on different computers, which limits greatly their utility when one
needs to run massively parallel calculations. In these cases – typically
when performing *ab initio* simulations – the force calculation becomes
the bottleneck, so there is no need for fast communication with the
server, and one can use internet sockets, that instead are specifically
designed for communication over a network.

Internet sockets are described by an address and a port number. The
address of the host is given as the IP address, or as a hostname that is
resolved to an IP address by a domain name server, and is specified by
the “address” variable of a object. The port number is an integer
between 1 and 65535 used to distinguish between all the different
sockets open on a particular host. As many of the lower numbers are
protected for use in important system processes or internet
communication, it is generally advisable to only use numbers in the
range 1025-65535 for simulations.

The object has two more parameters. The option “latency” specifies how
often i-PI polls the list of active clients to dispatch positions and
collect results: setting it to a small value makes the program more
responsive, which is appropriate when the evaluation of the forces is
very fast. In *ab initio* simulations, it is best to set it to a larger
value (of the order of 0.01 seconds), as higher latency will have no
noticeable impact on performance, but will reduce the cost of having
i-PI run in the background to basically zero.

Normally, i-PI can detect when one of the clients dies or disconnects,
and can remove it from the active list and dispatch its force
calculation to another instance. If however one of the client hangs
without closing the communication channel, i-PI has no way of
determining that something is going wrong, and will just wait forever.
One can specify a parameter “timeout”, that corresponds to the maximum
time – in seconds – that i-PI should wait before deciding that one of
the clients has become unresponsive and should be discarded.

Running i-PI over the network
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Understanding the network layout
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Running i-PI in any non-local configuration requires a basic
understanding of the layout of the network one is dealing with. Each
workstation, or node of a HPC system, may expose more than one network
interface, some of which can be connected to the outside internet, and
some of which may be only part of a local network. A list of the network
interfaces available on a given host can be obtained for instance with
the command

.. code-block::

   > /sbin/ip addr

which will return a list of interfaces of the form

Each item corresponds to a network interface, identified by a number and
a name (lo, eth0, eth1, …). Most of the interfaces will have an
associated IP address – the four numbers separated by dots that are
listed after “inet”, e.g. 192.168.1.254 for the eth0 interface in the
example above.

.. figure:: ../figures/ipi-network.*
   :width: 90.0%

   A schematic representation of the network layout one
   typically finds when running i-PI and the clients on a HPC system
   and/or on a local workstation.

Figure `3.1 <#fig:network>`__ represents schematically a typical network
layout for a HPC system and a local workstation. When running i-PI
locally on a workstation, one can use the loopback interface (that can
be referred to as “localhost” in the “address” field of both i-PI and
the client) for communication. When running both i-PI and the clients on
a HPC cluster, one should work out which of the the interfaces that are
available on the node where the i-PI server runs are accessible from the
compute nodes. This requires some trial and error, and possibly setting
the “address” field dynamically from the job that launches i-PI. For
instance, if one was running i-PI on the login node, and the clients on
different compute nodes, as in Figure `2.1 <#fig:running>`__\ b, then on
the HPC system described in Figure `3.1 <#fig:network>`__ one should set
the address to that of the *ib1* interface – :math:`111.111.111.111` in
the example above. If instead i-PI was launched in a job script, then
the submission script would have to check for the IP address associated
with the *ib0* interface on the node the job has been dispatched to, and
set that address (e.g. :math:`111.111.111.200`) in the inputs of both
i-PI and the clients that will be launched in the same (or separate)
jobs.

Running i-PI on a separate workstation
(Figure `2.1 <#fig:running>`__\ c) gives maximum flexibility, but is
also trickier as one has to reach the internet from the compute nodes,
that are typically not directly connected to it. We discuss this more
advanced setup in the next paragraph.

.. _ssh_sockets:

ssh tunnelling
^^^^^^^^^^^^^^

If i-PI is to be run in a distributed computing mode, then one should
make sure that the workstation on which the server will run is
accessible from the outside internet on the range of ports that one
wants to use for i-PI. There are ways to circumvent a firewall, but we
will not discuss them here, as the whole point of i-PI is that it can be
run on a low-profile PC whose security does not need to be critical.
Typically arrangements can be made to open up a range of ports for
incoming connections.

A more substantial problem – as it depends on the physical layout of the
network rather than on software settings of the firewall – is how to
access the workstation from the compute nodes, which in most cases do
not have a network interface directly connected to the outside internet.

The problem can be solved by creating a ssh tunnel, i.e. an instance of
the ssh secure shell that will connect the compute node to the login
node, and then forward all traffic that is directed to a designated port
on the compute node to the remote location that is running i-PI, passing
through the outbound network interface of the login node.

In the example above, if i-PI is running on a local workstation, one
should run:

from the job script that launches the client. For instance, with the
network layout of Figure `3.1 <#fig:network>`__, and if the i-PI server
is listening on port 12345 of the *eth0* interface, the tunnel should be
created as:

.. code-block::

   > ssh -f -N -L 54321:123.123.123.123:12345 -2 111.111.111.111

The client should then be configured to connect to *localhost* on port
54321. The connection with i-PI will be established through the tunnel,
and the data exchange can begin.

Note that, in order to be able to include the above commands in a
script, the login node and the compute nodes should be configured to
allow password-less login within the HPC system. This can be achieved
easily, and does not entail significant security risks, since it only
allows one to connect from one node to another within the local network.
To do so, you should log onto the HPC system, and create a pair of ssh
keys (if this has not been done already, in which case an id_rsa.pub
file should be present in the user’s ~/.ssh/ directory) by issuing the
command

.. code-block::

   > ssh-keygen -t rsa

The program will then prompt for a passphrase twice. Since we wish to
have use this in a job script where we will not be able to enter a
password, just hit enter twice.

This should now have created two files in the directory ~/.ssh, id_rsa
and id_rsa.pub. These should be readable only by you, so use the
following code to set up the correct file permissions:

Finally, copy the contents of the file id_rsa.pub and append them to the
file authorized_keys in the directory ~/.ssh of the user on the login
node, which is typically shared among all the nodes of a cluster and
therefore allows password-less login from all of the compute nodes.

.. _distrib:

Distributed execution
=====================

.. _transports:

Choosing how to talk to the driver
----------------------------------------

i-PI delegates the evaluation of energy, forces and virial to a *driver* 
(effectively connected as a *client* to the i-PI server). 
Several transport mechanisms are available, selected by the
type of ``<forcefield>`` block in the input and, for sockets, by its ``mode``
attribute. They all provide the same information to the rest of i-PI and differ
only in *how* the data is moved between the server and the driver.

.. list-table::
   :header-rows: 1
   :widths: 28 20 52

   * - Forcefield block
     - Where the driver runs
     - Typical use
   * - :ref:`ffdirect`
     - in-process (same Python interpreter)
     - Python potentials from ``ipi/pes`` (e.g. ML potentials); no inter-process
       communication, lowest overhead for inexpensive potentials.
   * - ``<ffsocket mode="unix">``
     - same machine
     - empirical/fast potentials run locally; low-latency UNIX-domain socket.
   * - ``<ffsocket mode="shm">``
     - same machine
     - as ``unix``, but the bulk position/force payload travels through shared
       memory — preferable for large systems.
   * - ``<ffsocket mode="inet">``
     - any machine on the network
     - *ab initio* clients, possibly on a different host or HPC cluster.
   * - :ref:`ffmpi`
     - co-launched MPI ranks
     - HPC runs that launch i-PI and the drivers in a single MPI job, or drivers
       that are themselves MPI-parallel.

Two further options apply on top of the socket and MPI transports:

* **Multiple clients** — several drivers can connect to the same i-PI server and
  evaluate different replicas (path-integral beads, or systems) in parallel; see
  :ref:`parallelization`.
* **Batching** (``batch_size``) — a single driver can evaluate several queued
  structures in one request, reducing the per-message overhead; see
  :ref:`batching`.

.. note::

   The transport only matters when the force evaluation itself is cheap (toy
   potentials, some ML models, or very large bead counts). For *ab initio*
   clients the cost of the force calculation dwarfs any communication overhead,
   and the choice should be driven by convenience (network reachability, how the
   job is launched) rather than by performance. As a rough guide in the cheap-PES
   regime: ``shm`` is the best socket transport for large systems, ``batch_size``
   helps whenever many replicas are queued, and ``ffdirect`` avoids communication
   altogether for Python potentials.

.. _communication-protocol-1:

Communication protocol
----------------------

i-PI is based on a clear-cut separation between the evolution of the
nuclear coordinates and the evaluation of energy and forces, which is
delegated to an external program. The two parts are kept as independent
as possible, to minimize the client-side implementation burden, and to
make sure that the server will be compatible with any empirical or *ab
initio* code that can compute inter-atomic forces for a given
configuration.

Once a communication channel has been established between the client and
the server, the two parties exchange minimal
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
      :ref:`trajectories`). The string can be formatted in the
      JSON format, in which case i-PI can extract and process individual
      fields, that can be printed separately to different files.

#. The server socket waits until the force data for each replica of the
   system has been calculated and returned, then the MD can be
   propagated for one more time step, and new force requests will be
   dispatched.


.. figure:: ../figures/ipi-comm.*
   :width: 90.0%

   A schematic simplified representation of the communication protocol.
   We note that most of the clients do not make use of the 'NEEDINIT' option,
   The communication is *asynchronous* but we have omitted the 'waiting' blocks
   for simplicity.

The description above is the *base* protocol. Two variants extend it: when
``batch_size`` is larger than one (see :ref:`batching`) the POSDATA and FORCEREADY
messages carry several structures at once; and in the ``shm`` mode (see
:ref:`shm-sockets`) the bulk arrays travel through shared memory rather than the
socket. The ``extras`` string returned after the forces can be any text, and is
most useful as a JSON dictionary (for example ``{"dipole": [...]}``), from which
individual fields can be selected for output (see :ref:`trajectories`). The
bundled Python driver, ``drivers/py/driver.py``, is a compact reference
implementation of all these variants and the best starting point for writing a
custom client.

.. _parallelization:

Parallelization
---------------

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
in a queue (see :ref:`hpc`), i-PI can start performing PIMD as
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
server, and the communication is typically performed by the main 
process of the client.

Sockets
-------

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
tag of in the xml input file and in the input of the client. By default
this file is created based on the address tag, with a `/tmp/ipi_` prefix.
This can be overridden setting the “sockets_prefix” attribute for the
:ref:`simulation` tag in the input file, or on the command-line using the
`-S` option. Note that several clients do not support changing the default
prefix.

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

.. _shm-sockets:

Shared-memory sockets (``mode="shm"``)
--------------------------------------------

When the driver runs on the same node as i-PI, the ``shm`` mode reduces the cost
of moving large position and force arrays. A UNIX-domain socket is still used for
the control handshake (status queries, the INIT string, the ``extras`` string),
but the bulk communication — positions, cell, forces and virial — is exchanged
through POSIX shared-memory segments owned by i-PI, avoiding a copy through the
socket. On the i-PI side this only requires choosing the mode::

   <ffsocket mode="shm">
      <address> driver </address>
   </ffsocket>

The client must be started in shared-memory mode: the bundled Python and Fortran
drivers both accept the ``--shm`` flag (see :ref:`runningclients`). Because it
relies on shared memory, ``shm`` only works when the server and the client run on
the same node; it is the recommended local transport for large systems, where the
position/force arrays dominate the communication cost.

.. _batching:

Batching and message consolidation
----------------------------------------

By default i-PI dispatches one structure per request and, for path-integral or
multiple-system runs, relies on several clients to evaluate the replicas in
parallel. An alternative is to let a *single* client evaluate several queued
structures at once, by setting ``batch_size`` greater than one::

   <ffsocket mode="unix">
      <address> driver </address>
      <batch_size> 8 </batch_size>
   </ffsocket>

i-PI then collects up to ``batch_size`` queued structures and sends them together;
the client evaluates them in a single call and returns all the results at once.
The batch size is announced to the client in the INIT string (as ``batch_size:N``),
and the bundled drivers pass the structures to the potential through the
``driver(cell_list, pos_list)`` interface. Batching removes most of the
per-message and per-request overhead, which is the dominant cost when the potential
is inexpensive; for expensive potentials it makes little difference. If fewer than
``batch_size`` structures are queued, the socket transport pads the batch by
replicating the last structure, so the run never stalls. The in-process
``ffdirect`` backend does not pad; instead, when run with ``threaded='True'`` it
flushes an incomplete final batch (evaluating fewer than ``batch_size`` structures)
once its poll loop goes idle. Without threading there is no idle flush, so the
number of structures evaluated must be a multiple of ``batch_size`` or the run
hangs.

Batching requires ``consolidate_messages`` (the default) and is incompatible with
``matching="lock"`` and with stateful clients that rely on receiving a continuous
trajectory for a given replica. ``consolidate_messages`` itself fuses the
STATUS/POSDATA/GETFORCE exchange into a single send and uses one thread to collect
the replies, lowering latency; it assumes clients follow the protocol strictly, so
turn it off only for clients that deviate from it.

Whether to use *many clients* or *one batched client* depends on the setup. Many
independent clients give true parallelism and resilience — i-PI keeps advancing the
simulation if one of them dies — and are the natural choice for *ab initio* codes
run as separate jobs. A single batched client minimises communication overhead and
is convenient for fast Python or compiled potentials, especially with large bead
counts.

.. _ffmpi-transport:

MPI transport (``ffmpi``)
-------------------------------

On HPC systems it is often more convenient to exchange data over MPI than over a
socket: there is no address or port to manage, the native high-speed interconnect
is used, and the drivers are launched by the same job scheduler as i-PI. The
:ref:`ffmpi` forcefield does exactly this::

   <ffmpi mode="mpi">
      <batch_size> 8 </batch_size>
   </ffmpi>

i-PI and the drivers are launched **together in a single MPI job**, with i-PI as
rank 0 of ``MPI_COMM_WORLD`` and the remaining ranks acting as drivers. With a
plain MPMD ``mpirun`` this reads::

   mpirun -n 1 i-pi input.xml : -n 4 i-pi-py_driver --mpi -m harmonic -o 1

How the processes are co-launched is left to the submit script — the colon syntax
above, ``srun --multi-prog``, heterogeneous-job steps, and so on — so that
platform-specific MPI quirks stay out of the i-PI input. i-PI discovers the driver
ranks automatically at start-up, and the number of drivers is set purely by the
launch command.

Everything that works over sockets works over MPI: multiple drivers evaluate
replicas in parallel, ``batch_size`` bundles several structures per request, and
the ``extras`` string (e.g. a dipole) is returned in the same way. In addition,
each driver may itself be an MPI-parallel program: passing ``--group-size G``
(``-g``) to the driver splits the driver ranks into groups of ``G`` ranks, where
only the root of each group communicates with i-PI and broadcasts the work to its
group.

The Python driver supports MPI out of the box (``i-pi-py_driver --mpi``); the
Fortran driver must be built with MPI support (``make -C drivers/f90 MPI=1``) and
then accepts ``--mpi``. Worked examples are in
``examples/clients/communication/mpi/harmonic`` and
``examples/clients/communication/mpi/water-dipole``.

Several ``<ffmpi>`` forcefields can run in the same MPMD job. Give each one a
distinct ``name`` and launch its drivers with the matching ``--mpi-name``; i-PI
routes the ranks to the right forcefield automatically::

   mpirun -n 1 i-pi input.xml \
     : -n 4 i-pi-py_driver --mpi --mpi-name solvent -m pesA \
     : -n 2 i-pi-py_driver --mpi --mpi-name solute  -m pesB

A single, unnamed ``<ffmpi>`` still works with untagged drivers (no
``--mpi-name``). When several forcefields each wrap an *MPI-parallel* driver
(``--group-size>1``), also pass a distinct ``--mpi-id`` per forcefield so a
multi-rank group never crosses a forcefield boundary. See
``examples/clients/communication/mpi/two-drivers``. As with sockets, i-PI sends an INIT string
(the ``<parameters>`` plus the batch size) to each driver, so a custom MPI client
can receive initialization data the same way a socket client does.

.. note::

   MPI ranks waiting for work typically *busy-wait*, so launching more driver
   ranks than there are cores available on a node can badly degrade performance
   (the idle drivers starve i-PI and one another). Keep the total number of ranks
   within the core count, or use ``batch_size`` with a single driver rank instead
   of one rank per bead, if the driver can evaluate multiple structures in parallel.

.. _ffdirect-eval:

In-process evaluation (``ffdirect``)
------------------------------------------

When the potential is implemented in Python it can be evaluated directly inside the
i-PI process, with no socket or MPI communication at all, using the :ref:`ffdirect`
forcefield::

   <ffdirect>
      <pes> harmonic </pes>
      <parameters> { k1: 1.0 } </parameters>
   </ffdirect>

The ``pes`` field selects one of the potentials in ``ipi/pes`` (the same classes
used by the Python driver, see :ref:`runningclients`), or ``custom`` together with
``pes_path`` pointing to your own ``.py`` file. ``ffdirect`` also supports
``batch_size``.

Running i-PI over the network
-----------------------------

Understanding the network layout
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

.. _fig-network:

.. figure:: ../figures/ipi-network.*
   :class: white-background
   :width: 90.0%

   A schematic representation of the network layout one
   typically finds when running i-PI and the clients on a HPC system
   and/or on a local workstation.

The figure represents schematically a typical network
layout for a HPC system and a local workstation. When running i-PI
locally on a workstation, one can use the loopback interface (that can
be referred to as “localhost” in the “address” field of both i-PI and
the client) for communication. When running both i-PI and the clients on
a HPC cluster, one should work out which of the the interfaces that are
available on the node where the i-PI server runs are accessible from the
compute nodes. This requires some trial and error, and possibly setting
the “address” field dynamically from the job that launches i-PI. For
instance, if one was running i-PI on the login node, and the clients on
different compute nodes, as in panel b of the :ref:`i-PI running figure <fig-running>`, then on
the HPC system described in this scheme one should set
the address to that of the *ib1* interface – :math:`111.111.111.111` in
the example above. If instead i-PI was launched in a job script, then
the submission script would have to check for the IP address associated
with the *ib0* interface on the node the job has been dispatched to, and
set that address (e.g. :math:`111.111.111.200`) in the inputs of both
i-PI and the clients that will be launched in the same (or separate)
jobs.

Running i-PI on a separate workstation (panel c of :ref:`this figure <fig-running>`)
gives maximum flexibility, but is
also trickier as one has to reach the internet from the compute nodes,
that are typically not directly connected to it. We discuss this more
advanced setup in the next paragraph.

.. _ssh_sockets:

ssh tunnelling
~~~~~~~~~~~~~~

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
network layout of :ref:`this figure <fig-network>`, and if the i-PI server
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

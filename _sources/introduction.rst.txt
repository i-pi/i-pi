Program Overview
================

i-PI is a force engine that operates with a client-server
paradigm. i-PI acts as a server in command of the evolution of the nuclear positions, while one or more instances
of the client code handle the evaluation of system-specific properties. Designed to be universal, i-PI’s architecture
is agnostic to the identity of the force providers, with a general and flexible backend that accommodates to a wide range of client codes.

The code is written in Python 3, a high-level, general-purpose interpreted programming language known for
its emphasis on code readability. This choice facilitates rapid prototyping of new ideas and relatively easy code
maintenance when compared to compiled languages traditionally used in scientific computing.

i-PI is structured in a modular way that represents the underlying physics of the target simulation as faithfully
as possible. To achieve this, the code is organized around the System class which encodes all the information related 
to the physical system, such as the number and identity of atoms, the ensemble to be sampled, the number 
of replicas in path integral simulations, the initialization procedure, and the algorithm for evolving the 
nuclear positions. Forces are managed through the Forces class, which integrates the individual force components,
each of which is computed following the strategy specified in a Forcefield class. This two-layer approach is particu-
larly advantageous for algorithms where the total forces and energies are obtained by a combination of these
quantities computed at different accuracy levels, or among different system portions.
Simulation techniques that require the evolution of many systems simultaneously make use of the SMotion (Systems Motion)
class. This class enables the definition of evolution steps that combine the (possibly many) different systems, facil-
itating, for example, exchanges between system replicas as done in replica exchange simulations. Finally, 
estimators can be defined within the code or computed by the provided post-processing tools or user-generated scripts.

A schematic representation of the code structure and the server-client communication is presented in Fig. 1. 

.. figure:: ../figures/ipi-structure-v3.*
   :width: 90.0%
   :class: white-background

   Figure 1. Schematic representation of the i-PI code structure and the server-client communication.
  

In Figure 1, the physical system is defined by one or more replicas (beads), and the sampling conditions (e.g. temperature and pressure) by an Ensemble class. The forces acting on the atoms are constructed based on a combination of energy contributions evaluated by one or more instances of the 
Forcefield class (three in this example). Each Forcefield instance communicates with a different client sending positions (*q*) and 
lattice vectors (*h*), and receiving energy (*V*), forces (*f*), stresses and possible extra properties (X) as a JSON formatted string. 
In simulations with multiple replicas, each Forcefield instance can accept connections from several clients, to achieve parallel 
evaluation. The Motion class determines the evolution of the atoms (e.g. time integration, or geometry optimization), while 
“system motion” classes (SMotion) can act on multiple systems, e.g. to carry out replica exchange simulations. The output 
class handles writing the output and restart files



Communication protocol
~~~~~~~~~~~~~~~~~~~~~~

Since i-PI is designed to be used with a wide range of codes and
platforms, it has to rely on a simple and robust method for
communicating between the server and the client. Even though other choices
are possible, and it should be relatively simple to implement other
means of communication, the preferred approach relies on sockets as the
underlying infrastructure. Both Internet and Unix domain sockets can be
used: the latter allows for fast communication on a single node, whereas
the former makes it possible to realise a distributed computing paradigm,
with clients running on different nodes or even on different HPC
facilities. In order to facilitate the implementation of the socket
communication in client codes, a simple set of C wrappers to the
standard libraries socket implementation is provided as part of the i-PI
distribution, that can be used in any programming language that can be
linked with C code.

As far as the communication protocol is concerned, the guiding principle
has been keeping it to the lowest common denominator, and avoiding any
feature that may be code-specific. Only a minimal amount of information
is transferred between the client and the server; the position of the
atoms and cell parameters in one direction, and the forces, virial and
potential in the other.

For more details about sockets and communication, see
:ref:`distrib`.


Force evaluation
~~~~~~~~~~~~~~~~

Within i-PI, the evaluation of the forces plays a crucial role, as it is
the step requiring communication with the client code. In order to have
a flexible infrastructure that makes it possible to perform simulations
with advanced techniques, the force evaluation
machinery in i-PI might appear complicated at first, and deserves a
brief discussion.

.. figure:: ../figures/ipi-forces.*
   :class: white-background   
   :width: 90.0%

   Schematic representation of the different objects that
   are involved in the evaluation of the forces. The multiple layers and
   complex structure are necessary to give the possibility of
   decomposing the evaluation of the forces between multiple different
   clients and using different imaginary time partitioning (e.g. one can
   compute the bonded interactions using one client, and use a different
   client to compute the long-range electrostatic interactions,
   contracted on a single bead :cite:`mark-mano08jcp`).


This figure provides an overall scheme of the objects involved in the calculation 
of the forces. The infrastructure comprises
a force provider class that deals with the actual subdivision of work
among the clients, and a sequence of objects that translate the request
of the overall force of the system into atomic evaluations of one
component of the force 
When running path integral simulations, the latter refers to the component of an individual bead: 
i-PI is built to hide the path integral infrastructure from the client, and so beads must be
transferred individually.

Let us discuss for clarity a practical example: a calculation of an
empirical water model where the bonded interactions are computed on 32
beads by the program A, and the non-bonded interactions are computed by
client B, ring-polymer contracted on 8 beads. Each client “type” is
associated with a :ref:`forcefield` object in the input. In the case of a
:ref:`ffsocket` interface, the
forcefield object specifies the address to which a client should
connect, and so multiple clients of type A or B can connect to i-PI at
the same time. Each forcefield object deals with queueing force
evaluation requests and computing them in a first-in-first-out fashion,
possibly executing multiple requests in parallel.

On the force evaluation side, the task of splitting the request of a
force evaluation into individual components and individual beads is
accomplished by a chain of three objects, Forces, ForceComponent and
ForceBead. is the main force Forces evaluator, that is built from the
prototypes listed within the :ref:`forces` field of the :ref:`system`. 
Each `forcecomponent`  item within the :ref:`forces` tag
describe one component of the force – in our example one ForceComponent
bound to a forcefield of type A, evaluated on 32 beads, and one
ForceComponent bound to type B, evaluated on 8 beads. Forces contains
the machinery that automatically contracts the actual ring polymer to
the number of beads required by each component, and combines the various
components with the given weights to provide the overall force, energy
and virial where required. Note that in order to support ring polymer
contraction (RPC), the RPC procedure is executed even if no contraction
was required (i.e. even if all clients contain the full amount of
beads). ForceComponent is a very simple helper class that associates
with each bead a ForceBead object, which  is the entity in charge of
filing force requests to the appropriate ForceField object and waiting
for completion of the evaluation.


Automated evaluation (depend objects)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

i-PI uses a caching mechanism with automatic value updating to make the
code used to propagate the dynamics as simply and clearly  as possible.
Every physical quantity that is referenced in the code is created using
a “depend” object class, which is given the parameters on which it
depends and a function used to calculate its value.

.. figure:: ../figures/ipi-depend.*
   :class: white-background   
   :width: 90.0%

   Schematic overview of the functioning of the
   *depend* class used as the base for properties and physical
   quantities in i-PI. A few “primitive” quantities – such as atomic
   positions or momenta – can be modified directly. For most properties,
   one defines a function that can compute that property based on the
   value of other properties. Whenever one property is modified, all the
   quantities that depend on it are marked as tainted, so that when the
   value of one of the properties is used, the function can be invoked
   and the updated value obtained. If a quantity is not marked as
   tainted, the cached value is returned instead.

“Depend” objects can be called to get the physical quantity they
represent. However, they have further functionality. Firstly, once the
value of a “depend” object has been calculated, its value is cached, so
further references to that quantity will not need to evaluate the
function that calculates it. Furthermore, the code keeps track of when
any of the dependencies of the variable are updated, and makes sure that
the quantity is automatically when it is needed (i.e., when
the quantity is assessed again).

This is a minimal example of how to implement dependencies in a class

.. code-block:: python

   from ipi.utils.depend import depend_value, dproperties

   class DObject:
      def __init__(self):
         # depend objects are created using an underscore prefix. 
         # this is a "primitive" value that doesn't depend on anything
         self._scalar = depend_value(value=1, name="scalar")

         # this is a dependent object. the definition contains a function that
         # is called to determine the value, and specification of what objects
         # it depend on. 
         self._double = depend_value(func=lambda: 2*self.scalar, name="double", 
                                     dependencies=[self._scalar])

   # after the definition of a class, this helper function should be called to
   # create property getters and setters (use the names with no leading underscore)
   # note that property accessors are added to the class, not to the instances
   dproperties(DObject, ["scalar", "double"])

   myobj = DObject()
   # "primitive values" can be set manually
   myobj.scalar = 4
   # dependent objects will be computed automatically on demand
   print(myobj.double) # prints `8`

This choice makes implementation slightly more complex when the physical
observables are first introduced as variables, as one has to take care
of stating their dependencies as well as the function that computes
them. However, the advantage is that when the physical quantities are
used, in the integrator of the dynamics or in the evaluation of physical
properties, one does not need to take care of book-keeping and the code
can be clean, transparent and readable.

It is also possible to define dependencies between different objects, in 
which case it's necessary to make sure that the compute function has access, 
at runtime, to the value of the other object. A typical usage pattern is

.. code-block:: python

   # NB: this is meant to be run after the previous code snippet
   from ipi.utils.depend import depend_array
   import numpy as np

   class DOther:
      def __init__(self):
         
         # depend arrays must be initialized with storage space
         self._vec = depend_array(value=np.ones(4), name="vec") 

      def bind(self, factor):

         self.factor = factor # stores a reference to the object holding the value
         self._scaled = depend_array(value=np.ones(4), name="vec",
                                    func=self.get_scaled, dependencies=[self._vec])

         # dependencies (or dependants) can also be added after object creation
         self._scaled.add_dependency(self.factor._double)

      def get_scaled(self):
         # computes a scaled version of the vector
         return self.vec*self.factor.double

   dproperties(DOther, ["vec", "scaled"])

   myoth = DOther() # creates the object
   myoth.bind(myobj) # makes connections
   
   myoth.vec = np.asarray([0,1,2,3])
   print(myoth.scaled) # prints [0,8,16,24]

   myoth.vec[3] = 0 # depend_arrays can be accessed as normal np.ndarray
   print(myoth.scaled) # prints [0,8,16,0]


Licence and credits
-------------------

The code is distributed under the dual GPL and MIT licence. For more details
see `www.gnu.org/licences/gpl.html <www.gnu.org/licences/gpl.html>`__ and 
https://fedoraproject.org/wiki/Licensing:MIT.

If you use this code in any future publications, please cite this using
:cite:`ceri+14cpc` for v1,
:cite:`Kapil:2019ju` for v2.
:cite:`litman2024ipi` for v3.

Contributors
~~~~~~~~~~~~

i-PI was originally written by M. Ceriotti and J. More at Oxford
University, together with D. Manolopoulos. The updated list of developers and 
contributors can be found 
`here <https://ipi-code.org/about/developers/>`__



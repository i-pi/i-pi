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

Internet and Unix domain sockets are the most general transports, but
they are not the only ones. On a single node a *shared-memory* socket
(``mode="shm"``) can move the bulk position/force arrays through POSIX
shared memory; many cheap evaluations can be *batched* into a single
exchange (``batch_size``); the ``ffmpi`` forcefield communicates over MPI
for clients co-launched with i-PI in one MPMD job; and ``ffdirect``
evaluates a Python potential in-process, with no inter-process
communication at all.

For more details about all of these, and how to choose between them, see
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

“Depend” objects store the physical quantity they represent, but they
have further functionality. Firstly, once the value of a “depend” object
has been calculated, its value is cached, so further references to that
quantity will not need to evaluate the function that calculates it.
Furthermore, the code keeps track of when any of the dependencies of the
variable are updated, and makes sure that the quantity is recomputed
automatically when it is needed (i.e., when the quantity is accessed again).

There are two classes of depend objects: ``depend_value`` holds a generic
Python object (typically a scalar), and ``depend_array`` holds an array.
Both can be "primitive", i.e. set directly, or computed by a function from
the values of other depend objects, which are declared as their
dependencies. Updates are lazy: modifying a value immediately marks all the
quantities that depend on it, directly or indirectly, as "tainted", but a
tainted quantity is only recomputed when it is read.

This is a minimal example of how to implement dependencies in a class

.. code-block:: python

   from ipi.utils.depend import depend_value, dproperties

   class DObject:
      def __init__(self):
         # depend objects are stored in attributes with a trailing underscore.
         # this is a "primitive" value that doesn't depend on anything
         self.scalar_ = depend_value(value=1, name="scalar")

         # this is a dependent object. the definition contains a function that
         # is called to determine the value, and specification of what objects
         # it depend on. 
         self.double_ = depend_value(func=lambda: 2*self.scalar, name="double", 
                                     dependencies=[self.scalar_])

   # after the definition of a class, this helper function should be called to
   # create property getters and setters (use the names with no trailing underscore)
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
   from ipi.utils.depend import depend_array, ddeps
   import numpy as np

   class DOther:
      def __init__(self):
         
         # depend arrays must be initialized with storage space
         self.vec_ = depend_array(value=np.ones(4), name="vec") 

      def bind(self, factor):

         self.factor = factor # stores a reference to the object holding the value
         self.scaled_ = depend_array(value=np.ones(4), name="vec",
                                    func=self.get_scaled, dependencies=[self.vec_])

         # dependencies can also be added after object creation
         ddeps(self.scaled_, [self.factor.double_])

      def get_scaled(self):
         # computes a scaled version of the vector
         return self.vec*self.factor.double

   dproperties(DOther, ["vec", "scaled"])

   myoth = DOther() # creates the object
   myoth.bind(myobj) # makes connections
   
   myoth.vec = np.asarray([0,1,2,3])
   print(myoth.scaled) # prints [0,8,16,24]

   # reading `myoth.vec` returns a read-only array. parts of a depend_array
   # are modified through the depend object, which updates the dependencies
   myoth.vec_[3] = 0
   print(myoth.scaled) # prints [0,8,16,0]
   myoth.vec_ += 1 # in-place operations are also done on the depend object
   print(myoth.scaled) # prints [8,16,24,8]

Indexing a depend object, e.g. ``myoth.vec_[1:3]``, returns a view of
part of the array that remains connected to the dependency graph: it
shares the memory of the full array, writing through it updates the
dependants of the full array, and reading a stale view recomputes the
full array. Views are used, for instance, to expose the positions of each
replica of a ring polymer as a separate ``Atoms`` object, or to bind a
thermostat to the momenta of one normal mode.

Rules of thumb
^^^^^^^^^^^^^^

For a depend member stored as ``obj.x_`` and exposed with ``dproperties``
as ``obj.x``:

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Operation
     - Syntax
   * - read the value (recomputed if needed)
     - ``obj.x``
   * - overwrite the whole value
     - ``obj.x = v``
   * - modify part of an array
     - ``obj.x_[i] = v``
   * - in-place update of a whole array
     - ``obj.x_ += v`` (faster than ``obj.x_[:] += v``)
   * - declare a dependency
     - ``dependencies=[other.y_]``, or ``ddeps(obj.x_, [other.y_])`` after
       creating a computed object
   * - read the value of a depend object that is not a class member
     - ``dep.get()``
   * - keep a view of part of an array that stays connected to the
       dependency graph (e.g. to bind a thermostat to some of the momenta)
     - ``obj.x_[i]``, ``obj.x_.flatten()``
   * - make an existing depend object computed from other ones
     - ``dcompute(obj.x_, func, [other.y_])``, ``dpipe(other.y_, obj.x_)``

Common pitfalls:

* Arrays read as ``obj.x`` are read-only, so ``obj.x[i] = v`` and
  ``obj.x += v`` raise an error: modify the array through ``obj.x_``.
* A value read as ``obj.x`` is only valid until ``x`` is written or
  recomputed: depending on the array backend, these may or may not change
  it. To keep track of ``x``, store the depend object ``obj.x_`` (or a view
  of it) instead.
* A ``depend_value`` that holds an array is not protected: modifying the
  array in place does not update its dependants.
* Depend objects are not arrays: use ``obj.x`` (or ``obj.x_.get()``) when
  the value is needed.
* Functions that compute a value should not modify other depend objects:
  they are called at unpredictable times, whenever a stale value is read.

Each read of ``obj.x`` checks whether the value is stale, and recomputes
it if needed. Recomputation is protected by a lock, so values can be read
from several threads, e.g. when the outputs of different systems are
written in parallel. Reading an up-to-date value is cheap (a fraction of a
microsecond), but in loops that do not modify it, it is still better to
read the value once into a local variable.

``dpipe`` makes a depend object follow the value of another one, e.g. to
share the temperature between the ensemble and the thermostats:

.. code-block:: python

   from ipi.utils.depend import dpipe

   class Thermo:
      def __init__(self):
         self.temp_ = depend_value(name="temp", value=0.0)

   dproperties(Thermo, ["temp"])

   thermo = Thermo()
   dpipe(myobj.double_, thermo.temp_) # thermo.temp follows myobj.double
   myobj.scalar = 5
   print(thermo.temp) # prints `10`

More generally, ``dcompute(dep, func, dependencies)`` makes an existing
depend object computed, as if it had been created with ``func`` and
``dependencies``, which is useful when the function is only known when
objects are bound together. Its current value is discarded, and it can no
longer be written. ``dpipe(dfrom, dto)`` is the same as
``dcompute(dto, dfrom.get, [dfrom])``, except that it can also make an object
that already follows another one follow ``dfrom`` instead, and
``dpipe(dfrom, dto, item=i)`` makes ``dto`` follow ``dfrom[i]``.

``dsynchro`` links two representations of the same data, that can both be
modified, such as the bead and normal-mode coordinates of a ring polymer.
It takes the two depend objects and the two functions that compute each
one from the other. Writing one of the two marks the other as stale, and
reading the stale one recomputes it from the one that was written last.
When they are linked, the first object is taken to hold the current data.
A view of a synchronized array (e.g. the positions of a single replica)
is synchronized as well.

.. code-block:: python

   from ipi.utils.depend import dsynchro

   class Length:
      def __init__(self):
         self.meters_ = depend_array(name="meters", value=np.zeros(3))
         self.feet_ = depend_array(name="feet", value=np.zeros(3))
         dsynchro(self.meters_, self.feet_,
                  lambda: self.feet * 0.3048,    # meters from feet
                  lambda: self.meters / 0.3048)  # feet from meters

   dproperties(Length, ["meters", "feet"])

   length = Length()
   length.meters = np.ones(3)
   print(length.feet) # prints [3.2808399 3.2808399 3.2808399]
   length.feet_[0] = 0.0
   print(length.meters) # prints [0. 1. 1.]

Array backends
^^^^^^^^^^^^^^

The arrays held by depend objects are created with the array backend
selected at launch (see :ref:`array-backends`), so the code that manipulates
them must work with numpy arrays as well as with torch tensors, possibly on a
GPU. The namespace of the current backend is ``xp``
(``from ipi.utils.array_backend import xp``): it is the numpy module with the
numpy backend, and a wrapper of torch that follows the
`array API standard <https://data-apis.org/array-api/>`_ (from
``array_api_compat``) with the torch backend. ``xp`` is updated in place when
the backend is selected, so it can be imported at module level, but arrays
must not be created at import time (module constants, default arguments).

* Arrays that are part of the simulation state, and the results of
  computations on them, use ``xp`` functions and operators (``a @ b`` rather
  than ``np.dot(a, b)``). The array methods ``.reshape``, ``.sum``,
  ``.mean`` and ``.T`` are portable too, and faster than the corresponding
  numpy functions. Not portable: ``.copy()`` (use
  ``xp.asarray(x, copy=True)``), ``.astype``, ``.flatten``, ``.size``,
  ``.max(axis=...)``, and slices with negative steps (use ``xp.flip``).
* Python scalars (parameters, time steps, temperatures) stay Python floats,
  and functions applied to them use ``math`` or ``np``: ``xp.sqrt(2.0)``
  fails with torch.
* The creation functions of ``xp`` produce double-precision arrays on the
  selected device: do not pass ``dtype=`` or ``device=`` for floating-point
  arrays.
* Data enters and leaves the backend only at the boundaries. Input values
  and files are converted with ``xp.asarray`` before being stored; outputs,
  checkpoints, sockets, forcefields and calls to scipy use ``to_numpy``,
  which may return a view of the stored data with numpy (copy it if it is
  kept). Depend objects never convert their argument: writing a numpy array
  into a torch-backed ``depend_array`` raises an error.
* Small arrays used for indexing or control on the host (atom names, index
  lists, the number of multiple-time-step sub-steps, the time scales of the
  normal-mode thermostats) are stored as numpy arrays with
  ``depend_array(..., host=True)``.
* With a GPU, ``float(x)``, ``bool(x)``, ``if x > 0`` and ``x.item()`` on a
  device array wait for all the queued operations to complete: avoid them in
  the MD step, and prefer array operations (e.g. ``xp.where``).
* ``prng.gvec`` and ``prng.gfill`` produce backend arrays (drawn on the device
  with torch), while scalar draws (``prng.u``, ``prng.g``, ``prng.gamma``, ...)
  are Python or numpy numbers.

The tests in ``ipi_tests/unit_tests/engine/test_backend_equivalence.py``
check that the torch backend reproduces the numpy trajectories and that an MD
step does not synchronize with the GPU. The regression tests and the examples
can be run with another backend with ``pytest --ipi-backend torch
[--ipi-device cuda]``. In these runs (and in the tests that use the
``strict_torch`` fixture) implicit conversions of torch tensors to numpy
arrays raise an error also on the CPU, as they do on GPUs, so that missing
``to_numpy`` calls are detected without a GPU.


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



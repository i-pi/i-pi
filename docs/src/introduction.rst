Introduction
============

Path Integral Molecular Dynamics
--------------------------------

Molecular dynamics (MD) is a technique used to study the properties of a
system of interacting particles by applying Newton’s equations of motion
to produce trajectories which can be used to efficiently explore the
phase space. This can be used to calculate many equilibrium and
dynamical properties and to study systems from isolated gas molecules to
condensed phase bulk materials.
 
However, while this technique has been very successful, in most MD
implementations the assumption is made that the nuclei behave as
classical particles, which for light nuclei such as hydrogen is often a
very poor approximation as the effect of zero-point energy (ZPE) and
quantum tunnelling can be large. For example, even at room temperature
the vibrational frequency of an OH stretch in water is over 15 times
larger than the available thermal energy, and so this motion will be
highly quantized. The current state-of-the-art method to include nuclear
quantum effects (NQE) in the calculation of static properties of
condensed phase systems is path integral molecular dynamics (PIMD).

PIMD generates the quantum-mechanical ensemble of a system of
interacting particles by using MD in an extended phase space. This is
derived from the path integral formalism
:cite:`feyn-hibb65book`, which relates the statistics of a
collection of quantum particles to those of a set of classical ring
polymers, a ring polymer being a number of replicas of a particle
coupled by harmonic springs. This so-called classical isomorphism is
exact in the limit as the number of replicas goes to infinity, but in
practice is converged numerically with only a finite number.

This then allows quantum phase space averages to be calculated from
classical trajectories, with only about an order of magnitude more
computing time than would be required for standard MD. Also, since PIMD
is simply classical MD in an extended phase space, many of the
techniques developed to improve the scope and efficiency of MD
simulations can be applied straightforwardly to the equivalent PIMD
calculations :cite:`ceri+10jcp,mart+99jcp`. Finally, several
techniques designed specifically for PIMD simulations are now available
to increase the rate of convergence with respect to the number of
replicas used
:cite:`mark-mano08jcp,ceri+11jcp,suzu95pla,chin97pla,ceri+12prsa,pere-tuck11jcp`,
further reducing the computational overhead of the method. All of these
facts mean that it is now feasible to do PIMD simulations with thousands
of molecules, or even to use *ab initio* electronic structure
calculations to propagate the dynamics for small systems.

Furthermore, the framework used to run PIMD simulations can be adapted
to generate approximate quantum dynamical information
:cite:`cao-voth93jcp,cao-voth94jcp,crai-mano04jcp,braa-mano06jcp`,
and so can also be used to calculate correlation functions. While
real-time quantum coherences cannot be captured, the inclusion of
quantum statistical information and the rapid decoherence observed in
condensed phase systems mean that in many cases very accurate results
can be obtained from such approximate treatments of quantum dynamics
:cite:`habe+13arpc`.

Implementation
--------------

Automated evaluation (depend objects)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

i-PI uses a caching mechanism with automatic value updating to make the
code used to propagate the dynamics as simple and clear as possible.
Every physical quantity that is referenced in the code is created using
a “depend” object class, which is given the parameters on which it
depends and a function used to calculate its value.

.. figure:: ../figures/ipi-depend.*
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

Force evaluation
~~~~~~~~~~~~~~~~

Within i-PI, the evaluation of the forces plays a crucial role, as it is
the step requiring communication with the client code. In order to have
a flexible infrastructure that makes it possible to perform simulations
with advanced techniques such as ring-polymer
contraction :cite:`mark-mano08jcp`, the force evaluation
machinery in i-PI might appear complicated at first, and deserves a
brief discussion.

.. figure:: ../figures/ipi-forces.*
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
of the forces. The infrastracture comprises
a force provider class that deals with the actual subdivision of work
among the clients, and a sequence of objects that translate the request
of the overall force of the system into atomic evaluations of one
component of the force for an individual bead: i-PI is built to hide the
path integral infrastructure from the client, and so beads must be
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
Each :ref:`forcecomponent`  item within the :ref:`forces` tag
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
with each bead a ForceBead object, that is the entity in charge of
filing force requests to the appropriate ForceField object and waiting
for completion of the evaluation.

Communication protocol
~~~~~~~~~~~~~~~~~~~~~~

Since i-PI is designed to be used with a wide range of codes and
platforms, it has to rely on a simple and robust method for
communicating between the server and client. Even though other choices
are possible, and it should be relatively simple to implement other
means of communication, the preferred approach relies on sockets as the
underlying infrastructure. Both Internet and Unix domain sockets can be
used: the latter allow for fast communication on a single node, whereas
the former make it possible to realise a distributed computing paradigm,
with clients running on different nodes or even on different HPC
facilities. In order to facilitate implementation of the socket
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

Licence and credits
-------------------

Most of this code is distributed under the GPL licence. For more details
see `www.gnu.org/licences/gpl.html <www.gnu.org/licences/gpl.html>`__.
So that they can easily be incorporated in other codes, the files in the
directory “drivers” are all held under the MIT licence. For more details
see https://fedoraproject.org/wiki/Licensing:MIT.

If you use this code in any future publications, please cite this using
:cite:`ceri+14cpc` for v1 and
:cite:`Kapil:2019ju` for v2.

Contributors
~~~~~~~~~~~~

i-PI was originally written by M. Ceriotti and J. More at Oxford
University, together with D. Manolopoulos. Several people contributed to
its further development. Developers who implemented a specific feature
are acknowledged above.

On-line resources
-----------------

Python resources
~~~~~~~~~~~~~~~~

For help with Python programming, see
`www.python.org <www.python.org>`__. For information about the NumPy
mathematical library, see `www.numpy.org <www.numpy.org>`__, and for
worked examples of its capabilities see
`www.scipy.org/Tentative_NumPy_Tutorial <www.scipy.org/Tentative_NumPy_Tutorial>`__.
Finally, see http://hgomersall.github.io/pyFFTW/ for documentation on
the Python FFTW library that is currently implemented with i-PI.

.. _librarywebsites:

Client code resources
~~~~~~~~~~~~~~~~~~~~~

Several codes provide out-of-the-box an i-PI interface, including CP2K,
DFTB+, Lammps, Quantum ESPRESSO, Siesta, FHI-aims, Yaff, deMonNano, TBE.
If you are interested in interfacing your code to i-PI please get in
touch, we are always glad to help!

There are several Fortran and C libraries that most client codes will
probably need to run, such as FFTW, BLAS and LAPACK. These can be found
at `www.fftw.org <www.fftw.org>`__,
`www.netlib.org/blas <www.netlib.org/blas>`__ and
`www.netlib.org/lapack <www.netlib.org/lapack>`__ respectively.

These codes do not come as part of the i-PI package, and must be
downloaded separately. See chapter :ref:`clientinstall` for more
details of how to do this.

i-PI resources
~~~~~~~~~~~~~~

For more information about i-PI and to download the source code go to
http://ipi-code.org/.

In http://gle4md.org/ one can also obtain colored-noise parameters to
run Path Integral with Generalized Langevin Equation thermostat
(PI+GLE/PIGLET) calculations.

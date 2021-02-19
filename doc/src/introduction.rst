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
the quantity is automatically recomputed when it is needed (i.e., when
the quantity is assessed again).

This choice makes implementation slightly more complex when the physical
observables are first introduced as variables, as one has to take care
of stating their dependencies as well as the function that computes
them. However, the advantage is that when the physical quantities are
used, in the integrator of the dynamics or in the evaluation of physical
properties, one does not need to take care of book-keeping and the code
can be clean, transparent and readable.

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

Force evaluation
~~~~~~~~~~~~~~~~

Within i-PI, the evaluation of the forces plays a crucial role, as it is
the step requiring communication with the client code. In order to have
a flexible infrastructure that makes it possible to perform simulations
with advanced techniques such as ring-polymer
contraction :cite:`mark-mano08jcp`, the force evaluation
machinery in i-PI might appear complicated at first, and deserves a
brief discussion.

A scheme of the objects involved in the calculation of the forces is
presented in Figure `1.3 <#fig:forces>`__. The infrastracture comprises
a force provider class that deals with the actual subdivision of work
among the clients, and a sequence of objects that translate the request
of the overall force of the system into atomic evaluations of one
component of the force for an individual bead: i-PI is built to hide the
path integral infrastructure from the client, and so beads must be
transferred individually.

Let us discuss for clarity a practical example – a calculation of an
empirical water model where the bonded interactions are computed on 32
beads by the program A, and the non-bonded interactions are computed by
client B, ring-polymer contracted on 8 beads. Each client “type” is
associated with a object in the input. In the case of a interface, the
forcefield object specifies the address to which a client should
connect, and so multiple clients of type A or B can connect to i-PI at
the same time. Each forcefield object deals with queueing force
evaluation requests and computing them in a first-in-first-out fashion,
possibly executing multiple requests in parallel.

On the force evaluation side, the task of splitting the request of a
force evaluation into individual components and individual beads is
accomplished by a chain of three objects, Forces, ForceComponent and
ForceBead. is the main force Forces evaluator, that is built from the
prototypes listed within the field of . Each item within the tag
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
`3.3 <#distrib>`__.

.. _units:

Internal units
~~~~~~~~~~~~~~

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
   atomic units, or angstrom for PDB files.

-  If units are specified in both input and structure file and they
   match, conversion happens just once. If they do not match, an error
   is raised and i-PI stops.

Core features
-------------

i-PI includes a large number of advanced molecular dynamics features,
with an obvious focus on path integral molecular dynamics, but also
several methods for sampling classical trajectories.

Features in version 1.0
~~~~~~~~~~~~~~~~~~~~~~~

-  molecular dynamics and PIMD in the *NVE*, *NVT* and *NPT* ensembles,
   with the high-frequency internal vibrations of the path propagated in
   the normal-mode representation :cite:`ceri+10jcp` to
   allow for longer time steps;

-  ring polymer
   contraction :cite:`mark-mano08jcp,mark-mano08cpl`,
   implemented by exposing multiple socket interfaces to deal with short
   and long-range components of the potential energy separately;
   treating different components that have different computational cost
   and characteristic time scale separately can reduce substantially the
   overall effort associated with a simulation

-  efficient stochastic velocity rescaling
   :cite:`buss+07jcp` and path integral Langevin equation
   thermostats :cite:`ceri+10jcp` to sample efficiently the
   canonical ensemble;

-  various generalized Langevin equation (GLE) thermostats, including
   the optimal sampling :cite:`ceri+09prl,ceri+10jctc`,
   quantum  :cite:`ceri+09prl2`, and
   :math:`\delta` :cite:`ceri-parr10pcs` thermostats; the
   parameters for different GLE flavors and the conditions in which they
   should be applied can be obtained from a separate
   website :cite:`gle4md`;

-  mixed path integral–generalized Langevin equation techniques for
   accelerated convergence, including both
   PI+GLE :cite:`ceri+11jcp` and the more recent and
   effective version PIGLET :cite:`ceri-mano12prl`; these
   techniques reduce the number of path integral replicas needed, while
   allowing for systematic convergence;

-  all the standard estimators for structural properties, the quantum
   kinetic energy, pressure, etc.;

-  more sophisticated estimators such as the scaled-coordinate heat
   capacity estimator :cite:`yama05jcp`, estimators to
   obtain isotope fractionation free energies by re-weighting a
   simulation of the most abundant
   isotope :cite:`ceri-mark13jcp`, and a displaced-path
   estimator for the particle momentum
   distribution :cite:`lin+10prl`;

-  the infrastructure that is needed to perform approximate quantum
   dynamics calculations such as ring polymer molecular dynamics
   (RPMD) :cite:`crai-mano04jcp,habe+13arpc` and centroid
   molecular dynamics
   (CMD) :cite:`cao-voth93jcp,cao-voth94jcp`.

Features added in version 2.0
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Further details can be found in Ref. :cite:`Kapil:2019ju`.

-  reweighted fourth-order path integral MD (M. Ceriotti, G.A.R.
   Brain) :cite:`ceri+12prsa,jang-voth01jcp`; this method
   makes it possible to obtain fourth-order statistics by re-weighting
   second-order trajectories; attention should be paid to avoid
   statistical inefficiencies;

-  finite-differences implementation of fourth-order path integrals (V.
   Kapil, M. Ceriotti) :cite:`kapi+16jcp2`; this schemes
   enables explicit fourth-order path integral simulations, that
   converge faster than conventional Trotter methods;

-  perturbed path integrals (I.
   Poltavsky) :cite:`polt-tkat16cs`; essentially, a
   truncated cumulant expansion of fourth-order reweighting, that often
   enables fast convergence avoiding statistical instability;

-  open path integrals and momentum distribution estimators (V. Kapil,
   A. Cuzzocrea, M. Ceriotti) :cite:`kapi+18jpcb`; makes it
   possible to compute the particle momentum distribution including
   quantum fluctuations of nuclei;

-  quantum alchemical transformations (B. Cheng)
   :cite:`liu+13jpcc,chen+16jpcl`; Monte Carlo exchanges
   between isotopes of different mass, useful to sample isotope
   propensity for different molecules or environments;

-  direct isotope fractionation estimators (B. Cheng, M. Ceriotti)
    :cite:`chen-ceri14jcp`; avoid thermodynamic integration
   to obtain isotope fractionation ratios;

-  spatially localized ring polymer contraction (M. Rossi, M. Ceriotti)
   :cite:`litm+17jcp`; simple contraction scheme for weakly
   bound molecules, e.g. on surfaces;

-  ring polymer instantons (Y. Litman, J.O. Richardson, M. Rossi);
   evaluation of reaction rates and tunnelling splittings for molecular
   rearrangements and chemical reactions;

-  thermodynamic integration (M. Rossi, M. Ceriotti)
   :cite:`ross+16prl`; classical scheme to compute free
   energy differences;

-  geometry optimizers for minimization and saddle point search (B.
   Helfrecht, R. Petraglia, Y. Litman, M. Rossi)
   :cite:`ross+16prl`

-  harmonic vibrations through finite differences (V. Kapil, S.
   Bienvenue) :cite:`ross+16prl`; simple evaluation of the
   harmonic Hessian;

-  multiple time stepping (V. Kapil, M. Ceriotti)
   :cite:`kapi+16jcp`; accelerated simulations by separating
   slow and fast degrees of freedom into different components of the
   potential energy;

-  metadynamics through a PLUMED interface (G. Tribello, M. Ceriotti);
   simulation of rare events and free energy calculations;

-  replica exchange MD (R. Petraglia, R. Meissner, M.
   Ceriotti) :cite:`petr+15jcc`; accelerated convergence of
   averages by performing Monte Carlo exchanges of configurations
   between parallel calculations

-  thermostatted RPMD :cite:`ross+14jcp`, including
   optimized-GLE TRPMD :cite:`ross+18jcp`; reduces
   well-known artifacts in the simulation of dynamical properties by
   path integral methods;

-  dynamical corrections to Langevin trajectories (M. Rossi, V. Kapil,
   M. Ceriotti) :cite:`ross+18jcp`; eliminates the artifacts
   introduced into dynamical properties by the presence of thermostats;

-  fast forward Langevin thermostat (M. Hijazi, D. M. Wilkins, M.
   Ceriotti); a simple scheme to reduce the impact of strongly-damped
   Langevin thermostats on sampling efficiency;
   :cite:`hija+18jcp`

-  Langevin sampling for noisy and/or dissipative forces (J. Kessler, T.
   D. Kühne); suitable to stabilize and correct the artifacts that are
   introduced in MD trajectories by different extrapolation schemes;

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
downloaded separately. See chapter `2.2 <#clientinstall>`__ for more
details of how to do this.

i-PI resources
~~~~~~~~~~~~~~

For more information about i-PI and to download the source code go to
http://ipi-code.org/.

In http://gle4md.org/ one can also obtain colored-noise parameters to
run Path Integral with Generalized Langevin Equation thermostat
(PI+GLE/PIGLET) calculations.

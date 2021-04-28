Getting started
===============

.. _install:

Installing i-PI
---------------

Requirements
~~~~~~~~~~~~

i-PI is Python code, and as such strictly speaking does not need to be
compiled and installed. The ``i-pi`` file in the ``bin`` directory of
the distribution is the main (executable) script, and can be run as long
as the system has installed:

-  Python version 2.4 or greater (i-PI up to version v2.0 is not Python3
   compatible)

-  The Python numerical library NumPy

See `2.3.1 <#runningsimulations>`__ for more details on how to launch
i-PI.

Additionally, most client codes will have their own requirements. Many
of them, including the test client codes given in the “drivers”
directory, will need a suitable Fortran compiler. A C compiler is
required for the sockets.c wrapper to the sockets standard library. Most
electronic structure codes will also need to be linked with some
mathematical libraries, such as BLAS, FFTW and LAPACK. Installation
instructions for these codes should be provided as part of the code
distribution and on the appropriate website, as given in
`1.6.2 <#librarywebsites>`__. Patching for use with i-PI should not
introduce further dependencies.

Using the setup.py module
^^^^^^^^^^^^^^^^^^^^^^^^^

While the ``i-pi`` file can be used to run any i-PI simulation, it is
often more convenient to install the package to the system’s Python
modules path, so that it is accessible by all users and can be run
without specifying the path to the Python script.

For this purpose we have included a module in the root directory of the
i-PI distribution, ``setup.py``, which handles creating a package with
the executable and all the modules which are necessary for it to run.
The first step is to build the distribution using:

.. code-block::

   > python setup.py build

Note that this requires the distutils package that comes with the
python-dev package.

This creates a "build" directory containing only the files that are used
to run an i-PI simulation, which can then be used to create the
executable. This can be done in two ways, depending on whether or not
the user has root access. If the user does have root access, then the
following command will add the relevant source files to the standard
Python library directory:

.. code-block::

   > python setup.py install

This will install the package in the /usr/lib/py_vers directory, where
py_vers is the version of Python that is being used. This requires
administrator priviledges.

Otherwise, one can install i-PI in a local Python path. If such path
does not exist yet, one must create directories for the package to go
into, using:

Next, you must tell Python where to find this library, by appending to
the Linux environment variable PYTHONPATH, using:

Finally, the code can be installed using:

.. code-block::

   > python setup.py install –prefix= 

Either way, it will now be possible to run the code automatically, using

i-PI download
~~~~~~~~~~~~~

You can find information how to download i-PI from
http://ipi-code.org/download/.

The i-PI executable, found in the ``bin`` folder, will run immediately,
without needing to be installed or compiled, but we include a setup.py
module in the main directory so it can be installed to the Python tree
if so desired.

Installing NumPy
~~~~~~~~~~~~~~~~

NumPy is the standard Python mathematics library, and is used for most
of the array manipulation and linear algebra in i-PI. It should be
installed alongside most standard Python environments on HPC facilities.
Otherwise, it is generally relatively straightforward to install it.

In any case you must first obtain the NumPy code, which can be
downloaded as a tar file from http://www.numpy.org. If the version of
NumPy being installed is given by “np_vers”, this can be extracted
using:

Before installing this code it first needs to be configured correctly.
Note that this requires the distutils package that comes with the
python-dev package. Assuming that the required software is installed,
the NumPy package is built using:

.. code-block::

   > python setup.py build

The next step is to install NumPy. By default the download is to the
directory /usr/local. If you have root access, and so can write to /usr,
then all that needs to be done to finish the install is:

.. code-block::

   > python setup.py install

If you do not have root access, then the next step depends on which
version of Python is beind used. With versions 2.6 or later there is a
simple command to automatically download into the directory $HOME/local:

.. code-block::

   > python setup.py install –user

With Python 2.4/2.5 the process is a little more involved. First you
must explicitly install the package in the directory of choice, “np_dir”
say, with the following command:

Next, you must tell Python where to find this library, by appending to
the Linux environment variable PYTHONPATH. If you are using Python
version “py_vers”, then the NumPy libraries will have been installed in
the directory “np_dir/lib/py_vers/site-packages”, or a close analogue of
this. In the above case the following command will allow the Python
interpreter to find the NumPy libraries:

Now Python scripts can import the NumPy libraries using:

.. code-block::

   import numpy

PyFFTW
~~~~~~

Some of the steps in the dynamics algorithm involve a change of
variables from the bead coordinates to the normal modes of the ring
polymers. Currently, this transformation is, at least by default,
computed using a fast-Fourier transform (FFT) library within the NumPy
distribution. This however is not the only distribution that could be
used, and indeed faster stand-alone versions exist. The gold-standard
FFT library is the FFTW library, which is a set of C libraries that have
been heavily optimized for a wide range of applications. There have been
a number of Python wrappers built around the FFTW library, one of which
is currently interfaced with i-PI. This code can be found at
https://github.com/hgomersall/pyFFTW, and has documentation at
http://hgomersall.github.io/pyFFTW/.

This code has the following dependencies:

-  Python version 2.7 or greater

-  Numpy version 1.6 or greater

-  FFTW version 3.2 or greater

This can be installed in the same way as NumPy, except using the code
distribution above, or using various installation packages as per the
instructions on the above documentation. Note that no other options need
to be specified in the input file; i-PI will check to see if this
library is available, and if it is it will be used by default. Otherwise
the slower NumPy version will be used.

.. _clientinstall:

Installing clients
------------------

As of today, the following codes provide out-of-the-box an i-PI
interface: CP2K, DFTB+, Lammps, Quantum ESPRESSO, Siesta, FHI-aims,
Yaff, deMonNano, TBE. Links to the webpages of these codes, including
information on how to obtain them, can be found in http://ipi-code.org/.

If you are interested in interfacing your code to i-PI please get in
touch, we are always glad to help. We keep some information below in
case you are interested in writing a patch to a code.

Writing a patch
~~~~~~~~~~~~~~~

If you have edited a client code, and wish to make a patch available for
the new version, then this can be done very simply. If your edited code
is in a directory “new”, and a clean distribution is held in a directory
“old”, then a patch “changes.patch” can be created using:

.. code-block::

   > diff -rupN old/ new/ > changes.patch

Running i-PI
------------

i-PI functions based on a client-server protocol, where the evolution of
the nuclear dynamics is performed by the i-PI server, whereas the energy
and forces evaluation is delegated to one or more instances of an
external program, that acts as a client. This design principle has
several advantages, in particular the possibility of performing PIMD
based on the forces produced by one’s favourite electronic
structure/empirical force field code. However, it also makes running a
simulation slightly more complicated, since the two components must be
set up and started independently.

.. _runningsimulations:

Running the i-PI server
~~~~~~~~~~~~~~~~~~~~~~~

i-PI simulations are run using the i-pi Python script found in the
“i-pi” directory. This script takes an xml-formatted file as input, and
automatically starts a simulation as specified by the data held in it.
If the input file is called “input_file.xml”, then i-PI is run using:

This reads in the input data, initializes all the internally used
objects, and then creates the server socket. The code will then wait
until at least one client code has connected to the server before
running any dynamics. Note that until this has happened the code is
essentially idle, the only action that it performs is to periodically
poll for incoming connections.

.. _runningclients:

Running the client code
~~~~~~~~~~~~~~~~~~~~~~~

Below we give examples on how to make different clients communicate with
i-PI. Most clients also include descriptions on how to do this from
their own documentation.

.. _driver.x:

Built-in, example client
^^^^^^^^^^^^^^^^^^^^^^^^

While i-PI is designed with *ab initio* electronic structure
calculations in mind, it also includes a Fortran empirical potential
client code to do simple calculations and to run the examples.

The source code for this is included in the directory “drivers”, and can
be compiled into an executable “i-pi-driver” using the UNIX utility
make.

This code currently has four empirical potentials hardcoded into it, a
Lennard-Jones potential, the Silvera-Goldman potential
:cite:`silv-gold78jcp`, a 1D harmonic oscillator potential,
and the ideal gas (i.e. no potential interaction).

How the code is run is based on what command line arguments are passed
to it. The command line syntax is:

.. code-block::

   > i-pi-driver [-u] -h hostname -p port -m [gas|lj|sg|harm] -o
   parameters [-v]

The flags do the following:

-u:
   Optional parameter. If specified, the client will connect to a unix
   domain socket. If not, it will connect to an internet socket.

-h:
   Is followed in the command line argument list by the hostname of the
   server.

-p:
   Is followed in the command line argument list by the port number of
   the server.

-m:
   Is followed in the command line argument list by a string specifying
   the type of potential to be used. “gas” gives no potential, “lj”
   gives a Lennard-Jones potential, “sg” gives a Silvera-Goldman
   potential and “harm” gives a 1D harmonic oscillator potential. Other
   options should be clear from their description.

-o:
   Is followed in the command line argument list by a string of comma
   separated values needed to initialize the potential parameters. “gas”
   requires no parameters, “harm” requires a spring constant, “sg”
   requires a cut-off radius and “lj” requires the length and energy
   scales and a cut-off radius to be specified. All of these must be
   given in atomic units.

-v:
   Optional parameter. If given, the client will print out more
   information each time step.

This code should be fairly simple to extend to other pair-wise
interaction potentials, and examples of its use can be seen in the
“examples” directory, as explained in `2.4 <#tests>`__.

CP2K
^^^^

To use CP2K as the client code using an internet domain socket on the
host address “host_address” and on the port number “port” the following
lines must be added to its input file:

If instead a unix domain socket is required then the following
modification is necessary:

The rest of the input file should be the same as for a standard CP2K
calculation, as explained at `www.cp2k.org/ <www.cp2k.org/>`__.

Quantum-Espresso
^^^^^^^^^^^^^^^^

To use Quantum-Espresso as the client code using an internet domain
socket on the host address “host_address” and on the port number “port”
the following lines must be added to its input file:

If instead a unix domain socket is required then the following
modification is necessary:

The rest of the input file should be the same as for a standard Quantum
Espresso calculation, as explained at
`www.quantum-espresso.org/ <www.quantum-espresso.org/>`__.

LAMMPS
^^^^^^

To use LAMMPS as the client code using an internet domain socket on the
host address “host_address” and on the port number “port” the following
lines must be added to its input file:

If instead a unix domain socket is required then the following
modification is necessary:

The rest of the input file should be the same as for a standard LAMMPS
calculation, as explained at http://lammps.sandia.gov/index.html. Note
that LAMMPS must be compiled with the ``yes-user-misc`` option to
communicate with i-PI. More information from
https://lammps.sandia.gov/doc/fix_ipi.html.

FHI-aims
^^^^^^^^

To use FHI-aims as the client code using an internet domain socket on
the host address “host_address” and on the port number “port” the
following lines must be added to its ``control.in`` file:

If instead a unix domain socket is required then the following
modification is necessary:

One can also communicate different electronic-structure quantities to
i-PI through the ``extra`` string from FHI-aims. In this case the
following lines can be added to the ``control.in`` file:

where option can be, e.g.,
``dipole, hirshfeld, workfunction, friction``.

.. _hpc:

Running on a HPC system
~~~~~~~~~~~~~~~~~~~~~~~

Running i-PI on a high-performance computing (HPC) system can be a bit
more challenging than running it locally using UNIX-domain sockets or
using the *localhost* network interface. The main problem is related to
the fact that different HPC systems adopt a variety of solutions to have
the different nodes communicate with each other and with the login
nodes, and to queue and manage computational jobs.

.. figure:: ../figures/ipi-running.*
   :width: 90.0%

   Different approaches to run i-PI and a number of
   instances of the forces code on a HPC system: a) running i-PI and the
   clients in a single job; b) running i-PI and the clients on the same
   system, but using different jobs, or running i-PI interactively on
   the login node; c) running i-PI on a local workstation, communicating
   with the clients (that can run on one or multiple HPC systems) over
   the internet.

Figure `2.1 <#fig:running>`__ represents schematically three different
approaches to run i-PI on a HPC system:

#. running both i-PI and multiple instances of the client as a single
   job on the HPC system. The job submission script must launch i-PI
   first, as a serial background job, then wait a few seconds for it to
   load and create a socket

   Then, one should launch with mpirun or any system-specific mechanism
   one or more independent instances of the client code. Note that not
   all queing systems allow launching several mpirun instances from a
   single job.

#. running i-PI and the clients on the HPC system, but in separate jobs.
   Since i-PI consumes very little resources, one should ideally launch
   it interactively on a login node

   or alternative on a queue with a very long wall-clock time. Then,
   multiple instances of the client can be run as independent jobs: as
   they start, they will connect to the server which will take care of
   adding them dynamically to the list of active clients, dispatching
   force calculations to them, and removing them from the list when
   their wall-clock time expires. This is perhaps the model that applies
   more easily to different HPC systems; however it requires having
   permission to run on the head node, or having access to a long
   wall-clock time queue that ensures that i-PI is always active.

#. running i-PI on a simple workstation, and performing communication
   over the internet with the clients that run on one or more HPC
   systems. This model exploits in full the distributed-computing model
   that underlies the philosophy of i-PI and is very robust – as the
   server can be always on, and the output of the simulation is
   generated locally. However, this is also the most complicated to set
   up, as the local workstation must accept in-coming connections from
   the internet – which is not always possible when behind a firewall –
   and the compute nodes of the HPC centre must have an outgoing
   connection to the internet, which often requires ssh tunnelling
   through a login node (see section `3.3 <#distrib>`__ for more
   details).

.. _tests:

Testing the install
-------------------

test the installation with ‘nose‘
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are several test cases included, that can be run automatically
with ‘nosetests‘ from the root directory.

.. code-block::

   > nosetests -v

test cases with input files
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Several test cases are distributed with the code to ensure that your
distribution is working correctly. There are also simple tests to see if
the client codes are working correctly.

All the input files are contained in the directory “examples”, which is
subdivided into the following directories:

tutorial:
   Contains the input files needed to run the tutorial in
   `5 <#tutorial>`__.

lj:
   This gives a simple classical Lennard-Jones simulation of Ne. The
   state points are given by (:math:`N`, :math:`\rho`, :math:`T`) =
   (864, 0.35, 1.62), (:math:`N`, :math:`\rho`, :math:`T`) = (864, 0.75,
   1.069) and (:math:`N`, :math:`\rho`, :math:`T`) = (864, 0.88, 1.095)
   in reduced Lennard-Jones units, so that the results can be compared
   to those in :cite:`lverlet67pr`.

ph2:
   This simulates para-hydrogen using the isotropic Silvera-Goldman pair
   potential :cite:`silv-gold78jcp`. There are three
   directories, “RPMD”, “nvt” and “Tuckerman”. “RPMD” and “nvt” have
   tests which can be compared to the results of
   :cite:`mill-mano05jcp`, and “Tuckerman” has tests which
   can be compared to the results of :cite:`mart+99jcp`.

qespresso:
   This has two simple examples to test to see if the Quantum-Espresso
   client is functioning correctly. There is one simple 4-atom lithium
   test, and a test using a single water molecule.

harmonic:
   This has a simple example of a 1D harmonic oscillator. This
   demonstrates the displaced path integral momentum distribution
   estimator as given in :cite:`lin+10prl`. As the momentum
   distribution is known analytically for this simple system, this
   provides an indication of how well the method is working.

lammps:
   This has a simple implementation of the q-TIP4P-F empirical water
   model of :cite:`habe+09jcp` using the classical molecular
   dynamics code LAMMPS. It demonstrates both the convergence of the
   PIGLET method :cite:`ceri-mano12prl`, as well as the use
   of ring-polymer contraction methods
   :cite:`mark-mano08jcp`.

   This also contains one example using LAMMPS to calculate the
   interactions between carbon atoms in graphene. This uses the
   optimized Tersoff parameters for carbon given in
   :cite:`lind-broi10prb`.

cp2k:
   Contains the tests for the CP2K client code. Holds input files to run
   the high-pressure water calculations presented in [cpc publication
   citation].

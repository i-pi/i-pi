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

-  Python version 3.6 or greater (starting from version 2.0, i-PI is not Python2
   compatible)
-  The Python numerical library NumPy

See :ref:`runningsimulations` for more details on how to launch i-PI.

Additionally, most client codes will have their own requirements. Many
of them, including the test client codes given in the “drivers”
directory, will need a suitable Fortran compiler. A C compiler is
required for the sockets.c wrapper to the sockets standard library. Most
electronic structure codes will also need to be linked with some
mathematical libraries, such as BLAS, FFTW and LAPACK. Installation
instructions for these codes should be provided as part of the code
distribution and on the appropriate website. 
Patching for use with i-PI should not
introduce further dependencies.


Quick Setup
~~~~~~~~~~~ 

To use i-PI with an existing driver, install and update using `pip`:

Last version:

.. code-block::
   
   python -m pip install git+https://github.com/i-pi/i-pi.git

Last Release:

.. code-block::

   pip install -U ipi


Full installation
~~~~~~~~~~~~~~~~~ 

i-PI
^^^^

To develop i-PI or test it with the self-contained driver, follow these
instructions. It is assumed that i-PI will
be run from a Linux environment, with a recent version of Python, Numpy and
gfortran, and that the terminal is initially in the i-pi package directory (the
directory containing this file), which you can obtain by cloning the repository

.. code-block::

   git clone https://github.com/i-pi/i-pi.git


Source the environment settings file `env.sh` as `source env.sh` or `.
env.sh`.  It is useful to put this in your `.bashrc` or other settings file if
you always want to have i-PI available.


Fortran built-in driver
^^^^^^^^^^^^^^^^^^^^^^^

The built-in driver requires a FORTRAN compiler, and can be built as

.. code-block::
   
   cd drivers/f90
   make
   cd ../..

Python driver and PES
^^^^^^^^^^^^^^^^^^^^^

In addition to the FORTRAN drive, the i-PI distribution contains also a Python 
driver, available in `drivers/py` and through the command-line command 
`i-pi-py_driver`, which evaluates potential energy surfaces evaluated by simple 
driver classes, that can be found in `ipi/pes`. 

These classes are particularly suitable to perform inference with machine-learning
potentials implemented in Python, and it is reasonably simple to add your own,
if you need to (see also the :ref:`contributing` section).

These PES files can also be used directly, without the need to go through a 
client-server interface, using a :ref:`ffdirect` forcefield, including in the
XML input a block similar to 

.. code-block::
   
   <ffdirect name="ff_name">
      <pes> harmonic </pes>
      <parameters> { k1: 1.0} </parameters>
   </ffdirect>


Alternative installation using the setup.py module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is sometimes  more convenient to install the package to the system’s Python
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
administrator privileges.

Otherwise, one can install i-PI in a local Python path. If such path
does not exist yet, one must create directories for the package to go
into, using:

.. code-block::

  > mkdir ~/bin
  > mkdir ~/lib/py_vers
  > mkdir ~/lib/py_vers/site-packages

Next, you must tell Python where to find this library, by appending it to
the Linux environment variable PYTHONPATH, using:

.. code-block::

  export PYTHONPATH=$PYTHONPATH:~/lib/py_vers/site-packages/

Finally, the code can be installed using:

.. code-block::

 > python setup.py install –prefix= 

Either way, it will now be possible to run the code automatically, using

.. code-block::

 > i-pi input-file.xml


Installing clients
------------------

As of today, the following codes provide out-of-the-box an i-PI
interface: CP2K, DFTB+, Lammps, Quantum ESPRESSO, Siesta, FHI-aims,
Yaff, deMonNano, TBE. Links to the web pages of these codes, including
information on how to obtain them, can be found at http://ipi-code.org/.

If you are interested in interfacing your code to i-PI please get in
touch, we are always glad to help. We keep some information below in
case you are interested in writing a patch to a code.

.. _runningsimulations:

Running i-PI
------------

i-PI functions based on a client-server protocol, where the evolution of
the nuclear dynamics is performed by the i-PI server, whereas the energy
and forces evaluation is delegated to one or more instances of an
external program, that acts as a client. This design principle has
several advantages, but it also makes running a
simulation slightly more complicated, since the two components must be
set up and started independently.

Running the i-PI server
~~~~~~~~~~~~~~~~~~~~~~~

i-PI simulations are run using the i-PI Python script found in the
“bin” directory. This script takes an xml-formatted file as input, and
automatically starts a simulation as specified by the data held in it.
If the input file is called “input_file.xml”, then i-PI is run using:

.. code-block::

    > python i-pi input_file.xml

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

Built-in, fortran client
^^^^^^^^^^^^^^^^^^^^^^^^

i-PI includes a Fortran empirical potential
client code to do simple calculations and to run the examples.

The source code for this is included in the directory “drivers/f90”, and can
be compiled into an executable “i-pi-driver” using the UNIX utility
make.

This code currently has several empirical potentials hardcoded into it, including
a Lennard-Jones potential, the Silvera-Goldman potential
:cite:`silv-gold78jcp`,
a primitive implementation of the  qtip4pf potential for water ,
:cite:`habe+09jcp`,
several toy model potentials,
the ideal gas (i.e. no potential interaction), and several more.

How the code is run is based on what command line arguments are passed
to it. The command line syntax is:

.. code-block::

   > i-pi-driver [-u] -a address [-p port] -m [model-name] -o [parameters] [-S sockets_prefix] [-v] 


The flags do the following:

-u:
   Optional parameter. If specified, the client will connect to a Unix
   domain socket. If not, it will connect to an internet socket.

-a:
   Is followed in the command line argument list by the hostname (address) of the
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
   Is followed in the command line argument list by a string of comma-separated values needed to initialize the potential parameters. “gas”
   requires no parameters, “harm” requires a spring constant, “sg”
   requires a cut-off radius and “lj” requires the length and energy
   scales and a cut-off radius to be specified. All of these must be
   given in atomic units.

-v:
   Optional parameter. If given, the client will print out more
   information each time step.

-S:
   Optional parameter. If given, overwrite the default socket prefix used in the creation of files for the socket communication.
   (default "/tmp/ipi\_")

This code should be fairly simple to extend to other pair-wise
interaction potentials, and examples of its use can be seen in the
“examples” directory, as explained in :ref:`tests`.

CP2K
^^^^

To use CP2K as the client code using an internet domain socket on the
host address “host_address” and on the port number “port” the following
lines must be added to its input file:

.. code-block::

    &GLOBAL
       ...
       RUN_TYPE DRIVER
       ...
    &END GLOBAL

    &MOTION
       ...
       &DRIVER
          HOST host_address
          PORT port
       &END DRIVER
       ...
    &END MOTION

If instead a Unix domain socket is required then the following
modification is necessary:

.. code-block::

    &MOTION
       ...
       &DRIVER
          HOST host_address
          PORT port
          UNIX
       &END DRIVER
       ...
    &END MOTION

The rest of the input file should be the same as for a standard CP2K
calculation, as explained at it the documentation of 
`CP2K <http://www.cp2k.org/>`__.

Quantum-Espresso
^^^^^^^^^^^^^^^^

To use Quantum-Espresso as the client code using an internet domain
socket on the host address “host_address” and on the port number “port”
the following lines must be added to its input file:

.. code-block::

    &CONTROL
       ...
       calculation=`driver'
       srvaddress=`host_address:port'
       ...
    /

If instead a Unix domain socket is required then the following
modification is necessary:

.. code-block::

    &CONTROL
       ...
       calculation=`driver'
       srvaddress=`UNIX:socket_name:port'
       ...
    /
    
The rest of the input file should be the same as for a standard Quantum
Espresso calculation, as explained at
`www.quantum-espresso.org/ <www.quantum-espresso.org/>`__.

LAMMPS
^^^^^^

To use LAMMPS as the client code using an internet domain socket on the
host address “host_address” and on the port number “port” the following
lines must be added to its input file:

.. code-block::

    fix  1 all ipi host_address port

If instead a unix domain socket is required then the following
modification is necessary:

.. code-block::

    fix  1 all ipi host_address port unix

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

.. code-block::

    use_pimd_wrapper host_address port

If instead a unix domain socket is required then the following
modification is necessary:

.. code-block::

    use_pimd_wrapper UNIX:host_address port

One can also communicate different electronic-structure quantities to
i-PI through the ``extra`` string from FHI-aims. In this case, the
following lines can be added to the ``control.in`` file:

.. code-block::

    communicate_pimd_wrapper option
    
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

.. _fig-running:

.. figure:: ../figures/ipi-running.*
   :class: white-background   
   :width: 90.0%

   Different approaches to run i-PI and a number of
   instances of the forces code on a HPC system: a) running i-PI and the
   clients in a single job; b) running i-PI and the clients on the same
   system, but using different jobs, or running i-PI interactively on
   the login node; c) running i-PI on a local workstation, communicating
   with the clients (that can run on one or multiple HPC systems) over
   the internet.

The figure represents schematically three different
approaches to run i-PI on a HPC system:

#. running both i-PI and multiple instances of the client as a single
   job on the HPC system. The job submission script must launch i-PI
   first, as a serial background job, then wait a few seconds for it to
   load and create a socket


    .. code-block::

        > python i-pi input_file.xml &> log & wait 10    

   Then, one should launch with mpirun or any system-specific mechanism
   one or more independent instances of the client code. Note that not
   all queing systems allow launching several mpirun instances from a
   single job.

#. running i-PI and the clients on the HPC system, but in separate jobs.
   Since i-PI consumes very little resources, one should ideally launch
   it interactively on a login node
   
   .. code-block::

        > nohup python i-pi input_file.xml < /dev/null &> log &

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
   through a login node (see section :ref:`distrib` for more
   details).

.. _tests:

Testing the install
-------------------

test the installation with ‘pytest‘
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are several test cases included, that can be run automatically
with ‘i-pi-tests‘ from the root directory.

.. code-block::

   > i-pi-tests

Note 1: pytest and pytest-mock python packages are required to run these tests, but they are required to run i-PI.
Note 2: please compile the fortran driver, as explained in :ref:`driver.x`.
Note 3: use the '-h' flag to see all the available tests

test cases and examples
~~~~~~~~~~~~~~~~~~~~~~~

The `examples/` folder contains a multitude of examples for i-PI, covering
most of the existing functionalities, and including also simple tests that
can be run with different client codes. 


The example folder is structured such that each sub-folder is focused on a different aspect of using i-PI:

- **clients**: 
    Contains examples that are code-specific, highlighting how the driver code should be set up
                    (client-specific syntax and tags) to run it properly with i-PI

- **features** :  
     Examples of different functionalities implemented in i-PI.
                    All examples can be run locally with the drivers provided with the code.

- **hpc_scripts** :  
      Examples of submission scripts on HPC platforms

- **temp**     :
     Temporary folder with historical examples that have not yet been adapted
                    to the current folder structure

- **init_files**: 
      repository of input files shared by many examples

We keep this folder updated as much as we can, and try to run automated tests on these inputs, but in some cases, e.g. when using external clients, we cannot run tests.
Please report a bug if you find something that is not working.

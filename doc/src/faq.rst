Frequently asked questions
==========================

Where to ask for help?
----------------------
There is distinction between questions and bugs.
If your question is related to running i-PI, using or understanding its features
(also if you think that something is not documented enough),
the best way to ask is to post your question on the i-PI forum: https://groups.google.com/g/ipi-users.
In the forum, everything discussed will be available to other members of the community in the future!

If you are seeing an actual bug,
please consider reporting an issue on Github: https://github.com/i-pi/i-pi/issues.

Which units does i-PI use?
---------------------------
Atomic units are used everywhere inside the code.
The exact values are listed in the manual section "Internal units".
For I/O operations, various units are accepted.
The XML tags responsible for physical settings accept the attribute ``units=<unit>``.
If nothing is specified, i-PI will assume atomic units.
Many of the provided examples have these tags,
and the reference section of the manual contains a full list.
Another way to use units is related to input ``.xyz`` files. i-PI can accept declarations of these units in the ``.xyz`` comment line. One example is:

.. code-block:: bash 

  # CELL(abcABC): a b c alpha beta gamma cell{units} positions/forces/velocities{units}

Another place worth looking at is ``<i-pi-root>/ipi/utils/units.py``,
which contains definitions and values for non-default units,
such as ``electronvolt``, ``inversecm`` etc.

How to build i-PI?
------------------
i-PI is a Python code, and as such, strictly speaking, does not need to be compiled and installed.
``<i-pi-root>/bin/i-pi`` file is an executable.
``source <i-pi-root>/env.sh`` ensures that necessary paths are added to PATH, PYTHONPATH etc.
However, it is often more convenient to install the package to the system’s Python modules path,
so that it is accessible by all users and can be run without specifying the path to the Python script.
For this purpose we have included a module in the root directory of the i-PI distribution, ``setup.py``,
which handles creating a package with the executable and all the modules which are necessary for it to run.
Detailed explanation can be found in the manual.

In addition, i-PI can be installed through ``pip`` with the command

.. code-block:: bash

   $ pip install -U i-PI

And to test it with pytest

.. code-block:: bash

   $ pip install pytest
   $ pytest --pyargs ipi.tests


How to run i-PI with the client code \<name\>?
----------------------------------------------
i-PI communicates with electronic structure and MD codes via socket communication protocol.
For many popular codes (CP2K, FHI-aims, LAMMPS, Quantum ESPRESSO, VASP, etc.)
we provide examples in ``<i-pi-root>/examples`` folder.

Another way of connecting to client codes is using the ASE client.
This way you get access to wide variety of codes that are connected to ASE,
but for some of them current implementation requires restarting the client code after each step,
which may lead to a significant overhead in case of electronic structure calculations.
We recommend using the socket connection from the client code to ASE and then from ASE to i-PI, when possible. An example of a "double-socket" setup
can be found in:

``https://github.com/i-pi/i-pi/tree/master/examples/ASEClient``

How to run i-PI on a cluster?
-----------------------------
There are different ways of running i-PI on HPC systems,
described in details in the manual.
One simple setup, which assumes that i-PI and all client codes
are launched from the same Slurm script, is given in ``examples/slurm/sl-rpc``.
In this example, i-PI is launched in the background first,
and then its hostname is passed to the input files of the clients.
Clients can communicate to the same socket (identified by the same port number)
or to different sockets (different port numbers) that i-PI has opened.
The clients are subsequently launched from the same slurm script.

Other setups are possible, please consult the manual.

How to setup colored-noise thermostats with meaningful parameters?
------------------------------------------------------------------
If you have an application that will profit from colored noise thermostats, you can download the parameters from http://gle4md.org.
This website contains sets of input parameters optimized for a variety of conditions. You can choose i-PI syntax, such that you will only have
to copy-paste that section inside the i-PI input file.
Also the ``<i-pi-root>/examples`` (especially withing ``lammps``) folder demonstrates the syntax for various thermostats.

How to perform geometry optimization?
-------------------------------------
Several examples of geometry optimizaion, including fixing certain atoms,
are present in the ``<i-pi-root>/examples`` folder.
The folders with names \*geop\* correspond to geometry optimization.
Note that the optimization algorithms that involve line search have two sets of convergence parameters:
one for line search and another for the main algorithm.

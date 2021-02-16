Frequently asked questions
==========================

Where to ask for help?
----------------------
There is distinction between questions and bugs.
If your question is related to running i-PI, using or understanding its features
(also if you think that something is not documented enough),
the best way to ask is to post your question on the i-PI forum: https://groups.google.com/g/ipi-users.
We would like to discourage you from writing personal e-mails to developers,
so that other members of the community could read the conversation.

If you are sure that the output you get is different from what is declared,
please consider reporting an issue on Github: https://github.com/i-pi/i-pi/issues.

How to run i-PI with the client code \<name\>?
----------------------------------------------
i-PI communicates with electronic structure and MD codes via socket communication protocol.
For many popular codes (CP2K, FHI-aims, LAMMPS, Quantum ESPRESSO, VASP etc)
we provide examples in ``<i-pi-root>/examples/`` folder.
Another way of connecting to client codes is using ASE client.
This way you get access to wide variety of codes that are connected to ASE,
but for some of them current implementation requires restarting force code after each step,
which may lead to significant overhead in case of electronic structure calculations.
We recommend using direct connections when possible.

How to setup colored-noise thermostats with meaningful parameters?
------------------------------------------------------------------
We warn you against using complicated thermostats before reading
the papers describing underlying physics.
If you already know what you are doing, you may find http://gle4md.org useful.
This website contains sets of input parameters optimized for a variety of conditions.
Also ``<i-pi-root>/examples/`` folder demonstrates the syntax for various thermostats.

How to perform geometry optimization?
-------------------------------------
Several examples of geometry optimizaion including fixing certain atoms
are present in the ``<i-pi-root>/examples/`` folder.
The folders with names \*geop\* correspond to geometry optimization.
Note that the optimization algorithms that involve line search have two sets of convergence parameters:
one for line search and another for the main algorithm.

How to build i-PI?
------------------
i-PI is Python code, and as such, strictly speaking, does not need to be compiled and installed.
``<i-pi-root>/bin/i-pi`` file is an executable.
However, it is often more convenient to install the package to the system’s Python modules path,
so that it is accessible by all users and can be run without specifying the path to the Python script.
For this purpose we have included a module in the root directory of the i-PI distribution, ``setup.py``,
which handles creating a package with the executable and all the modules which are necessary for it to run.
Detailed explanation can be found in the manual.

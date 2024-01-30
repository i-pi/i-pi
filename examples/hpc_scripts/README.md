i-PI on HPC
===========

This folder gives examples of setup to run i-PI on high-performance computing systems.
Running on HPC tends to be a sore point with i-PI, because one has to handle both
the driver code (which is the one usually consuming serious CPU resources), i-PI
(which is usually inexpensive, but still has to run somewhere) and the communication,
that may depend on the network setup on the cluster.

You may get some ideas on the difficulties involved in the 
[distributed execution section](https://ipi-code.org/i-pi/user-guide.html#distributed-execution)
of the manual.  This folder is meant to provide examples of submission scripts
that can help you set up calculations on a typical HPC environment. 
You will probably have to make several adjustments to compensate for the quirks
of your setup: if you succeed and think that your script could help others, you
are encouraged to contribute it.


Single-node SLURM setup
=======================

This folder provides an example of how to set up an i-PI calculation to run 
on a single node using the SLURM queue manager. The simulation is just a PIMD
run for the Zundel cation, and hardly needs HPC resources, but it illustrates
the mechanism and provides a template for more substantial calculations.

The defining feature for this setup is that communication happens using 
unix-domain sockets, which are a shared-memory inter-process communication
mechanism, and so require all the jobs to run on the same node.

While these scripts should run on a basic SLURM queue setup, keep in mind that
you may have to customize them to use them in your own case - from setting the
appropriate environment variables, to changing the number of cores per job. 
If you need particular workarounds to make it function on your system, feel
free to open an issue so we can perhaps cover more use cases. 

The i-PI run 
------------

The basic i-PI execution is not dissimilar to running in a normal setting:
`i-pi` is launched, and then several instances of the force driver 
(in this case the FORTRAN `i-pi-driver`) are executed so that they can
connect. Everything is run in the background, and then the process waits
for all the spawned runs to terminate.

```bash
i-pi input.xml &
sleep 5
for nbead in `seq 1 8`; do
   i-pi-driver -a slurm-one-node -m zundel -u -v &
done
wait 
```

Note the `sleep 5` command, that waits 5 seconds to ensure that i-PI has the
time to load and to open the socket before the drivers start connecting. 

Simple SLURM script
-------------------

The `job_base.sh` script contains an example on how to run this job
in the most basic setting: i-PI runs for just one time, and there is 
no naming conflict between the socket used by this job and other jobs
potentially running on the same node. 

In this case, the only additional commands are those used to set up 
the SLURM job, and the `srun` call prepended to the different commands, 
which dispatches the jobs to the compute node.

The main caveat here is that i-PI should end *before* the allocated time
finishes. To do that, it is important to add the `<total_time>` block
inside `input.xml`, specifying a number of seconds that is smaller than 
the time requested to SLURM, so that the calculation can be terminated
gracefully even if it has not ended: in this case i-PI will write a 
`RESTART` file that can be used to continue the run.

Advanced SLURM script
---------------------

The `job_advanced.sh` script attempts to cover some of the possible 
issues that can arise when running on an HPC system.  This somewhat
obfuscates the workflow, so it is good to first make sure that
the basic version works smoothly.

In particular, this script takes care of:

1. Restarting automatically the calculation if a `RESTART` file is found
2. Renaming the socket address to be a UUID, to avoid potential conflicts
   with  previous or concurrent runs

This means one can run calculations that do not finish in a single run, but
are executed in a chain. Remember to always run a test with a dummy potential
before setting up a complicated SLURM setup that burns serious CPU time. 
 

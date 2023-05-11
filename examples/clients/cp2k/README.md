CP2K - i-PI example
===================

*******************************************************************************
  NOTE: The examples in this directory correspond to the calculations 
  performed in Ceriotti, More, Manolopoulos, Comp. Phys. Comm. 2013.
  Contrary to most examples, the parameters used here are suitable 
  for production calculations (except nvt-cl), and so these examples
  are computationally demanding and should not be run on a workstation.
*******************************************************************************

Running the CP2K example
------------------------

 * First, it is necessary to patch and compile CP2K, and create a make.in
   file containing the path to the cp2k executable, e.g.

```
CP2K:=~/bin/cp2k.popt
```

** Run the examples automatically:
 
 * You can run classical and PIGLET simulations of 64 water molecules at 
   750K and 10GPa, both in the NPT and NVT ensembles. The make targets are 
   self-explanatory -- npt-classical, nvt-classical, npt-piglet, nvt-piglet. 
   To run the NPT PIGLET example, for instance, just type:

```bash
$ make npt-piglet
```

 * To clean up output files:

```bash
$ make clean
```

** Run the examples manually:
 
 * In the main example directory run 

```bash
$ python path/bin/i-pi input.xml
 ```

 * In another terminal, create up to 4 temporary directories to run cp2k,
 and launch it from within:

```bash 
$ mkdir run_1
$ cd run_1
$ $CP2K -i ../in.cp2k
 ```

 * Repeat several times to get multiple instances
 
 * To run CP2K remotely, you should:
    1. make sure that the machine you run i-pi on has a port accessible 
       from outside
    2. change the <interface> option in nptpgl_ipi.xml so that <address> points
       to the external network interface (i.e. give a fully-qualified domain 
       name or IP address rather than localhost) and that the <port>
       selected is accessible from outside
    3. change the ADDRESS and PORT input parameter in nptpgl_cp2k.in, 
       making them point to the same address and port
    4. copy the cp2k inputs to the remote machine (do not forget the basis sets!) 
       and launch n copies of $CP2K remotely

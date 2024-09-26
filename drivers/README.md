Driver code directory
=====================

This folder contains simple (or not so simple...) driver codes.
These serve both as examples of how to build an interface with i-PI,
and to run some of the examples. 

f90
---

This is a F90 implementation, with a `sockets.c` module that contains
an interface to the standard socket library. On a reasonably set-up
system containing a fortran compiler and libraries, it should be 
possible to generate the executable `driver.x` by just typing

```
make
```

py
--

A Python version of the driver. Some simple potentials, and interfaces
with codes that have a Python API can be found in the `pes/` 
folder. 



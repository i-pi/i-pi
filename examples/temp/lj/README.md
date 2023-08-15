Example with the rudimentary LJ code
====================================

This gives an example of Neon with a Lennard-Jones potential pair potential,
which can be compared to the following data from Loup Verlet, Phys. Rev. 159,
98, (1967).

- low_density:
  State point (N, rho, T) = (864, 0.35, 1.620)
- mid_density:
  State point (N, rho, T) = (864, 0.75, 1.069)
- high_density:
  State point (N, rho, T) = (864, 0.88, 1.095)
- nst:
  Constant stress example for high density (rho=0.88) and low temperature
  (T=0.54305 in reduced units, or 20 K) - produces a solid, in agreement with
  Verlet.

Note that the "high_density" case appears to illustrate a typo in the paper,
since the potential energy as found by this wrapper and by the LAMMPS MD code
are both -5.86 rather than -5.66 as reported in the original paper.


Run the examples automatically
------------------------------

These can be run automatically using the Makefile provided. The make targets
are self-explanatory - high_density, mid_density and low_density.  To run the
high density example, for instance, just type:

```
$ make high_density
```

To clean up output files:

```
$ make clean
```


Run the examples manually
-------------------------

Go back to the example directory and run

```
$ i-pi input.xml
```

the wrapper will start and sit waiting on the UDS `/tmp/ipi`.

Open a separate terminal and run the LJ driver code. Each of the tests runs
using a different server address, which are identical to their directory names.
For example, the high density problem is run using:

```
$ i-pi-driver -u -h high_density -m lj -o 5.270446,1.1663e-4,13.176115
```

You can run multiple instances of the code; it is so fast that parallel scaling
won't be appreciable.

If your system does not support Unix domain sockets, just set in `input.xml`
`<socket mode="inet"> <port> port_no </port>`

then run the driver as:

```
$ i-pi-driver -h high_density -m lj -o 5.270446,1.1663e-4,13.176115 -p port_no
```

Using the i-PI driver with (high quality) potential energy surfaces
===================================================================

Zundel cation
-------------

The `zundel/` folder contains a simple example of how to run the i-PI driver
together with a potential fit for the Zundel cation in the gas phase based on
high-end quantum chemistry methods. 

Once the env.sh file in the i-pi root directory has been sourced with

```bash
source <i-pi-root>/env.sh
```

One simply needs to run i-pi with

```bash
i-pi input.xml
```

and then run one or more instances of the driver code that contains a routine
by Xinchuan Huang, Bastiaan J. Braams, and Joel M. Bowman, [J. Chem. Phys. 122,
044308 (2005)] to compute the energy for H5O2+. Note that the data files
`h5o2.pes4B.coeff.dat` and `h5o2.dms4B.coeff.com.dat` should be present in the
folder one wants to run the driver code. 

```bash
i-pi-driver -u -h zundel -m zundel
```


qTIP4P/f water model
--------------------

The `qtip4pf/` folder contains an example of how to run the i-PI driver with a
simplistic implementation of the q-TIP4P/f water model.

One runs i-PI as usual, followed by one or more instances of the driver:

```bash
i-pi input.xml &
i-pi-driver -u -h driver -m qtip4pf
```

Remember that the `<i-pi-root>/env.sh` file must be sourced before running
`i-pi` or `i-pi-driver`:

```bash
source <i-pi-root>/env.sh
```

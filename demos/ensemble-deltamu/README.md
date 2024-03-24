Solid-liquid free energy with uncertainty quantification
========================================================

This example demonstrates a rather sophisticated setup to compute the $\Delta\mu_{SL}$ 
for the ice/water system using an ensemble of ML potentials, including an evaluation of
the uncertainty using a direct propagation of the ensemble (cf. [Imbalzano2021, Kellner2024]).
NB: while these simulations are realistic, the system size is too small to be converged.

Prerequisites
-------------

Besides a recent version of i-PI, this demonstration requires some external dependencies.
In particular you will need to install [Plumed](https://plumed.org), making sure to enable
the Python interface and the `crystallization` module.

```
./configure ---enable-python --enable-modules +crystallization
make
```

You may need to install additional packages, follow the documentation of Plumed.

The MLIP used in this example is XXXXXXXXXX MATTHIAS INCLUDE HERE INSTRUCTION TO DOWNLOAD
AND INSTALL/CONFIGURE


Interface pinning
-----------------

[Interface pinning](https://doi.org/10.1103/physrevb.88.094101) is a method to compute the 
free-energy difference between two phases (typically liquid and solid, which we consider here) 
by performing simulations of a supercell containing two regions of the two phases, separated 
by planar interfaces.
The total free energy of the system is then proportional to the fraction $x$ of one of the two
phases, the proportionality constant being the molar difference in free energy between 
solid and liquid, $\Delta \mu$. 
The solid fraction does not have have to be computed explicitly, but can be assessed by using
any atom-centered descriptor that distinguishes the two phases. This example uses the 
Steinhardt parameter $Q_6$, computed by the [Plumed package](https://plumed.org).
The total $Q_6$ computed over the entire box for a configuration $A$, $Q(A)$ is then used 
to compute a _pinning potential_ $W(A)=k/2(Q(A)-\bar{Q})^2$. 
  
Assuming that the mean $Q(A)$ is proportional to the solid fraction (see e.g. the discussion
in [Cheng et al.](http://doi.org/10.1103/PhysRevB.92.180102) for why that would be the case)
the deviation of $\langle Q(A) \rangle$ from the  target value is equal to the 
_total_ difference in free-energy for the box when the mean order parameter changes, so

$$
\Delta \mu_{SL} \propto k (\langle Q(A) \rangle - \bar{Q})
$$

the proportionality constant depending on the order parameter and the box geometry.
In this example we use a simple trick to obtain a quantity that yields directly the chemical
potential: we define a modified order parameter that counts the number of atoms above
a set threshold, which - given that $Q_6$ separates well between the two phases - 
provides a good proxy for the number of solid-like atoms. 

This setup is implemented using the `<ffplumed>` forcefield, that requires installing 
the [plumed library](https://plumed.org). Follow the installation instructions for 
plumed, and make sure to also compile and install the python module. 
The relevant section in the `input.xml` file is

```
  <ffplumed name="plumed">
      <file mode="xyz">init.xyz</file>
      <plumeddat>plumed.dat </plumeddat>
  </ffplumed>
```

that recalls the plumed input file `plumed.dat`. Besides the definition of the order parameter, 
the pinning restraint is defined as 

```
RESTRAINT ARG=m6.morethan-1 AT=165 KAPPA=0.01 LABEL=res
```

and the values of the order parameter and the bias, that at present 
are not communicated back to i-PI, are printed out for further processing. 
Note that the stride for this output should match that used in i-PI to print 
out trajectory statistics, as the two pieces of information must be combined
to complete the analysis

```
PRINT STRIDE=50 ARG=m6.*,res.bias FILE=COLVAR
FLUSH STRIDE=50
```

(Shallow) ensembles
-------------------

This example uses a committee of machine-learning potentials. This means that the potential energy
(and the forces) used to propagate the dynamics are obtained by averaging over multiple models, 
that differ by the training details and/or the architecture, and are calibrated to yield a spread
in the predictions that estimates the uncertainty in the energy prediction. 
You can read [Musil2019] for an early application to atomistic modeling, and [Kellner2024] for
a description of the model used here, that uses weight sharing to reduce the cost of evaluation. 

[** installation instructions for the model **]








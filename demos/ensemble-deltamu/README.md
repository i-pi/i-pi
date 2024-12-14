Solid-liquid free energy with uncertainty quantification
========================================================

Authors: 
`Matthias Kellner <matthias.kellner@epfl.ch>`
`Michele Ceriotti <michele.ceriotti@epfl.ch>`

This example demonstrates a rather sophisticated setup to compute the $\Delta\mu_{SL}$ 
for the ice/water system using an ensemble of ML potentials, including an evaluation of
the uncertainty using a direct propagation of the ensemble (cf. [Imbalzano2021, Kellner2024]).
NB: while these simulations are realistic, the system size is too small to be converged.

Simulations can be performed with either a deep ensemble of 4 separate Behler-Parrinello
neural network, computed with the n2p2 driver in `lammps`, or with a shallow ensemble of 5
weight-sharing multi-layer perceptrons based on SOAP features, computed with `rascaline`

Installation of the dependencies
-------------------------------
For convenience we provide a script to install the necessary dependencies for both the "full" committee `n2p2`/`lammps`/`plumed` version of this tutorial, but also for the shallow ensemble `rascaline`/`plumed` version. 

### Important information for Mx MacOS users:

- conda prebuilt versions of plumed for Aarch64 do not include the `crystallization` module, you will have to compile PLUMED from source.

Installation of the rascaline model
-----------------------------------

Besides a recent version of i-PI, this demonstration requires some external dependencies.
In particular you will need to install [Plumed](https://plumed.org), making sure to enable
the Python interface and the `crystallization` module - prebuilt versions from conda do not
include this module.

```
./configure ---enable-python --enable-modules +crystallization
make
```

You may need to install additional packages, follow the documentation of Plumed.
The install script discussed below will try to install Plumed in the conda environment
but you mileage may vary. You may also need to set manually `PLUMED_PATH` to detect it. 

The MLIP used in this example is a [BPNN](https://doi.org/10.1103/PhysRevLett.98.146401)-type neural network, 
using [SOAP](https://doi.org/10.1103/PhysRevB.87.184115) descriptors and Radial Spectrum descriptors as input.

In order to install the necessary packages for the MLIP we provide a setup script.
IMPORTANT: We recommend to install the packages in a fresh virtual environment.

```
./install_rascaline.sh
```

The install script will install [rascaline](https://github.com/Luthaf/rascaline) the atomic 
descriptor library, [metatensor](https://github.com/lab-cosmo/metatensor) for metadata handling 
and [torch](https://pytorch.org/) for the neural network training. 

Install from within the `./demos/ensemble-deltamu` directory.

Installation of the n2p2 model
------------------------------

You will need a version of LAMMPS patched with the N2P2 driver, and to download the models.\
Alternatively, we provide a script to install the necessary dependencies for the n2p2 model, including LAMMPS/n2p2 and PLUMED.\
The dependency resolution might take a while.

```
./install_n2p2.sh
```

Running the simulations
----------------------

The simulations can be run using the provided `run_rascaline.sh` script, from within the 
`./demos/ensemble-deltamu` directory.

```
./run_rascaline.sh
```

From within the `./demos/ensemble-deltamu` directory you can run the provided run-script to run the n2p2 demo after the conda-based installation.

```
./run_n2p2.sh
```

Alternatively, more geared towards running the n2p2 example on an HPC environment, set the environment variables in 
`run_n2p2_hpc.sh` `PLUMED_PATH` and `IPI_PATH` to the 
paths were the source of Plumed and i-PI are located and execute it.

```
./run_n2p2_hpc.sh
```


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
      <plumed_dat>plumed.dat </plumed_dat>
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

This demo uses the functionality of the i-pi `ff-committe` to interface with the MLIP, specifically 
using a driver that returns ALL the committee members as entries in a JSON-formatted
string. In here we highlight the utility of this functionality,
using a "Shallow Ensemble" Machine Learning interatomic potential.
A Shallow Ensemble is a committee of neural networks, where each
committee member implicitely shares the weights up to the last layer.
Passing the committee members predictions as a JSON string allows for
the efficient computation of committee predictions, as the recomputation of 
intermediate results (in particular the embeddings of the atomic environments 
and their gradients) can be avoided using only one calculator.


Interfacing a custom MLIP code with the `ff-committee` driver
-------------------------------------------------------------

To interface a custom MLIP code that yields efficiently an ensemble of predictions,
with the `ff-committee` driver, the JSON extras dict of the respective i-pi calculator 
must be formatted as follows:

```json
{
    "committee_pot": [pot_committee_1, pot_committee_2, ..., pot_committee_N],
    "committee_force": [pot_force_atom_1_x_committee_1, ..., pot_force_atom_1_x_committee_N, 
                        pot_force_atom_1_y_committee_1, ..., ...,
                        pot_force_atom_N_z_committee_N],
    "committee_virial": [virial_xx_committee_1, ..., virial_xx_committee_N, 
                        ..., virial_xy_committee_1, ..., ...,virial_zz_committee_N]
}
```

Or, suppose the committee are obtained as numpy arrays of shape (N_atoms, 3),

```python

# force_list is a list of numpy arrays of shape (N_atoms, 3)
# virial_list is a list of numpy arrays of shape (3, 3)
# pots is a numpy array of shape (N_committee, 1)

forces = np.stack([f.flatten() for f in force_list ] )
virials = np.stack([v.flatten() for v in virial_list ] )

# new arrays of shape (N_committee, N_components)

extras = {
    "committee_pot": pot_list.tolist(),
    "committee_force": forces.tolist(),
    "committee_virial": virials.tolist()
}
```

where `pot1`, `pot2`, ..., `potN` are the energies predicted by each committee,

NB: Since the committee predictions in the extras dict no longer pass through the
unit conversion of the driver, make sure that the calculator instance of your MLIP returns
i-pi units (Hartree, Bohr, etc.) and not ASE units (eV, Angstrom, etc.).

NB: Dependent on the exact architecture details of the MLIP, obtaining committee
forces and virials from a weight-shared ensemble might introduce a computational overhead.
It is therefore also possible to only return the committee energies and compute
the forces from the mean energies, which is not more expensive than computing the forces of
a single committee when using a weight-shared ensemble.

The calculator then should only return the committee energies in the extras dict.
Note that in this case force errors cannot be computed. 

```python

extras = {
    "committee_pot": pot_list.tolist()
}
```


Post-processing to evaluate the error on $\Delta \mu$
-----------------------------------------------------

The cumulant expansion approximation (CEA) provides a statistically-stable (but biased) 
expression for the thermodynamic average of a quantity $y$ over the ensemble generated by
the committee potential member $V^{(k)}$, starting from a group of configurations 
collected from a simulation performed with the mean potential $\bar{V}$. 

$$
\
\langle y \rangle_{V^{(k)}} \approx 
\langle y \rangle_{\bar{V}} - \beta 
\left[
\langle y (V^{(k)}-\bar{V}) \rangle_{\bar{V}} -
\langle y \rangle_{\bar{V}} \langle V^{(k)}-\bar{V} \rangle_{\bar{V}}
\right]
\
$$

It is important to keep in mind that the CEA is based on the assumption that $y$ and 
$V^{(k)}-\bar{V}$ are approximately joint Gaussian, and major deviations from this assumption
would lead to a very large bias, and possibly a meaningless estimate of the correction factors.

In our case, we want to compute the correction to $y= (Q(A) - \bar{Q})$; discarding the first 
1000 steps for equilibration, we can manipulate the `COLVAR` file to extract the correct entry
(look at `plumed.dat`, column 3 is the `morethan-1` entry, and 165 the $\bar{Q}$ reference 
value for the pinning potential). 

```bash
tail -n +2000 COLVAR | awk '{print $3-165}' > QQbar
```

The ensemble of physical potentials is printed out by i-PI into 
`pinning.committee_pot_0`, with some comment lines that can be eliminated for simplicity

```bash
grep -v \# pinning.committee_pot_0 | tail -n +2000 > POTS
```

NB: it is important that the two files are "synchronized", i.e. they use the same stride 
and start at the same point.

You can inspect the two files separately, or run the data analysis using the `i-pi-committee-reweight`
script, that performs a CEA (or direct reweighting) analysis based on these two files

```bash
i-pi-committee-reweight POTS QQbar --kt 0.00091837535
```

This will write out average and standard deviation of $\langle Q(A) - \bar{Q} \rangle$, as well
as the values for each ensemble.
By repeating this analysis for several temperatures, one can fit a linear relation linking $\Delta \mu$ 
to $T$, obtain the intercept and therefore an ensemble of melting temperatures, from which one can obtain
mean and standard error. 

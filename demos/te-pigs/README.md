# Te PIGS demonstration

Author: Venkat Kapil [v.kapil@ucl.ac.uk](v.kapil@ucl.ac.uk)

Here is a demonstration of the Te PIGS methods [1] on estimating first-principles quality IR and Raman spectra of liquid water at 300 K. It comprises two main steps:

1. <b>Training</b>: We will train a quantum effective potential corresponding to the potential of mean force of the centroid of the path integral using the Te PIGS method. The training data will be generated from a path-integral molecular dynamics simulation at an elevated temperature of 500 K. The effective potential will be fit using the MACE [2] machine learning framework. 
2. <b>Prediction</b>: This is as easy as performing molecular dynamics using the Te PIGS quantum effective potential at 300 K. IR and Raman spectra will be estimated as time correlation functions of the total polarization and polarizability within linear response theory. 

We employ the MACE-OFF(S) potential [3], a recently developed hybrid-functional DFT-level machine learning potential for general organic molecules, as the potential energy surface. To account for the IR and Raman selection rules, we use a [machine-learning model](https://github.com/venkatkapil24/ML-quantum-vibrational-spectroscopy) to predict the polarization and polarizability of bulk water [4].  

## Training 

### Step 1: Generating reference PIMD training data

Go to the directory ```0_reference_pimd``` and check out the ```input.xml``` file. We perform a path integral molecular dynamics simulation at 500 K using 8 imaginary time slices and a timestep of 0.5 fs. The total training data corresponds to 10 ps of simulation time. 

To fit the Te PIGS effective potential, we regress the centroid force against the centroid positions, but to make the fitting task easier, we can learn the difference between the centroid force and the physical force acting on the centroid position. The centroid positions and forces can be printed out using:

```xml
<trajectory filename='xc' stride='20' format='xyz' cell_units='ase'> x_centroid{ase} </trajectory>
<trajectory filename='centroid_force' stride='20' format='ase'> f_centroid </trajectory>
```

The physical force acting on the centroid can be printed out on the fly in two steps. First, we define a 'dummy' component of the system's force with a zero weight that receives the centroid position (set by ```nbeads='1'```) using:
```xml
<forces>
    <force forcefield='maceoff23' weight='1'> </force>
    <force forcefield='maceoff23' weight='0' nbeads='1'> </force>
</forces>
```
Since the weight of this component is zero, it doesn't affect the dynamics of the system. Then, we print out the raw value (i.e. without applying the weight) of the force component as 
```xml
<trajectory filename='physical_force' stride='20' format='ase'> forces_component_raw(1) </trajectory>
```

We conform to standard ```ASE``` units and print out the forces in extended xyz format that is easily read by ```ASE```. 

This simulation can be run by setting off the bash script:

```bash
bash run.sh
```

> **Note:**
> This step assumes that you have a working installation of i-PI using `pip` or have sourced the `env.sh` file in the top-level source directory.

### Step 2: Curating the dataset

Go to the directory ```1_dataset_curation``` and checkout the ```get_pigs_dataset.py``` file. This is a simple `ASE` based I/O parser which creates a set of atoms objects with positions set to centroid positions and the centroid and physical forces saved as ```centroid_force``` and ```physical_force``` arrays, respectively. Since we want to regress on the difference between the centroid and the physical forces, we store this quantity as a ```delta_force``` array. We save the dataset in an extended xyz format. 

> **Note:**
> This step assumes that your python environment has access to ```ase``` and ```ipi```. 


### Step 3: Fitting the Te PIGS potential

Go to the directory ```1_dataset_curation``` and check out the ```train.sh``` file. It's a standard ```MACE``` training script for training on forces only using  the  ```--energy_weight=0.0``` and ```--scaling='no_scaling'``` flags. To ensure we fit the difference between the centroid and physical forces, we use the flag ```--forces_key='delta_force'```. 

Since quantum nuclear effects are localized to short ranges, we use a relatively small cutoff using ```--r_max=3.0```. Since the regression task is easy, the model can be made simpler by reducing the value associated with the flag ```--hidden_irreps='64x0e'```.

> **Note:**
> This step assumes that your have access to ```mace``` python bindings and cli.  

## Prediction


### Step 4: Production simulations

Go to the directory ```3_production_simulations``` and checkout the ```input.xml``` file. Here, we run a simple molecular dynamics simulation at 300 K using a potential that adds a correction to the physical potential. The two potentials can be added by defining two forcefield sockets using

```xml
<ffsocket name='maceoff23' mode='unix' pbc='false'>
    <address> driver </address>
</ffsocket>
<ffsocket name='maceoff23-pigs' mode='unix' pbc='false'>
    <address> driver-pigs </address>
</ffsocket>
```

and setting up two force components using

```xml
<forces>
    <force forcefield='maceoff23' weight='1'> </force>
    <force forcefield='maceoff23-pigs' weight='1'> </force>
</forces>
```

To estimate the vibrational spectrum, we print out the positions every 2 fs using: 

```xml
<trajectory filename='pos' stride='4' flush='100' format='ase'> positions </trajectory>
```


The Te PIGS potential is run using an ASEClient setup in the ```run-ase-pigs.py``` file. 

The simulation can be run using the bash script. 

```bash
bash run.sh
```

> **Note:**
> This step assumes that you have a working installation of i-PI using `pip` or have sourced the `env.sh` file in the top-level source directory and your Python environment has access to ```mace```.


### Step 5: Predicting polarization and polarizability

To estimate IR and Raman spectra, we need the system's total polarization and polarizability over its dynamical trajectory. Go to the directory ```4_dielectric_response_prediction``` and run the ``` get_dielectric_response.sh ``` script. This step uses the MACE dipole and polarizability model developed in Ref. [4]. This step is optional if you only want to estimate the vibrational density states. This script produces the following files: 

| File Name             | Description                                                                 |
|-----------------------|-----------------------------------------------------------------------------|
| `MACE_mu.txt`         | contains the 3-dimensional polarization vector                              |
| `MACE_alpha.txt`      | contains the 9-dimensional flattened polarizability tensor                  |
| `MACE_alpha_sh.txt`   | contains the six (L=0,m=0 L=2,m=-2,-1,0,1,2) spherical harmonic components of the polarizability tensor |


> **Note:**
> The dipole-polarizability feature is not yet in the main branch of `mace`. You can use this fork branch for now ```https://github.com/venkatkapil24/mace/tree/feat/mu_alpha```. This demo will be updated when this feature is merged into the main repository. 

### Step 6: Predicting polarization and polarizability

Go to ```5_final_spectra``` and run the bash script `get_spectra.sh` to estimate the vibrational density of states, IR, isotropic Raman and anisotropic Raman spectra. The script produces the following files:

| File Name              | Description                   |
|------------------------|-------------------------------|
| `xx_der_facf.data`     | Vibrational density of states   |
| `mm_der_facf.data`     | IR   |
| `L0L0_der_facf.data`   | Isotropic Raman   |
| `L2L2_der_facf.data`   | Anisotropic Raman   |

The first column of these files is the frequency in atomic units while the second column is the intensity modulo a multiplicative constant. Multiply with 219474.63,  to convert the frequencies into 'cm-1' units. 

> **Note:**
> This step assumes that you have a working installation of i-PI using `pip` or have sourced the `env.sh` file in the top-level source directory.

## References 

1. Musil, F., Zaporozhets, I., Noé, F., Clementi, C., & Kapil, V. (2022). Quantum dynamics using path integral coarse-graining. *The Journal of Chemical Physics, 157*(18), 181102. [https://doi.org/10.1063/5.0120386](https://doi.org/10.1063/5.0120386)

2. Batatia, I., Kovacs, D. P., Simm, G. N. C., Ortner, C., & Csanyi, G. (2022). MACE: Higher Order Equivariant Message Passing Neural Networks for Fast and Accurate Force Fields. *Advances in Neural Information Processing Systems*. [https://openreview.net/forum?id=YPpSngE-ZU](https://openreview.net/forum?id=YPpSngE-ZU)

3. Kovács, D. P., Moore, J. H., Browning, N. J., Batatia, I., Horton, J. T., Kapil, V., Witt, W. C., Magdău, I.-B., Cole, D. J., & Csányi, G. (2023). Mace-off23: Transferable machine learning force fields for organic molecules. [https://doi.org/10.48550/ARXIV.2312.15211](https://doi.org/10.48550/ARXIV.2312.15211)

4. Kapil, V., Kovács, D. P., Csányi, G., & Michaelides, A. (2023). First-principles spectroscopy of aqueous interfaces using machine-learned electronic and quantum nuclear effects. *Faraday Discussions*. [https://doi.org/10.1039/D3FD00113J](https://doi.org/10.1039/D3FD00113J)


## Contact

Feel free to email Venkat Kapil [v.kapil@ucl.ac.uk](v.kapil@ucl.ac.uk) for any questions. 

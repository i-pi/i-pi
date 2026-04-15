# Permutationally Invariant Polynomial (PIP) Potential Example
Runs an example using a PIP potential of CH<sub>5</sub><sup>+</sup> for a simple gas-phase MD simulation at 100 K. Other potentials and corresponding drivers are also available and runs in a similar way. See detailed instructions in each separate potential.

## Installation
(1) Run the i-PI `env.sh` script to set up the environment variables for i-PI, specifically, the `IPI_ROOT` variable.

(2) Download the CH<sub>5</sub><sup>+</sup> potential and the associated Fortran driver from the repository pippes at https://github.com/szquchen/pippes. The potential and the inteface for CH<sub>5</sub><sup>+</sup> are located in the folder `03_ch5p`. This can be downloaded to a location separate from i-PI.

(3) Go to `03_ch5p/src/pes_shell.f90` and modify the variable `coef_path` to the actual path of the `coeff` directory on your machine.

(4) A Makefile is provided in `03_ch5p` to compile and install the client. It assumes gcc and gfortran compiler; change it if you use other compilers. Use command `make` to compile and a symbolic link to the executable will be automatically created in `IPI_ROOT/bin/`.

## Running the Example
Run the following to start the i-PI simulation:  
```bash
i-pi input.xml &
sleep 5  
i-pi-driver -u -h driver -m ch5p
```

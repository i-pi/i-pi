# Permutationally Invariant Polynomial (PIP) Potential Example
Runs an example using a PIP potential of CH5<sup>+</sup> for a simple gas-phase MD simulation at 100 K.

## Installation
(1) Run the i-PI `env.sh` script to set up the environment variables for i-PI, specifically, the `IPI_ROOT` variable.

(2) Download the CH5<sup>+</sup> potential and the associated interface from the repository pippes at https://github.com/szquchen/pippes. The potential and the inteface for CH5<sup>+</sup> are located in the folder `03_ch5p`.

(3) Go to `03_ch5p/src/pes_shell.f90` and modify the variable `coef_path` to the actual path of the `coeff` directory on your machine.

(4) A Makefile is provided in `03_ch5p` to compile and install the client. It assumes gcc and gfortran compiler; change it if you use other compilers. Simply use command `make` to compile and install.  

## Running the Example
Run the following to start the i-PI simulation:  
```bash
i-pi input.xml &
sleep 5  
i-pi-driver -u -h driver -m ch5p
```

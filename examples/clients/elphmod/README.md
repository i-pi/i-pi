elphmod - i-PI example
======================

This example connects i-PI with [elphmod](https://github.com/janberges/elphmod),
which provides energy surfaces for distorted systems based on *ab initio* data
for the symmetric one, using the approximation of a linearized electron-lattice
coupling. This realizes *model III* by Schobert *et al.*, [SciPost Phys. **16**,
046 (2024)](https://doi.org/10.21468/SciPostPhys.16.2.046).

The Python packages `elphmod` and `matplotlib` (for plotting only) are required.
Both can be installed with `pip` or from `conda-forge`.

This example performs the structural optimization of a charge-density wave in
monolayer tantalum disulfide. Running `make` performs the following steps:

1. Usually, the calculation is based on a Wannier Hamiltonian from Wannier90
   (`model_hr.dat`), force constants from the PHonon code of Quantum ESPRESSO
   (`model.ifc`), and the corresponding representation of the coupling from a
   [patched version](https://github.com/janberges/elphmod/tree/master/patches)
   of the EPW code (`model.epmatwp`, `model.wigner`). To avoid the *ab initio*
   step, here we use a nearest-neighbor model instead, which can be created with
   `make model` or `python3 model.py`.

2. Based on these data, we can build a supercell model and set up the driver
   with `make driver` or `python3 driver.py`. This will create files with the
   driver object (`driver.pickle`) and initial ionic positions (`driver.xyz`).
   Note that the dummy `driver.xyz` required for testing will be overwritten.
   Git will show the file as modified (or deleted after running `make clean`).

3. Now we are ready to run i-PI and the driver, which will communicate via a
   Unix domain socket, and perform the optimization with `make run`. This will
   start i-PI in the background, wait until it is ready, and start the driver
   with `i-pi-driver-py -u -m elphmod -o driver.pickle`. i-PI will output the
   potential energies (`run.out`), positions (`run.pos_0.xyz`), and forces
   (`run.for_0.xyz`) for all steps.

4. Finally, we can view the trajectory with `make plot` or `python3 plot.py`.
   This should open an interactive figure and finally create `plot.pdf`, where
   arrows indicate atomic displacements. Since we have used a simplistic model,
   the periodic lattice distortion differs from the expected 3-by-3 pattern.

Note that saving the driver to an intermediate file `driver.pickle` has two
advantages: First, when running it through `i-pi-driver-py`, only a single
parameter (the filename) has to be specified on the command line, rather than
all details of the flexible setup in `driver.py`. Second, setting up the driver
can be quite time consuming, so it is more efficient to do this only once when
launching it multiple times later.

It is alternatively possible to run the driver through the module `ipi._driver`,
which is however only available if i-PI has been installed using `setuptools`,
e.g., through `pip`. This can be tested with `make run2` or `python3 run.py`.

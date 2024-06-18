Kinetic Monte Carlo for vacancy diffusion in an Al-Si-Mg alloy
==============================================================

Author: `Michele Ceriotti <michele.ceriotti@gmail.com>`

This example is somewhat atypical, in that it demonstrates a very ad hoc implementation of 
KMC, designed specifically to run simulations of vacancy diffusion in an FCC alloy of aluminum
with ~1% Si and Mg. 
The example is interesting because it demonstrates several tricks that rely on the way i-PI treats
atoms in communicating with the driver code, and uses a machine-learning potential to compute
energy and forces. 
The implmentation was first used in 
[A. C. P. Jain, D. Marchand, A. Glensk, M. Ceriotti, and W. A. Curtin, "Machine learning for metallurgy III: A neural network potential for Al-Mg-Si," Phys. Rev. Materials 5(5), 053805 (2021)](https://doi.org/10.1103/PhysRevMaterials.5.053805), 
and could probably be generalized with relatively little effort, so if you need to perform similar
simulations for a different systems, read on. 

Quasi-on-lattice KMC
--------------------

In general terms, KMC involves enumerating the possible reaction pathways for a given configuration
of the system, estimating the associated rates, choosing one mechanism at random with the appropriate
probability, and evolving the configuration accordingly. The process is repeated, giving a coarse-grained
dynamics that can span much larger time scales than explicit molecular dynamics. 

The implementation used here is peculiar, in that it considers an ideal fcc lattice decorated with
atoms *and* vacancies, which define the configurations. However, the energies and barriers are estimated
based on the energy of *relaxed* structures obtained starting from the ideal-lattice configuration.
This approach allows incorporating the contributions from relaxation into the energetics of the KMC, while 
still being able to enumerate possible reaction pathways based on the simple topology of the lattice. 

Vacancies are assumed to be able to diffuse in all directions, and the barrier is taken to be a fixed 
shift on top of the mean between reactant and product energies. 
Thus, at each step the code has to estimate the energy of all the possible outcomes for vacancy 
diffusion, and compute their relaxed energies. 


Required software
-----------------

This simulation requires the use of [LAMMPS](https://www.lammps.org/), including the 
`fix_ipi`, as well as compiling the [N2P2](https://compphysvienna.github.io/n2p2) library
and its LAMMPS interface, to compute the Behler-Parrinello potential for Al-Mg-Si alloys. 

Usage
-----

After compiling LAMMPS with the N2P2 module, You will have to launch i-PI

```bash
i-pi input.xml &> log.ipi
```

and then run one or more instances of LAMMPS

```bash
$LAMMPS_EXECUTABLE < in.lmp &> log.lammps
```

What to look for
----------------

The `N2P2` parameters (and the network weights) are in the folder `potential_n2p2_v2ag`. 
The file `input.nn` contains all the symmetry-function parameters, as well as the parameters that 
were used for training.


The `motion` class that implements KMC is very verbose, and reflects the highly ad-hoc nature 
of the implementation. It is important to understand that i-PI in this case will essentially ignore
the initialization structure, and create a new one, with a `ncell×ncell×ncell` replication of a primitive 
fcc cell, with the stated (conventional cell) lattice parameter `a0`. The cell defaults to a pure Al
matrix, and the first atoms are set to be `nsi` Si atoms and `nmg` Mg atoms. The last `nvac` atoms are 
set to vacancies (V), and we will discuss later how these are treated. 

Then, the motion class specifies what geometry optimization algorithm should be run (these can be any 
of the `optimizer` classes uses with `motion mode='geop'`) and how many steps of relaxation (`nstep`) should be performed after each KMC reaction. It then defines the estimated diffusion barrier and prefactor.

Finally, the class implements some optimizations. First, `<neval>` potential reactive pathways are
evaluated in parallel (provided enough calculators are available). Second, reactions are cached, 
i.e. the relaxed positions and energies are saved and re-used in case the system visits again the same
configuration. This is useful when the system gets trapped into a loop of states connected by low 
reaction barriers. 

```xml
<motion mode='al-kmc'>
    <al6xxx_kmc>
       <ncell> 5 </ncell>
       <a0 units='angstrom'> 4.057 </a0>       
       <nsi> 6 </nsi>
       <nmg> 6 </nmg>
       <nvac> 2 </nvac>
       <optimizer mode='lbfgs'> 
            <ls_options> <iter> 3 </iter> </ls_options>
       </optimizer>
       <nstep> 5 </nstep>
       <diffusion_barrier_al units='electronvolt'> 0.52  </diffusion_barrier_al> <!-- Mantina2009 dHm -->
       <diffusion_prefactor_al units='terahertz'> 16.6 </diffusion_prefactor_al> <!-- Mantina2009 v*  -->
       <neval> 8 </neval>
       <ecache_file> KMC_ECACHE </ecache_file>
       <qcache_file> KMC_QCACHE </qcache_file>
       <max_cache_len> 100 </max_cache_len>   
   </al6xxx_kmc>
</motion>
```

This application is also a good example of how the communication between i-PI and the driver code can
be exploited to realize complicated simulations. Here, we want to define some "virtual atoms" that 
act as placeholders for vacancies, but LAMMPS (at least with the N2P2 plugin) has no mechanism to 
interpret these as anything but actual atoms. We get around this problem by using the 
`activelist` option to the forcefield command, that specifies the list of atoms that are
actually passed on to the driver for calculation. So, even if the fcc box has 125 atoms,
and on the i-PI side the last to 'V' atoms get printed out and manipulated, only the first
123 "real" atoms are passed to LAMMPS. 

```xml
<activelist> [
       0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,
33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,
65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,
97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122] </activelist>
```

Incidentally, this means that on the LAMMPS side we should have only 123 atoms: 
6 Si, 6 Mg and then 111 Al atoms. 

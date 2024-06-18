Alchemical isotope exchanges
============================

Author: `Michele Ceriotti <michele.ceriotti@gmail.com>`

This example demonstrates the use of the alchemical exchange method implemented in i-PI
(cf. the theory introduced in [Jian Liu, Richard S Andino, Christina M Miller, Xin Chen, David M Wilkins, Michele Ceriotti, David E Manolopoulos, “A surface-specific isotope effect in mixtures of light and heavy water”, J. Phys. Chem. C 117(6), 2944-2951 (2013)](http://dx.doi.org/10.1021/jp311986m) [bibtex](https://www.doi2bib.org/bib/10.1021/jp311986m), 
and the implementation discussed in [Cheng, Bingqing, Jörg Behler, Michele Ceriotti, “Nuclear Quantum Effects in Water at the Triple Point: Using Theory as a Link Between Experiments.” J. Phys. Chem. Lett. 7(12), 2210-2215 (2016)](https://ipi-code.org/about/features/dx.doi.org/10.1021/acs.jpclett.6b00729) [bibtex](http://www.doi2bib.org/bib/10.1021%2Facs.jpclett.6b00729) ).

The masses of different isotopes (and the corresponding labels) are exchanged with a 
Monte Carlo acceptance criterion based on the full path integral Hamiltonian, leading to
equilibration between different molecular positions and/or thermodynamic phases. 

In this example, the simulation starts with 32 water HOD molecules, and rapidly reaches
an equilibrium with 1:2:1 H2O:HOD:D2O ratio; simulations at a lower temperature (and 
using a higher number of PI replicas) will lead to a small deviation from this ratio,
due to the presence of equilibrium isotope effects. 


Required software
-----------------

This simulation requires the use of [LAMMPS](https://www.lammps.org/), including the 
`fix_ipi` module and all the modules needed to use a TIP4P water potential.


Usage
-----

Much as with any PIMD simulation, you will have to launch i-PI

```bash
i-pi input.xml &> log.ipi
```

and then run one or more instances of the driver code

```bash
$LAMMPS_EXECUTABLE < in.lmp &> log.lammps
```

Note that the simulation uses some common files from the `init_files` folder,
so you may have to copy those and update paths if you want to run it outside
the examples folder.


What to look for
----------------

The `input.xml` file describes a run-of-the-mill PIMD calculation. The key section for
this example is 

```xml
<motion mode='alchemy'>
  <alchemy>
    <names> [ H, D ] </names>
    <nxc> 10 </nxc>
  </alchemy>
</motion>
```

That requires attempting on average 10 exhanges between H and D atoms per step. Given
that the exchange is assumed to only change the mass matrix, the acceptance probability
can be computed quickly, without having to re-evaluate the potential energy. 

On the other hand, note that this assumes that the atoms being swapped are bona fide 
isotopes, and that they behave chemically in the same way: in this example, the driver
code does not need to know about the swap having happened, because H and D behave the same 
from the point of view of the forcefield.

If you are looking for actual chemical swaps, these can be realized using a different
class, [atomswap](https://ipi-code.org/i-pi/input-reference.html#atomswap).

See in the output positions file how the atomic labels are swapped between structures.

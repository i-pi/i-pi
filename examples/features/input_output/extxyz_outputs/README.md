Extended XYZ formatted outputs
=================================================

The example is for a water molecule. It can be run with the pswater driver mode. The positions, forces, and the kinetic energy tensor are printed in standard and extended XYZ formatted file. When non position quantities are printed out, the position arrays will correspond to zero. The extxyz files can be combined into a single one using:

```
from ase.io import iread, write
  
for atoms, atoms_f, atoms_k in zip(iread('simulation.pos-extxyz_0.extxyz'), iread('simulation.frc-extxyz_0.extxyz'), iread('simulation.kin-extxyz.extxyz')):
    for key, value in atoms_f.arrays.items():
        if key not in atoms.arrays:
            atoms.arrays[key] = value

    for key, value in atoms_k.arrays.items():
        if key not in atoms.arrays:
            atoms.arrays[key] = value

    write('-', atoms, format='extxyz')
```



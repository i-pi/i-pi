# Syntax on

set term qt 0 position 100,200
set key l t

p 'energies.bfgs.dat' w lp lw 2, \
  'energies.bfgstrm.dat' w lp lw 2, \
  'energies.cg.dat' w lp lw 2, \
  'energies.lbfgs.dat' w lp lw 2, \
  'energies.sd.dat' w lp lw 2, \
  'ase_comparison/energies-ase-relaxed.bfgs.dat' w p ps 3

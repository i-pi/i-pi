# Syntax on

set term qt 1 position 800,200
set key c c
set title "'far-away' test"

p 'energies.bfgs.dat' w lp lw 2, \
  'energies.bfgstrm.dat' w lp lw 2, \
  'energies.cg.dat' w lp lw 2, \
  'energies.lbfgs.dat' w lp lw 2, \
  'energies.sd.dat' w lp lw 2

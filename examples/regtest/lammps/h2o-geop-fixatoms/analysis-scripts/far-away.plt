# Syntax on

set term qt 1 position 800,200
set key c c
set title "'far-away' test. All lines must be horizontal."

p 'energies.far-away.bfgs.dat' w lp lw 2, \
  'energies.far-away.bfgstrm.dat' w lp lw 2, \
  'energies.far-away.cg.dat' w lp lw 2, \
  'energies.far-away.lbfgs.dat' w lp lw 2, \
  'energies.far-away.sd.dat' w lp lw 2

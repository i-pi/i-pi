# Syntax on

set term qt 0 position 100,200
set key l t
set title "Main test. 0 and 1 should be equal, 2 above them.\nDifference between methods should be below threshold."

p 'energies.bfgs.dat' w lp lw 2, \
  'energies.bfgstrm.dat' w lp lw 2, \
  'energies.cg.dat' w lp lw 2, \
  'energies.lbfgs.dat' w lp lw 2, \
  'energies.sd.dat' w lp lw 2

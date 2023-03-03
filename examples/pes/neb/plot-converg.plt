#! gnuplot

set term qt font 'Sans, 12'
set key outside

set title "Zundel, 1 atom fixed"

p for [i=2:17] 'dbfgs.neb'  u 1:(column(i)) t sprintf("new %i", i-1)  w lp lw 2 dt 2

pause mouse close

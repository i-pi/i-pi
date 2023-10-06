#! gnuplot
set term qt 1 size 800,600 font 'Sans, 14'
unset key
set grid

set xlabel 'Bead (counted from 1)'
set ylabel 'Potential (Hartree)'

#set palette negative
unset colorbox

p 'dbfgs.neb' matrix u 1:3:2 every 1:10:1 w lp lw 2 palette

pause mouse close

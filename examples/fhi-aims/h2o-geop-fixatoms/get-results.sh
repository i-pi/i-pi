#!/bin/bash

mode=$1
if [ $mode == '' ]
then
    echo -e "ERROR: 1 argument needed:\n\tmode - bfgs, bfgstrm, cg or sd."
    exit
fi

rm energies-ipi-relaxed.${mode}.dat ase_comparison/energies-ase-relaxed.bfgs.dat
for i in `seq 0 3`
do
    grep 'Total energy uncorrected' aims.${i}-fixed.${mode}.out| tail -n 1| awk "{print $i,\$6}" >> energies-ipi-relaxed.${mode}.dat
    grep -o 'energy=-[.0-9]*' ase_comparison/relaxation_${i}.xyz|tail -n 1 |grep -o [-.0-9]* >> ase_comparison/energies-ase-relaxed.bfgs.dat
done

gnuplot -p << EOF
p 'energies-ipi-relaxed.${mode}.dat' w lp lw 2, 'ase_comparison/energies-ase-relaxed.bfgs.dat' w lp ps 2
EOF

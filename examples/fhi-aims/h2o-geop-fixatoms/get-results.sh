#!/bin/bash

# Get numbers for all optimizers:
if [ $# != 1 ]
then
    echo "Performing all: bfgs, bfgstrm, lbfgs, cg, sd."

    rm ase_comparison/energies-ase-relaxed.bfgs.dat
    for i in `seq 0 2`
    do
        grep -o 'energy=-[.0-9]*' ase_comparison/relaxation_${i}.xyz | tail -n 1 | sed 's/energy=//g' >> ase_comparison/energies-ase-relaxed.bfgs.dat
    done

    for mode in bfgs lbfgs bfgstrm cg sd
    do
        rm energies.${mode}.dat
        for i in `seq 0 2`
        do
            grep 'Total energy uncorrected' aims.${i}-fixed.${mode}.out | tail -n 1 | awk "{print $i,\$6}" >> energies.${mode}.dat
        done
    done

    gnuplot -p all.plt
    exit

# Get numbers for a particular optimizer:
elif [$# == 1]; then mode=$1
else 
    echo -e "ERROR: 1 argument needed:\n\tmode - bfgs, bfgstrm, lbfgs, cg or sd."
    exit
fi

rm energies-ipi-relaxed.${mode}.dat ase_comparison/energies-ase-relaxed.bfgs.dat
for i in `seq 0 2`
do
    grep 'Total energy uncorrected' aims.${i}-fixed.${mode}.out | tail -n 1 | awk "{print $i,\$6}" >> energies-ipi-relaxed.${mode}.dat
    grep -o 'energy=-[.0-9]*' ase_comparison/relaxation_${i}.xyz | tail -n 1 | grep -o [-.0-9]* >> ase_comparison/energies-ase-relaxed.bfgs.dat
done

gnuplot -p all.plt
#gnuplot -p << EOF
#p 'energies-ipi-relaxed.${mode}.dat' w lp lw 2, \
#  'ase_comparison/energies-ase-relaxed.bfgs.dat' w lp ps 2
#EOF

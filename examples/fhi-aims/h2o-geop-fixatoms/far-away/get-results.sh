#!/bin/bash

if [ $# != 1 ]
then
    echo "Performing all: bfgs, bfgstrm, lbfgs, cg, sd."
    for mode in bfgs lbfgs bfgstrm cg sd
    do
        rm energies.${mode}.dat
        for i in `seq 0 1`
        do
            grep 'Total energy uncorrected' aims.${i}-fixed.${mode}.out| tail -n 1| awk "{print $i,\$6}" >> energies.${mode}.dat
        done
    done

    gnuplot -p all.plt
    exit
elif [$# != 1]; then mode=$1
else 
    echo -e "ERROR: 1 argument needed:\n\tmode - bfgs, bfgstrm, lbfgs, cg or sd."
    exit
fi

rm energies.${mode}.dat
for i in `seq 0 1`
do
    grep 'Total energy uncorrected' aims.${i}-fixed.${mode}.out| tail -n 1| awk "{print $i,\$6}" >> energies.${mode}.dat
done

gnuplot -p all.plt
wait

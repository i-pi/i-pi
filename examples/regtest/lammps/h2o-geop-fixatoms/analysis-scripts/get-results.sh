#!/bin/bash

echo "Collecting all: bfgs, bfgstrm, lbfgs, cg, sd."
for mode in bfgs lbfgs bfgstrm cg sd
do
    echo $mode
    if [ -e energies.${mode}.dat ]; then rm energies.${mode}.dat; fi
    for i in `seq 0 2`
    do
        tail -n 1 regtest-run/${mode}-${i}fixed/sim-${i}fixed.${mode}.out | awk "{print $i,\$2}" >> energies.${mode}.dat
    done
    cat energies.${mode}.dat

    # far-away:
    echo -e "\t\t\t'far-away' test:"
    if [ -e energies.far-away.${mode}.dat ]; then rm energies.far-away.${mode}.dat; fi
    for i in `seq 0 1`
    do
        tail -n 1 regtest-run/${mode}-far-away-${i}fixed/sim-far-away-${i}fixed.${mode}.out | awk "{print $i,\$2}" >> energies.far-away.${mode}.dat
    done
    cat energies.far-away.${mode}.dat | sed 's/^/\t\t\t/'
    echo "-----------------------------------------"
done

gnuplot -p analysis-scripts/012fixed.plt
gnuplot -p analysis-scripts/far-away.plt

exit

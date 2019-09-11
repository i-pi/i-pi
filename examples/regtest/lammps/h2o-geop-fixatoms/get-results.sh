#!/bin/bash

if [ $# != 1 ]
then
    echo "Performing all: bfgs, bfgstrm, lbfgs, cg, sd."
    for mode in bfgs lbfgs bfgstrm cg sd
    do
        rm energies.${mode}.dat
        for i in `seq 0 2`
        do
            tail -n 1 sim-${i}-fixed.${mode}.out | awk "{print $i,\$2}" >> energies.${mode}.dat
        done
    done

    gnuplot -p all.plt

    cd far-away
        ./get-results.sh
    cd ..

    exit
elif [$# != 1]; then mode=$1
else 
    echo -e "ERROR: 1 argument needed:\n\tmode - bfgs, bfgstrm, lbfgs, cg or sd."
    exit
fi

rm energies.${mode}.dat
for i in `seq 0 2`
do
    tail -n 1 sim-${i}-fixed.${mode}.out | awk "{print $i,\$2}" >> energies.${mode}.dat
done

gnuplot -p all.plt

cd far-away
    ./get-results.sh ${mode}
cd ..

#!/bin/bash

if [ $# != 1 ]
then
#    echo "Total number of force calls reported by i-PI:"
    for mode in bfgs bfgstrm lbfgs cg sd
    do
        echo $mode
        echo "What i-pi says:"
        for file in `ls -1 log.*.${mode}.ipi`
        do
            grep 'Number of force calls' $file | grep -o '[0-9]*$' | awk "BEGIN{s=0} {s+=\$1}END{print \"\t$file     \", s}"
        done
        echo "What aims says (grep 'Begin self-consistency loop'):"
        for file in `ls -1 aims.*.${mode}.out`
        do
            grep 'Begin self-consistency loop' $file | wc -l | awk "{print \"\t$file    \", \$1 }"
        done
        echo
    done
    exit

elif [ $# == 1 ]; then mode=$1
else 
    echo -e "ERROR: 1 argument needed:\n\tmode - bfgs, bfgstrm, lbfgs, cg or sd."
    exit
fi

echo "Total number of force calls for '${mode}':"
echo "What i-pi says:"
for file in `ls -1 log.*.${mode}.ipi`
do
    grep 'Number of force calls' $file | grep -o '[0-9]*$' | awk "{s+=\$1}END{print \"$file    \", s}"
done
echo "What aims says (grep 'Begin self-consistency loop'):"
for file in `ls -1 aims.*.${mode}.out`
do
    grep 'Begin self-consistency loop' $file | wc -l | awk "{print \"$file    \", \$1 }"
done

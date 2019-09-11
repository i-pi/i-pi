#!/bin/bash

if [ $# != 1 ]
then
    for mode in bfgs bfgstrm lbfgs cg sd
    do
        echo $mode
        echo "What i-pi says:"
        for file in `ls -1 log.*.${mode}.ipi`
        do
            grep 'Number of force calls' $file | grep -o '[0-9]*$' | awk "BEGIN{s=0} {s+=\$1}END{print \"\t$file     \", s}"
        done
        echo
    done
    exit


elif [$# == 1]; then mode=$1
else 
    echo -e "ERROR: 1 argument needed:\n\tmode - bfgs, bfgstrm, lbfgs, cg or sd."
    exit
fi

echo "Total number of force calls for '${mode}':"
echo "What i-pi says:"
grep 'Number of force calls' log.*${mode}* | grep -o '[0-9]*$' | awk '{s+=$1}END{print s}'

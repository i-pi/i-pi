#!/bin/bash

if [ $# -ne 1 ]
then
    echo -e "ERROR: 1 argument needed:\n\tmode - bfgs, bfgstrm, cg or sd."
    exit
fi
mode=$1

echo "Total number of force calls for '${mode}':"
echo "What i-pi says:"
grep 'Number of force calls' log.*${mode}.ipi | grep -o '[0-9]*$' | awk '{s+=$1}END{print s}'
echo "What aims says (number of [re]initializations):"
grep 'Begin self-consistency loop' aims.*${mode}.out|wc -l

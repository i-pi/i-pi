#!/bin/bash

if [ $# -ne 1 ]
then
    echo -e "ERROR: 1 argument needed:\n\tmode - bfgs, bfgstrm, cg or sd."
    exit
fi
mode=$1

echo "Total number of force calls for '${mode}':"
echo "What i-pi says:"
for file in `ls -1 log.*${mode}.faraway.ipi`
do
    echo -e "$file  "
    grep 'Number of force calls' $file | grep -o '[0-9]*$' | awk '{s+=$1}END{print s}'
done

echo
echo "What aims says (grep 'Begin self-consistency loop'):"
for file in `ls -1 aims.*${mode}.faraway.out`
do
    echo -e "$file  "
    grep 'Begin self-consistency loop' $file | wc -l
done

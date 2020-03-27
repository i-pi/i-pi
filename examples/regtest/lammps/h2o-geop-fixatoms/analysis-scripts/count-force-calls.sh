#!/bin/bash

echo "Counting force calls"
for mode in bfgs lbfgs bfgstrm cg sd
do
    echo $mode
#    echo "What i-pi says:"
    for file in `ls -1 regtest-run/${mode}-?fixed/ipi_output.out 2>/dev/null`
    do
        grep 'Number of force calls' $file | grep -o '[0-9]*$' | awk "BEGIN{s=0} {s+=\$1}END{print \"\t$file     \", s}"
    done
    echo
done
exit

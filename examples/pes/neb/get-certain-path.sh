#!/bin/bash
# You can choose which path to extract:
# By defaut, the last path is taken, but if you give an integer parameter,
# corresponding step will be extracted.

natoms=`head -n 1 dbfgs.pos_00.xyz`

if [ $# == 1 ]   # by step index
then
    indx=$1

    if [ -f path.${indx}.xyz ]
    then
        rm path.${indx}.xyz
    fi
    for f in dbfgs.pos_*.xyz
    do
        echo "$natoms" >> path.${indx}.xyz
        grep -A $natoms "Step:\s*$indx\>" $f >> path.${indx}.xyz
    done
    sed -i -e 's/^\s\+/ /g' -e 's/:\s\+/: /g' path.${indx}.xyz
    echo "path.${indx}.xyz written."

elif [ $# == 0 ]  # the last step
then
    if [ -f path.xyz ]
    then
        rm path.xyz
    fi
    tail -q -n $((natoms + 2)) dbfgs.pos_*.xyz >> path.xyz
    sed -i -e 's/^\s\+/ /g' -e 's/:\s\+/: /g' path.xyz
    echo "path.xyz written."

else
    echo "Error: 0 or 1 argument (step) expected."
fi

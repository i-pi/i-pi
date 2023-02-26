#!/bin/bash
# A script to extract certain NEB path from the output of beads' positions.
# Made for NEB, but can be useful to analize other bead-based methods, e.g. instantons.
# Arguments:
# 1. base of the name of a files with beads coordinates, e.g. sim.pos_ for sim.pos_*.xyz files
#    CAREFUL: no consistency checks are done -
#    the files are expected to be correct and to have equal length.
# 2. (optional) a step index to extract that particular step, or nothing to extract the last step.

if [ $# == 2 ]   # extract a given step
then
    fnamebase="$1"
    indx="$2"
    natoms=`head -n 1 ${fnamebase}00.xyz`

    if [ -f path.${indx}.xyz ]
    then
        rm path.${indx}.xyz
    fi
    for f in ${fnamebase}*.xyz
    do
        echo "$natoms" >> path.${indx}.xyz
        grep -A $natoms "Step:\s*$indx\>" $f >> path.${indx}.xyz
    done
    sed -i -e 's/^\s\+/ /g' -e 's/:\s\+/: /g' path.${indx}.xyz
    echo "path.${indx}.xyz written."

elif [ $# == 1 ]  # extract the last step
then
    echo "Extracting the last step..."
    fnamebase="$1"
    natoms=`head -n 1 ${fnamebase}00.xyz`
    echo "$natoms atoms"

    if [ -f path.xyz ]
    then
        rm path.xyz
    fi
    tail -q -n $((natoms + 2)) ${fnamebase}*.xyz >> path.xyz
    sed -i -e 's/^\s\+/ /g' -e 's/:\s\+/: /g' path.xyz
    echo "path.xyz written."

else
    echo "Error: 1 (filename base) or 2 arguments (filename base, step) are expected."
    exit
fi

echo 'Done.'

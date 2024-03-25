#!/bin/bash

export OMP_NUM_THREADS=1

PATH_TO_PLUMED=/home/username/plumed2

source PATH_TO_PLUMED/sourceme.sh
rm /tmp/ipi*

if test -f "simulation.restart";
then
        echo "continue last simulation"
        i-pi simulation.restart >> log.i-pi &
else
        i-pi input-pin.xml > log.i-pi &
fi

sleep 10
srun --hint=nomultithread --exclusive --mem=16G -n 4  ~/source/lammps/src/lmp_intel_universal < in-1.lmp > log.lmp  &
srun --hint=nomultithread --exclusive --mem=16G -n 4  ~/source/lammps/src/lmp_intel_universal < in-2.lmp > log.lmp  &
srun --hint=nomultithread --exclusive --mem=16G -n 4  ~/source/lammps/src/lmp_intel_universal < in-3.lmp > log.lmp  &
srun --hint=nomultithread --exclusive --mem=16G -n 4  ~/source/lammps/src/lmp_intel_universal < in-4.lmp > log.lmp  &

wait 
exit 0
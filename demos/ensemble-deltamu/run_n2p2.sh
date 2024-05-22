#!/bin/bash

export OMP_NUM_THREADS=1

PATH_TO_PLUMED=/home/username/plumed2
source PATH_TO_PLUMED/sourceme.sh
rm /tmp/ipi*

LMP=~/source/lammps/src/lmp_intel_universal

if test -f "pinning-n2p2.restart";
then
        echo "continue last simulation"
        i-pi pinning-n2p2.restart >> log.i-pi &
else
        i-pi input_n2p2.xml > log.i-pi &
fi

sleep 10
srun --hint=nomultithread --exclusive --mem=16G -n 4  $LMP < in-1.lmp > log.lmp  &
srun --hint=nomultithread --exclusive --mem=16G -n 4  $LMP < in-2.lmp > log.lmp  &
srun --hint=nomultithread --exclusive --mem=16G -n 4  $LMP < in-3.lmp > log.lmp  &
srun --hint=nomultithread --exclusive --mem=16G -n 4  $LMP < in-4.lmp > log.lmp  &

wait 
exit 0

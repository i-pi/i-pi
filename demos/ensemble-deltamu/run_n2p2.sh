#!/bin/bash

source ./miniconda3/bin/activate n2p2_install

export CONDA_BASE="$PWD/miniconda3/bin"

ipi=i-pi
driver=i-pi-py_driver

plumed_path=$(which plumed)
echo "Plumed is located at: $plumed_path"

sleep_time=2

rm -f /tmp/ipi_nnp*

if test -f "pinning.restart";
then
        echo "continue last simulation"
        ${ipi} pinning.restart >> log.i-pi &
else
        ${ipi} input_n2p2.xml > log.i-pi & 
fi

sleep 10
# Run LAMMPS simulations
lmp_serial < in-1.lmp > log.lmp1 &
lmp_serial < in-2.lmp > log.lmp2 &
lmp_serial < in-3.lmp > log.lmp3 &
lmp_serial < in-4.lmp > log.lmp4 &
wait # hold on until simulations are finished

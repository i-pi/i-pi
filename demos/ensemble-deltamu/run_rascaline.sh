#!/bin/bash

source ./miniconda3/bin/activate rascaline_install

export CONDA_BASE="$PWD/miniconda3/bin"

export OMP_NUM_THREADS=1

ipi=i-pi
driver=i-pi-py_driver

rm /tmp/ipi*
plumed_path=$(which plumed)
echo "Plumed is located at: $plumed_path"

export MLIP_DIR=$PWD/H2O/driver/

sleep_time=2

if test -f "simulation.restart";
then
        echo "continue last simulation"
        ${ipi} simulation.restart >> log.i-pi &
else
        ${ipi} input_rascaline.xml > log.i-pi & 
fi

echo "# i-PI is running"
echo "# Waiting for ${sleep_time} (s) before executing driver"
sleep ${sleep_time}

${driver} -m lightning -a pinning -u -o ./rascaline_potential/example_5_converted.ckpt,template.xyz &
wait # ... until it's all done 
